#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define ADIABATIC_GAMMA (5.0 / 3.0)
#define min2(a, b) (a) < (b) ? (a) : (b)
#define max2(a, b) (a) > (b) ? (a) : (b)

#ifdef SINGLE
typedef float real;
#define square_root sqrtf
#define power powf
#else
typedef double real;
#define square_root sqrt
#define power pow
#endif

#define NDIM 2
#define DG_ORDER 1
#define NFACE 1
#define NCELL 1
#define NK    1  // number of basis polynomials
#define NCONS 4  // number of conserved variables

struct Cells cell;

struct Cells
{
    // Legendre polynomials scaled by sqrt(n^2 + 1)

    // n = 0
    
    /*
    real nhat;
    real leg0   = 1.0;
    real dxleg0 = 0.0; // x derivative of leg0
    real dyleg0 = 0.0; // y derivative of leg0
    */
    struct Nodes
    {
        real phi[NK];
        real dphidx[NK];
        real dphidy[NK];
        real gw;
    };

    struct Nodes node[NCELL];
    struct Nodes faceli[NFACE]; // left face nodes
    struct Nodes faceri[NFACE]; // right face nodes
    struct Nodes facelj[NFACE]; // bottom face nodes
    struct Nodes facerj[NFACE]; // top face nodes   

};

struct Cells set_cell(void)
{
    struct Cells cell;

    // Legendre polynomials scaled by sqrt(n^2 + 1)

    // n = 0
    real nhat;
    real leg0   = 1.0;
    real dleg0dx = 0.0; // x derivative of leg0
    real dleg0dy = 0.0; // y derivative of leg0

    // cell nodes

    int nc = 0;

    for (int i = 0; i < 1; ++i)
    {
        for (int j = 0; j < 1; ++j)
        {   
            cell.node[nc].phi[0]    = leg0;
            cell.node[nc].dphidx[0] = dleg0dx;
            cell.node[nc].dphidy[0] = dleg0dy;
            cell.node[nc].gw = 1.0;   // 2D Gaussian weight

            ++nc;
        }
    }

    // face nodes

    for (int n = 0; n < 1; ++n)
    {
        // left face
        cell.faceli[n].phi[0]    = leg0;
        cell.faceli[n].dphidx[0] = dleg0dx;
        cell.faceli[n].dphidy[0] = dleg0dy;
        cell.faceli[n].gw = 1.0 * (-1.0); // 1D Gaussian weight * nhat
        // right face
        cell.faceri[n].phi[0]    = leg0;
        cell.faceri[n].dphidx[0] = dleg0dx;
        cell.faceri[n].dphidy[0] = dleg0dy;
        cell.faceri[n].gw = 1.0 *  (1.0); // 1D Gaussian weight * nhat   
        // bottom face
        cell.facelj[n].phi[0]    = leg0;
        cell.facelj[n].dphidx[0] = dleg0dx;
        cell.facelj[n].dphidy[0] = dleg0dy;
        cell.facelj[n].gw = 1.0 * (-1.0); // 1D Gaussian weight * nhat
        // top face
        cell.facerj[n].phi[0]    = leg0;
        cell.facerj[n].dphidx[0] = dleg0dx;
        cell.facerj[n].dphidy[0] = dleg0dy;
        cell.facerj[n].gw = 1.0 *  (1.0); // 1D Gaussian weight * nhat 
    }
    return cell;
}

void conserved_to_primitive(const real *cons, real *prim)
{
    const real rho = cons[0];
    const real px = cons[1];
    const real py = cons[2];
    const real energy = cons[3];

    const real vx = px / rho;
    const real vy = py / rho;
    const real kinetic_energy = 0.5 * rho * (vx * vx + vy * vy);
    const real thermal_energy = energy - kinetic_energy;
    const real pressure = thermal_energy * (ADIABATIC_GAMMA - 1.0);

    prim[0] = rho;
    prim[1] = vx;
    prim[2] = vy;
    prim[3] = pressure;
}

void primitive_to_conserved(const real *prim, real *cons)
{
    const real rho = prim[0];
    const real vx = prim[1];
    const real vy = prim[2];
    const real pressure = prim[3];

    const real px = vx * rho;
    const real py = vy * rho;
    const real kinetic_energy = 0.5 * rho * (vx * vx + vy * vy);
    const real thermal_energy = pressure / (ADIABATIC_GAMMA - 1.0);

    cons[0] = rho;
    cons[1] = px;
    cons[2] = py;
    cons[3] = kinetic_energy + thermal_energy;
}

real primitive_to_velocity_component(const real *prim, int direction)
{
    switch (direction)
    {
        case 0: return prim[1];
        case 1: return prim[2];
        default: return 0.0;
    }
}

void primitive_to_flux_vector(const real *prim, real *flux, int direction)
{
    const real vn = primitive_to_velocity_component(prim, direction);
    const real pressure = prim[3];
    real cons[4];
    primitive_to_conserved(prim, cons);

    flux[0] = vn * cons[0];
    flux[1] = vn * cons[1] + pressure * (direction == 0);
    flux[2] = vn * cons[2] + pressure * (direction == 1);
    flux[3] = vn * cons[3] + pressure * vn;
}

real primitive_to_sound_speed_squared(const real *prim)
{
    const real rho = prim[0];
    const real pressure = prim[3];
    return ADIABATIC_GAMMA * pressure / rho;
}

void primitive_to_outer_wavespeeds(const real *prim, real *wavespeeds, int direction)
{
    const real cs = square_root(primitive_to_sound_speed_squared(prim));
    const real vn = primitive_to_velocity_component(prim, direction);
    wavespeeds[0] = vn - cs;
    wavespeeds[1] = vn + cs;
}

void riemann_hlle(const real *pl, const real *pr, real *flux, int direction)
{
    real ul[4];
    real ur[4];
    real fl[4];
    real fr[4];
    real al[2];
    real ar[2];

    primitive_to_conserved(pl, ul);
    primitive_to_conserved(pr, ur);
    primitive_to_flux_vector(pl, fl, direction);
    primitive_to_flux_vector(pr, fr, direction);
    primitive_to_outer_wavespeeds(pl, al, direction);
    primitive_to_outer_wavespeeds(pr, ar, direction);

    const real am = min2(0.0, min2(al[0], ar[0]));
    const real ap = max2(0.0, max2(al[1], ar[1]));

    for (int i = 0; i < 4; ++i)
    {
        flux[i] = (fl[i] * ap - fr[i] * am - (ul[i] - ur[i]) * ap * am) / (ap - am);
    }
}

void initial_weights(real *weights, int ni, int nj, real x0, real x1, real y0, real y1)
{
    real dx = (x1 - x0) / ni;
    real dy = (y1 - y0) / nj;

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
            real x = (i + 0.5) * dx;
            real y = (j + 0.5) * dy;
            real *wij = &weights[NCONS * NK * (i * nj + j)];
            real r2 = power(x - 0.5, 2) + power(y - 0.5, 2);

            /*
            cw[0 * NK] = 1.0;
            cw[1 * NK] = 0.0;
            cw[2 * NK] = 0.0;
            cw[3 * NK] = 1.0 * exp(-r2/0.01);
            */

            //if (square_root(r2) < 0.125)
            if (y < 0.5)    
            {
                wij[0] = 1.0;
                wij[1] = 0.0;
                wij[2] = 0.0;
                wij[3] = 1.0;
            }
            else
            {
                wij[0] = 0.1;
                wij[1] = 0.0;
                wij[2] = 0.0;
                wij[3] = 0.125;
            }
        }
    }
}


struct UpdateStruct
{
    int ni;
    int nj;
    real x0;
    real x1;
    real y0;
    real y1;

    real *weights;
};

struct UpdateStruct update_struct_new(int ni, int nj, real x0, real x1, real y0, real y1)
{
    struct UpdateStruct update;
    update.ni = ni;
    update.nj = nj;
    update.x0 = x0;
    update.x1 = x1;
    update.y0 = y0;
    update.y1 = y1;

    update.weights = (real*) malloc(ni * nj * NCONS * NK * sizeof(real));

    return update;
}

void update_struct_del(struct UpdateStruct update)
{
    free(update.weights);
}

void update_struct_set_weights(struct UpdateStruct update, const real *weights_host)
{
    int ni = update.ni;
    int nj = update.nj;
    int num_zones = ni * nj;

    memcpy(
        update.weights,
        weights_host,
        num_zones * NCONS * NK * sizeof(real)
    );
}

void update_struct_get_weights(struct UpdateStruct update, real *weights_host)
{
    int num_zones = update.ni * update.nj;
    memcpy(weights_host,
        update.weights,
        num_zones * NCONS * NK * sizeof(real)
    );
}


void update_struct_do_advance_weights(struct UpdateStruct update, real dt)
{
    int ni = update.ni;
    int nj = update.nj;
    const real dx = (update.x1 - update.x0) / update.ni;
    const real dy = (update.y1 - update.y0) / update.nj;

    real consm[NCONS];
    real consp[NCONS];
    real primm[NCONS];
    real primp[NCONS];
    real flux [NCONS];
    
    real *delta_weights = (real*) malloc(ni * nj * NCONS * NK * sizeof(real));
    real dw[ni][nj][NCONS][NK];

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
            int il = i - 1;
            int ir = i + 1;
            int jl = j - 1;
            int jr = j + 1;

            if (il == -1)
                il += 1;

            if (ir == ni)
                ir -= 1;

            if (jl == -1)
                jl += 1;

            if (jr == nj)
                jr -= 1;

            /* */ real *wij = &update.weights[NCONS * NK * (i  * nj + j )];

            const real *wli = &update.weights[NCONS * NK * (il * nj + j )];
            const real *wri = &update.weights[NCONS * NK * (ir * nj + j )];
            const real *wlj = &update.weights[NCONS * NK * (i  * nj + jl)];
            const real *wrj = &update.weights[NCONS * NK * (i  * nj + jr)];
   
            real *dwij = &delta_weights[NCONS * NK * (i  * nj + j )];

            for (int q = 0; q < NCONS; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    dwij[NK * q + l] = 0.0;
                    dw[i][j][q][l] = 0.0;
                }
            }

            // Left Face
            for (int n = 0; n < NFACE; ++n)
            {
                for (int q = 0; q < NCONS; ++q)
                {
                    consm[q] = 0.0;
                    consp[q] = 0.0;

                    for (int l = 0; l < NK; ++l)
                    {
                        consm[q] += wli[NK * q + l] * cell.faceri[n].phi[l]; // right face of zone i -1 
                        consp[q] += wij[NK * q + l] * cell.faceli[n].phi[l]; // left face of zone i                     
                    }
                }

                conserved_to_primitive(consm, primm);
                conserved_to_primitive(consp, primp);
                
                riemann_hlle(primm, primp, flux, 0);

                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux[q] * cell.faceli[n].phi[l] * cell.faceli[n].gw;  
                    }
                }
            }

            // Right Face
            for (int n = 0; n < NFACE; ++n)
            {
                for (int q = 0; q < NCONS; ++q)
                {
                    consm[q] = 0.0;
                    consp[q] = 0.0;

                    for (int l = 0; l < NK; ++l)
                    {
                        consm[q] += wij[NK * q + l] * cell.faceri[n].phi[l]; // right face of zone i  
                        consp[q] += wri[NK * q + l] * cell.faceli[n].phi[l]; // left face of zone i + 1                     
                    }
                }

                conserved_to_primitive(consm, primm);
                conserved_to_primitive(consp, primp);
                
                riemann_hlle(primm, primp, flux, 0);

                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux[q] * cell.faceri[n].phi[l] * cell.faceri[n].gw;
                    }
                }
            }

            // Bottom Face
            for (int n = 0; n < NFACE; ++n)
            {
                for (int q = 0; q < NCONS; ++q)
                {
                    consm[q] = 0.0;
                    consp[q] = 0.0;

                    for (int l = 0; l < NK; ++l)
                    {
                        consm[q] += wlj[NK * q + l] * cell.facerj[n].phi[l]; // top face of zone j-1 
                        consp[q] += wij[NK * q + l] * cell.facelj[n].phi[l]; // bottom face of zone j                     
                    }
                }

                conserved_to_primitive(consm, primm);
                conserved_to_primitive(consp, primp);
                
                riemann_hlle(primm, primp, flux, 0);

                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux[q] * cell.facelj[n].phi[l] * cell.facelj[n].gw;  
                    }
                }
            }

            // Top Face
            for (int n = 0; n < NFACE; ++n)
            {
                for (int q = 0; q < NCONS; ++q)
                {
                    consm[q] = 0.0;
                    consp[q] = 0.0;

                    for (int l = 0; l < NK; ++l)
                    {
                        consm[q] += wij[NK * q + l] * cell.facerj[n].phi[l]; // top face of zone j  
                        consp[q] += wrj[NK * q + l] * cell.facelj[n].phi[l]; // bottom face of zone j + 1                     
                    }
                }

                conserved_to_primitive(consm, primm);
                conserved_to_primitive(consp, primp);
                
                riemann_hlle(primm, primp, flux, 0);

                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux[q] * cell.facerj[n].phi[l] * cell.facerj[n].gw;  
                    }
                }
            } 
        }
    }

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
            real *wij = &update.weights[NCONS * NK * (i  * nj + j )];
            real *dwij = &delta_weights[NCONS * NK * (i  * nj + j )];

            for (int q = 0; q < NCONS; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    wij[ NK * q + l] += dwij[NK * q + l] * dt / dx; 
                }
            }
        }

    }
    free(delta_weights);
}

int main()
{
    const int ni = 128;
    const int nj = 128;
    const int fold = 1;
    const real x0 = 0.0;
    const real x1 = 1.0;
    const real y0 = 0.0;
    const real y1 = 1.0;
    const real dx = (x1 - x0) / ni;
    const real dy = (y1 - y0) / nj;

    cell = set_cell();

    real *weights_host = (real*) malloc(ni * nj * NCONS * NK * sizeof(real));
    struct UpdateStruct update = update_struct_new(ni, nj, x0, x1, y0, y1);

    initial_weights(weights_host, ni, nj, x0, x1, y0, y1);
    update_struct_set_weights(update, weights_host);

    int iteration = 0;
    real time = 0.0;
    real dt = min2(dx, dy) * 0.05;

    while (time < 0.1)
    {
        clock_t start = clock();

        for (int i = 0; i < fold; ++i)
        {
            update_struct_do_advance_weights(update, dt);

            time += dt;
            iteration += 1;
        }
        clock_t end = clock();

        real seconds = ((real) (end - start)) / CLOCKS_PER_SEC;
        real mzps = (ni * nj / 1e6) / seconds * fold;
        printf("[%d] t=%.3e Mzps=%.2f\n", iteration, time, mzps);
    }

    update_struct_get_weights(update, weights_host);
    update_struct_del(update);

    FILE* outfile = fopen("euler2d.dat", "w");

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
            real *wij = &weights_host[NCONS * NK * (i * nj + j)];
            real x = (i + 0.5) * dx;
            real y = (j + 0.5) * dy;
            fprintf(outfile, "%f %f %f %f %f %f\n", x, y, wij[0*NK], wij[1*NK], wij[2*NK], wij[3*NK]);
        }
    }
    fclose(outfile);
    free(weights_host);
    return 0;
}
