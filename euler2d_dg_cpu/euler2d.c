#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define ADIABATIC_GAMMA (5.0 / 3.0)
#define PI 3.14159265359
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

// Legendre polynomials scaled by sqrt(2n+1), and their derivatives

real p0(const real xsi)
{
    return 1.0;
}

real p0_prime(const real xsi)
{
    return 0.0;
}

real p1(const real xsi)
{
    return sqrt(3.0) * xsi;
}

real p1_prime(const real xsi)
{
    return sqrt(3.0);
}

real p2(const real xsi)
{
    return sqrt(5.0) * 0.5 * (3.0 * xsi * xsi - 1.0);
}

real p2_prime(const real xsi)
{
    return sqrt(5.0) * 0.5 * (6.0 * xsi);
}

real p3(const real xsi)
{
    return sqrt(7.0) * 0.5 * (5.0 * xsi * xsi * xsi - 3.0 * xsi);
}

real p3_prime(const real xsi)
{
    return sqrt(7.0) * 0.5 * (15.0 * xsi * xsi - 3.0);
}

struct Cells set_cell(void)
{
    struct Cells cell;

    real nhat;

    #if (DG_ORDER == 1)
        real xsi[NFACE] = {0.0};
        real gw[NFACE]  = {2.0}; // 1D Gaussian weight       
    #elif (DG_ORDER == 2)
        real xsi[NFACE] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
        real gw[NFACE]  = {1.0, 1.0}; // 1D Gaussian weight
    #elif (DG_ORDER == 3)
        real xsi[NFACE] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
        real gw[NFACE]  = {5.0/9.0, 8.0/9.0, 5.0/9.0}; // 1D Gaussian weight
    #endif
    // cell nodes

    int nc = 0;

    for (int i = 0; i < NFACE; ++i)
    {
        for (int j = 0; j < NFACE; ++j)
        {   
            // nk = 0
            cell.node[nc].phi[0]    = p0(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidx[0] = p0_prime(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidy[0] = p0(xsi[i]) * p0_prime(xsi[j]);
            #if (DG_ORDER >= 2)
            // nk = 1
            cell.node[nc].phi[1]    = p0(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidx[1] = p0_prime(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidy[1] = p0(xsi[i]) * p1_prime(xsi[j]);
            // nk = 2
            cell.node[nc].phi[2]    = p1(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidx[2] = p1_prime(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidy[2] = p1(xsi[i]) * p0_prime(xsi[j]);
            #endif
            #if (DG_ORDER >= 3)
            // nk = 3
            cell.node[nc].phi[3]    = p0(xsi[i]) * p2(xsi[j]);
            cell.node[nc].dphidx[3] = p0_prime(xsi[i]) * p2(xsi[j]);
            cell.node[nc].dphidy[3] = p0(xsi[i]) * p2_prime(xsi[j]);
            // nk = 4
            cell.node[nc].phi[4]    = p1(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidx[4] = p1_prime(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidy[4] = p1(xsi[i]) * p1_prime(xsi[j]); 
            // nk = 5
            cell.node[nc].phi[5]    = p2(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidx[5] = p2_prime(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidy[5] = p2(xsi[i]) * p0_prime(xsi[j]);                       
            #endif
            cell.node[nc].gw = gw[i] * gw[j];   // 2D Gaussian weight

            nc = nc + 1;
        }
    }

    // face nodes

    for (int n = 0; n < NFACE; ++n)
    {
        // left face (x = -1)
        // nk = 0
        cell.faceli[n].phi[0]    = p0(-1.0) * p0(xsi[n]);
        cell.faceli[n].dphidx[0] = p0_prime(-1.0) * p0(xsi[n]);
        cell.faceli[n].dphidy[0] = p0(-1.0) * p0_prime(xsi[n]);
        #if (DG_ORDER >= 2)
        // nk = 1
        cell.faceli[n].phi[1]    = p0(-1.0) * p1(xsi[n]);
        cell.faceli[n].dphidx[1] = p0_prime(-1.0) * p1(xsi[n]);
        cell.faceli[n].dphidy[1] = p0(-1.0) * p1_prime(xsi[n]);
        // nk = 2
        cell.faceli[n].phi[2]    = p1(-1.0) * p0(xsi[n]);
        cell.faceli[n].dphidx[2] = p1_prime(-1.0) * p0(xsi[n]);
        cell.faceli[n].dphidy[2] = p1(-1.0) * p0_prime(xsi[n]);
        #endif
        #if (DG_ORDER >= 3)
        // nk = 3
        cell.node[nc].phi[3]    = p0(-1.0) * p2(xsi[n]);
        cell.node[nc].dphidx[3] = p0_prime(-1.0) * p2(xsi[n]);
        cell.node[nc].dphidy[3] = p0(-1.0) * p2_prime(xsi[n]);
        // nk = 4
        cell.node[nc].phi[4]    = p1(-1.0) * p1(xsi[n]);
        cell.node[nc].dphidx[4] = p1_prime(-1.0) * p1(xsi[n]);
        cell.node[nc].dphidy[4] = p1(-1.0) * p1_prime(xsi[n]); 
        // nk = 5
        cell.node[nc].phi[5]    = p2(-1.0) * p0(xsi[n]);
        cell.node[nc].dphidx[5] = p2_prime(-1.0) * p0(xsi[n]);
        cell.node[nc].dphidy[5] = p2(-1.0) * p0_prime(xsi[n]);                       
        #endif
        cell.faceli[n].gw = gw[n] * (-1.0); // 1D Gaussian weight * nhat

        // right face (x = +1)

        cell.faceri[n].phi[0]    = p0(1.0) * p0(xsi[n]);
        cell.faceri[n].dphidx[0] = p0_prime(1.0) * p0(xsi[n]);
        cell.faceri[n].dphidy[0] = p0(1.0) * p0_prime(xsi[n]);
        #if (DG_ORDER >= 2)
        cell.faceri[n].phi[1]    = p0(1.0) * p1(xsi[n]);
        cell.faceri[n].dphidx[1] = p0_prime(1.0) * p1(xsi[n]);
        cell.faceri[n].dphidy[1] = p0(1.0) * p1_prime(xsi[n]);

        cell.faceri[n].phi[2]    = p1(1.0) * p0(xsi[n]);
        cell.faceri[n].dphidx[2] = p1_prime(1.0) * p0(xsi[n]);
        cell.faceri[n].dphidy[2] = p1(1.0) * p0_prime(xsi[n]);
        #endif
        #if (DG_ORDER >= 3)
        // nk = 3
        cell.node[nc].phi[3]    = p0(1.0) * p2(xsi[n]);
        cell.node[nc].dphidx[3] = p0_prime(1.0) * p2(xsi[n]);
        cell.node[nc].dphidy[3] = p0(1.0) * p2_prime(xsi[n]);
        // nk = 4
        cell.node[nc].phi[4]    = p1(1.0) * p1(xsi[n]);
        cell.node[nc].dphidx[4] = p1_prime(1.0) * p1(xsi[n]);
        cell.node[nc].dphidy[4] = p1(1.0) * p1_prime(xsi[n]); 
        // nk = 5
        cell.node[nc].phi[5]    = p2(1.0) * p0(xsi[n]);
        cell.node[nc].dphidx[5] = p2_prime(1.0) * p0(xsi[n]);
        cell.node[nc].dphidy[5] = p2(1.0) * p0_prime(xsi[n]);                       
        #endif

        cell.faceri[n].gw = gw[n] * (1.0); // 1D Gaussian weight * nhat

        // bottom face (y = -1)

        cell.facelj[n].phi[0]    = p0(xsi[n]) * p0(-1.0);
        cell.facelj[n].dphidx[0] = p0_prime(xsi[n]) * p0(-1.0);
        cell.facelj[n].dphidy[0] = p0(xsi[n]) * p0_prime(-1.0);
        #if (DG_ORDER >= 2)
        cell.facelj[n].phi[1]    = p0(xsi[n]) * p1(-1.0);
        cell.facelj[n].dphidx[1] = p0_prime(xsi[n]) * p1(-1.0);
        cell.facelj[n].dphidy[1] = p0(xsi[n]) * p1_prime(-1.0);

        cell.facelj[n].phi[2]    = p1(xsi[n]) * p0(-1.0);
        cell.facelj[n].dphidx[2] = p1_prime(xsi[n]) * p0(-1.0);
        cell.facelj[n].dphidy[2] = p1(xsi[n]) * p0_prime(-1.0);
        #endif
        #if (DG_ORDER >= 3)
        // nk = 3
        cell.node[nc].phi[3]    = p0(xsi[n]) * p2(-1.0);
        cell.node[nc].dphidx[3] = p0_prime(xsi[n]) * p2(-1.0);
        cell.node[nc].dphidy[3] = p0(xsi[n]) * p2_prime(-1.0);
        // nk = 4
        cell.node[nc].phi[4]    = p1(xsi[n]) * p1(-1.0);
        cell.node[nc].dphidx[4] = p1_prime(xsi[n]) * p1(-1.0);
        cell.node[nc].dphidy[4] = p1(xsi[n]) * p1_prime(-1.0); 
        // nk = 5
        cell.node[nc].phi[5]    = p2(xsi[n]) * p0(-1.0);
        cell.node[nc].dphidx[5] = p2_prime(xsi[n]) * p0(-1.0);
        cell.node[nc].dphidy[5] = p2(xsi[n]) * p0_prime(-1.0);                       
        #endif  

        cell.facelj[n].gw = gw[n] * (-1.0); // 1D Gaussian weight * nhat

        // top face (y = +1)

        cell.facerj[n].phi[0]    = p0(xsi[n]) * p0(1.0);
        cell.facerj[n].dphidx[0] = p0_prime(xsi[n]) * p0(1.0);
        cell.facerj[n].dphidy[0] = p0(xsi[n]) * p0_prime(1.0);
        #if (DG_ORDER >= 2)
        cell.facerj[n].phi[1]    = p0(xsi[n]) * p1(1.0);
        cell.facerj[n].dphidx[1] = p0_prime(xsi[n]) * p1(1.0);
        cell.facerj[n].dphidy[1] = p0(xsi[n]) * p1_prime(1.0);

        cell.facerj[n].phi[2]    = p1(xsi[n]) * p0(1.0);
        cell.facerj[n].dphidx[2] = p1_prime(xsi[n]) * p0(1.0);
        cell.facerj[n].dphidy[2] = p1(xsi[n]) * p0_prime(1.0);
        #endif
        #if (DG_ORDER >= 3)
        // nk = 3
        cell.node[nc].phi[3]    = p0(xsi[n]) * p2(1.0);
        cell.node[nc].dphidx[3] = p0_prime(xsi[n]) * p2(1.0);
        cell.node[nc].dphidy[3] = p0(xsi[n]) * p2_prime(1.0);
        // nk = 4
        cell.node[nc].phi[4]    = p1(xsi[n]) * p1(1.0);
        cell.node[nc].dphidx[4] = p1_prime(xsi[n]) * p1(1.0);
        cell.node[nc].dphidy[4] = p1(xsi[n]) * p1_prime(1.0); 
        // nk = 5
        cell.node[nc].phi[5]    = p2(xsi[n]) * p0(1.0);
        cell.node[nc].dphidx[5] = p2_prime(xsi[n]) * p0(1.0);
        cell.node[nc].dphidy[5] = p2(xsi[n]) * p0_prime(1.0);                       
        #endif

        cell.facerj[n].gw = gw[n] * (1.0); // 1D Gaussian weight * nhat
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
            real x = x0 + (i + 0.5) * dx;
            real y = y0 + (j + 0.5) * dy;
            real xmid = 0.5 * (x0 + x1);
            real ymid = 0.5 * (y0 + y1);
            real *wij = &weights[NCONS * NK * (i * nj + j)];
            real r2 = power(x - xmid, 2) + power(y - ymid, 2);

            /*real amag = 1e-2;

            wij[0 * NK] = 1.0 + amag * 1.0 * sin(2.0*PI*x);
            wij[1 * NK] = 1.0 + amag * (-1.0) * sin(2.0*PI*x);
            wij[2 * NK] = amag *  (1.0) * sin(2.0*PI*x);
            wij[3 * NK] = (1.0 / ADIABATIC_GAMMA) / (ADIABATIC_GAMMA -1.0) + amag * 1.5 * sin(2.0*PI*x); */
            //if (j==0) printf("%f %e %e\n",x,sin(2.0*PI*x),wij[1 * NK]);
            
            real prim[NCONS];
            real cons[NCONS];
            real flux[NCONS];
            real rho, vx, vy, pressure;

            // For |y| > 0.25, we set Vx = -0.5 and ρ = 1, for |y| ≤ 0.25, Vx = 0.5 and ρ = 2. 

            /*if ( y < 0.75 * (y1-y0) && y > 0.25*(y1-y0) )
            {
                vx = 0.5;
                rho = 2.0;
            }
            else
            {
                vx = -0.5;
                rho = 1.0;
            }
            vy = 0.01*sin(6*x);
            */
            /*rho      = 2.0 + 0.5*sin(2.0*PI*x);
            vx       = 3.0;
            vy       = 0.0;
            pressure = 2.0;

            prim[0] = rho;
            prim[1] = vx;
            prim[2] = vy;
            prim[3] = pressure;

            primitive_to_conserved(prim,cons);
            primitive_to_flux_vector(prim,flux,0);

            //printf("%f %f %f %f\n", cons[0],cons[1],cons[2],cons[3]);
            //printf("%f %f %f %f\n", flux[0],flux[1],flux[2],flux[3]);

            wij[0 * NK] = cons[0];
            wij[1 * NK] = cons[1];
            wij[2 * NK] = cons[2];
            wij[3 * NK] = cons[3];
*/
            //if (square_root(r2) < 0.125)
            if (x < xmid)    
            {
                wij[0*NK] = 1.0;
                wij[1*NK] = 0.0;
                wij[2*NK] = 0.0;
                wij[3*NK] = 1.0;
            }
            else
            {
                wij[0*NK] = 0.1;
                wij[1*NK] = 0.0;
                wij[2*NK] = 0.0;
                wij[3*NK] = 0.125;
            }
            /*
            if (y < ymid)    
            {
                wij[0*NK] = 1.0;
                wij[1*NK] = -0.1;
                wij[2*NK] = 0.02*sin(6*x);
                wij[3*NK] = 1.0;
            }
            else
            {
                wij[0*NK] = 1.0;
                wij[1*NK] = 0.1;
                wij[2*NK] = 0.02*sin(6*x);
                wij[3*NK] = 1.0;
            }*/
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

    real cons[NCONS];
    real prim[NCONS];
    real consm[NCONS];
    real consp[NCONS];
    real primm[NCONS];
    real primp[NCONS];

    real flux_x[NCONS];
    real flux_y[NCONS];
    
    real *delta_weights = (real*) malloc(ni * nj * NCONS * NK * sizeof(real));

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
            int il = i - 1;
            int ir = i + 1;
            int jl = j - 1;
            int jr = j + 1;

            // Outflow BC
            if (il == -1)
                il += 1;

            if (ir == ni)
                ir -= 1;

            if (jl == -1)
                jl += 1;

            if (jr == nj)
                jr -= 1;

 /*           //Periodic BC
            if (il == -1)
                il = ni-1;

            if (ir == ni)
                ir = 0;

            if (jl == -1)
                jl = nj-1;

            if (jr == nj)
                jr = 0; */

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
                
                riemann_hlle(primm, primp, flux_x, 0);
                //printf("%d left face flux : %f %f %f %f\n", i, flux_x[0],flux_x[1],flux_x[2],flux_x[3]);
                
                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux_x[q] * cell.faceli[n].phi[l] * cell.faceli[n].gw;  
                        //printf("i=%d q=%d l=%d left face: %f %f %f\n\n", i, q, l, flux_x[q] * cell.faceli[n].phi[l] * cell.faceli[n].gw, cell.faceli[n].phi[l], cell.faceli[n].gw);
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
                
                riemann_hlle(primm, primp, flux_x, 0);
                //printf("%d right face flux : %f %f %f %f\n\n", i, flux_x[0],flux_x[1],flux_x[2],flux_x[3]);
            
                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux_x[q] * cell.faceri[n].phi[l] * cell.faceri[n].gw;
                        //printf("i=%d q=%d l=%d right face: %f %f %f\n\n", i,q, l,flux_x[q] * cell.faceri[n].phi[l] * cell.faceri[n].gw , cell.faceri[n].phi[l], cell.faceri[n].gw);
                    }
                }
            }
/*
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
                
                riemann_hlle(primm, primp, flux_y, 1);

                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux_y[q] * cell.facelj[n].phi[l] * cell.facelj[n].gw;  
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
                
                riemann_hlle(primm, primp, flux_y, 1);

                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux_y[q] * cell.facerj[n].phi[l] * cell.facerj[n].gw;  
                    }
                }
            } 
*/
            // Volume term
/*
            for (int n = 0; n < NCELL; ++n)
            {
                for (int q = 0; q < NCONS; ++q)
                {
                    cons[q] = 0.0;

                    for (int l = 0; l < NK; ++l)
                    {
                        cons[q] += wij[NK * q + l] * cell.node[n].phi[l];
                        //printf("%d %d %d %f %f %f\n",n, q, l, cons[q],wij[NK * q + l],cell.node[n].phi[l]);        
                    }
                }

                conserved_to_primitive(cons, prim);

                primitive_to_flux_vector(prim, flux_x, 0);
                primitive_to_flux_vector(prim, flux_y, 1);
                //if(i==2) printf("%f %f %f %f\n", flux_x[0],flux_x[1],flux_x[2],flux_x[3]);

                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] += flux_x[q] * cell.node[n].dphidx[l] * cell.node[n].gw;
                        dwij[NK * q + l] += flux_y[q] * cell.node[n].dphidy[l] * cell.node[n].gw;   
                    }
                }

            }*/
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
                    wij[ NK * q + l] += 0.5 * dwij[NK * q + l] * dt / dx; 
                }
            }
        }

    }
    free(delta_weights);
}

int main()
{
    const int ni = 128;
    const int nj = 1;
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
    real dt = dx * 0.01;

    while (time < 0.001)
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
