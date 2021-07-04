#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

//#define ADIABATIC_GAMMA (5.0 / 3.0) 
#define ADIABATIC_GAMMA (2.0) 
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

#define NDIM 1
#define DG_ORDER 2
#define NFACE 1
#define NNODE 4
#define NK    2  // number of basis polynomials
#define NCONS 4  // number of conserved variables

double node[NNODE][NK + NDIM * NK + 1];

void set_cell()
{
    real xnode[NNODE] = {-1.0, -1/square_root(3.0), 1.0/square_root(3.0), 1.0};

    for (int nn = 0; nn < NNODE; ++nn)
    {
        node[nn][0] = 1.0;                          // phi_0 
        node[nn][1] = square_root(3.0) * xnode[nn]; // phi_1 
        node[nn][2] = 0.0 * 1.0;                    // d(phi_0)/dx * Gaussian weight
        node[nn][3] = square_root(3.0) * 1.0;       // d(phi_1)/dx * Gaussian weight
                          
        if (nn == 0) {
            node[nn][4] = -1.0; // nhat at left face
        }
        else 
        {
            node[nn][4] = 1.0; // nhat at right face, also set to 1.0 for cell nodes
        }
        /*
        node[nn][0] = 1.0;                          // phi_0 
        node[nn][1] = 0.0 * 1.0;                    // d(phi_0)/dx * Gaussian weight
                          
        if (nn == 0) {
            node[nn][2] = -1.0; // nhat at left face
        }
        else 
        {
            node[nn][2] = 1.0; // nhat at right face, also set to 1.0 for cell nodes
        } */
    }
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
    real cons[NCONS];
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
    real ul[NCONS];
    real ur[NCONS];
    real fl[NCONS];
    real fr[NCONS];
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

    for (int i = 0; i < NCONS; ++i)
    {
        flux[i] = (fl[i] * ap - fr[i] * am - (ul[i] - ur[i]) * ap * am) / (ap - am);
    }
}

void initial_conserved(real *conserved, int ni, real x0, real x1)
{
    real dx = (x1 - x0) / ni;

    for (int i = 0; i < ni; ++i)
    {
        {
            real x = x0 + (i + 0.5) * dx;
            real xmid = (x1 + x0) / 2.0;

            real *cons = &conserved[NK * NCONS * (i)];
                
            cons[0*NK] = 1.0;
            cons[1*NK] = 0.0;
            cons[2*NK] = 0.0;
            cons[3*NK] = 1.0* exp(-(x-xmid)*(x-xmid)/0.01) / (ADIABATIC_GAMMA - 1.0);

            cons[0*NK+1] = 0.0;
            cons[1*NK+1] = 0.0;
            cons[2*NK+1] = 0.0;
            cons[3*NK+1] = 0.0;

            /*if ( x < 0.5)
            {
                cons[0*NK] = 1.0;
                cons[1*NK] = 0.0;
                cons[2*NK] = 0.0;
                cons[3*NK] = 1.0 / (ADIABATIC_GAMMA - 1.0);

                cons[0*NK+1] = 0.0;
                cons[1*NK+1] = 0.0;
                cons[2*NK+1] = 0.0;
                cons[3*NK+1] = 0.0;
            }
            else
            {
                cons[0*NK] = 0.1;
                cons[1*NK] = 0.0;
                cons[2*NK] = 0.0;
                cons[3*NK] = 0.125 / (ADIABATIC_GAMMA - 1.0);

            
                cons[0*NK+1] = 0.0;
                cons[1*NK+1] = 0.0;
                cons[2*NK+1] = 0.0;
                cons[3*NK+1] = 0.0;
            }*/
            
        }
    }
}

struct UpdateStruct
{
    int ni;
    real x0;
    real x1;
    real *conserved;
};

struct UpdateStruct update_struct_new(int ni, real x0, real x1)
{
    struct UpdateStruct update;
    update.ni = ni;
    update.x0 = x0;
    update.x1 = x1;

    update.conserved = (real*) malloc(ni * NK * NCONS * sizeof(real));

    return update;
}

void update_struct_del(struct UpdateStruct update)
{
    free(update.conserved);
}

void update_struct_set_conserved(struct UpdateStruct update, const real *conserved_host)
{
    int ni = update.ni;
    int num_zones = ni;

    memcpy(
        update.conserved,
        conserved_host,
        num_zones * NK * NCONS * sizeof(real)
    );
}

void update_struct_get_conserved(struct UpdateStruct update, real *conserved_host)
{
    int num_zones = update.ni;
    memcpy(conserved_host,
        update.conserved,
        num_zones * NK * NCONS * sizeof(real)
    );
}

void update_struct_do_advance_cons(struct UpdateStruct update, real dt)
{
    int ni = update.ni;

    const real dx = (update.x1 - update.x0) / update.ni;

    real cons[NCONS];
    real prim[NCONS];
    real consm[NCONS];
    real consp[NCONS];
    real primm[NCONS];
    real primp[NCONS];

    real flux[NCONS];

    real surface[NCONS][NK];
    real volume [NCONS][NK];

    for (int i = 0; i < ni; ++i)
    {
        memset(surface, 0.0, NCONS * NK * sizeof(real));

        int il = i - 1;
        int ir = i + 1;

        if (il == -1)
                il = 0;
        if (ir == update.ni)
                ir = update.ni - 1;

        real *cweight = &update.conserved[NK * NCONS * i];

        const real *cli = &update.conserved[NK * NCONS * il];
        const real *cri = &update.conserved[NK * NCONS * ir];

        // Sum over face nodes

        // Left Face

        for (int q = 0; q < NCONS; ++q)
        {
            consm[q] = 0.0;
            consp[q] = 0.0;

            for (int j = 0; j < NK; ++j)
            {
                // conserved variables at minus side of left face
                consm[q] += cli[NK * q + j] * node[NNODE - 1][j]; // right face of i-1
                // conserved variables at plus side of left face
                consp[q] += cweight[NK * q + j] * node[0][j];         // left face of i
            }
        }

        conserved_to_primitive(consm, primm);
        conserved_to_primitive(consp, primp);

        riemann_hlle(primm,primp,flux,0);

        for (int q = 0; q < NCONS; ++q)
        {
            for (int j = 0; j < NK; ++j)
            {
                surface[q][j] += flux[q] * node[0][j] * node[0][NK + NDIM * NK];
                //printf("i = %d q = %d j = %d surface = %f\n",i,q,j,surface[q][j]);
            }
        }

        // Right Face

        for (int q = 0; q < NCONS; ++q)
        {
            consm[q] = 0.0;
            consp[q] = 0.0;

            for (int j = 0; j < NK; ++j)
            {
                // conserved variables at minus side of right face
                consm[q] += cweight[NK * q + j] * node[NNODE - 1][j]; // right face of i
                // conserved variables at plus side of right face
                consp[q] += cri[NK * q + j] * node[0][j];         // left face of i + 1
            }
        }

        conserved_to_primitive(consm, primm);
        conserved_to_primitive(consp, primp);

        riemann_hlle(primm,primp,flux,0);

        for (int q = 0; q < NCONS; ++q)
        {
            for (int j = 0; j < NK; ++j)
            {
                surface[q][j] += flux[q] * node[NNODE-1][j] * node[NNODE-1][NK + NDIM * NK];
            }
        }   


        memset(volume,  0.0, NCONS * NK * sizeof(real));

        // Sum over cell quadrature nodes
        for (int nc = 1; nc < (NNODE - 1); ++nc)
        {
            // compute conserved variables at cell nodes

            for (int q = 0; q < NCONS; ++q)
            {
                cons[q] = 0.0;

                for (int j = 0; j < NK; ++j)
                {
                    cons[q] += cweight[NK * q + j] * node[nc][j];
                }
            }

            conserved_to_primitive(cons, prim);

            primitive_to_flux_vector(prim, flux, 0);

            for (int q = 0; q < NCONS; ++q)
            {
                for (int j = 0; j < NK; ++j)
                {
                    volume[q][j] += flux[q] * node[nc][NK + j];
                }
            }
        }

        for (int q = 0; q < NCONS; ++q)
        {
            for (int j = 0; j < NK; ++j)
            {
                cweight[NK * q + j] += (volume[q][j] - surface[q][j]) * dt / dx;
            }
        }

    }
}


int main()
{
    const int ni   = 2048;
    const int fold = 1;
    const real x0 = 0.0;
    const real x1 = 1.0;
    const real dx = (x1 - x0) / ni;

    real *conserved = (real*) malloc(ni * NK * NCONS * sizeof(real));
    struct UpdateStruct update = update_struct_new(ni, x0, x1);

    set_cell();

    initial_conserved(conserved, ni, x0, x1);
    update_struct_set_conserved(update, conserved);

    int iteration = 0;
    real time = 0.0;
    real dt = dx * 0.05;

    while (time < 0.01)
    {
        clock_t start = clock();

        for (int i = 0; i < fold; ++i)
        {
            update_struct_do_advance_cons(update, dt);

            time += dt;
            iteration += 1;
        }
        clock_t end = clock();

        real seconds = ((real) (end - start)) / CLOCKS_PER_SEC;
        real mzps = (ni / 1e6) / seconds * fold;
        printf("[%d] t=%.3e Mzps=%.2f\n", iteration, time, mzps);
    }

    update_struct_get_conserved(update, conserved);
    update_struct_del(update);

    FILE* outfile = fopen("dg.dat", "w");

    for (int i = 0; i < ni; ++i)
    {
        real *cweight = &conserved[NK * NCONS * i];
        real x = (i + 0.5) * dx;

        fprintf(outfile, "%f %f %f %f %f\n", x, cweight[0], cweight[NK], cweight[2*NK], cweight[3*NK]);
        //fprintf(outfile, "%f %f %f %f %f\n", x, cweight[0+1], cweight[NK+1], cweight[2*NK+1], cweight[3*NK+1]);
    }
    fclose(outfile);
    free(conserved);
    return 0;
}
