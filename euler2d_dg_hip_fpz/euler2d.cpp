#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <hip/hip_runtime.h>

#define ADIABATIC_GAMMA (5.0 / 3.0)
#define min2(a, b) (a) < (b) ? (a) : (b)
#define max2(a, b) (a) > (b) ? (a) : (b)

#define NK 3     // number of basis polynomials
#define NQFACE 2 // number of quadrature points per face
#define NQCELL NQFACE * NQFACE // number of cell quadrature points

#ifdef SINGLE
typedef float real;
#define square_root sqrtf
#define power powf
#else
typedef double real;
#define square_root sqrt
#define power pow
#endif

__device__ void conserved_to_primitive(const real *cons, real *prim)
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

__device__ __host__ void primitive_to_conserved(const real *prim, real *cons)
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

__device__ real primitive_to_velocity_component(const real *prim, int direction)
{
    switch (direction)
    {
        case 0: return prim[1];
        case 1: return prim[2];
        default: return 0.0;
    }
}

__device__ void primitive_to_flux_vector(const real *prim, real *flux, int direction)
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

__device__ real primitive_to_sound_speed_squared(const real *prim)
{
    const real rho = prim[0];
    const real pressure = prim[3];
    return ADIABATIC_GAMMA * pressure / rho;
}

__device__ void primitive_to_outer_wavespeeds(const real *prim, real *wavespeeds, int direction)
{
    const real cs = square_root(primitive_to_sound_speed_squared(prim));
    const real vn = primitive_to_velocity_component(prim, direction);
    wavespeeds[0] = vn - cs;
    wavespeeds[1] = vn + cs;
}

__device__ void riemann_hlle(const real *pl, const real *pr, real *flux, int direction)
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

void initial_conserved(real *conserved, int ni, int nj, real x0, real x1, real y0, real y1)
{
    real dx = (x1 - x0) / ni;
    real dy = (y1 - y0) / nj;

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
            real x = (i + 0.5) * dx;
            real y = (j + 0.5) * dy;
            real *cons = &conserved[NK * 4 * (i * nj + j)];
            real r2 = power(x - 0.5, 2) + power(y - 0.5, 2);

            if (square_root(r2) < 0.125)
            {
                cons[0] = 1.0; // mass density
                cons[1] = 0.0;
                cons[2] = 0.0;
                cons[3] = 0.0; // x momentum density
                cons[4] = 0.0;
                cons[5] = 0.0;
                cons[6] = 0.0; // y momentum density
                cons[7] = 0.0;
                cons[8] = 0.0;
                cons[9] = 1.0 / (ADIABATIC_GAMMA - 1.0); // energy density
                cons[10] = 0.0;
                cons[11] = 0.0;

            }
            else
            {
                cons[0] = 0.1;
                cons[1] = 0.0;
                cons[2] = 0.0;
                cons[3] = 0.0;
                cons[4] = 0.0;
                cons[5] = 0.0;
                cons[6] = 0.0;
                cons[7] = 0.0;
                cons[8] = 0.0;
                cons[9] = 0.125 / (ADIABATIC_GAMMA - 1.0);
                cons[10] = 0.0;
                cons[11] = 0.0;

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
    real *conserved;

    // point values of basis polynomials and their derivatives at the cell quadrature points

    real phi   [NQCELL][NK];
    real dphidx[NQCELL][NK];
    real dphidy[NQCELL][NK];
    real cell_weight[NQCELL]; // weights for Gaussian quadrature

    // point values of basis polynomials at the face quadrature points

    real phi_left  [NQFACE][NK];
    real phi_right [NQFACE][NK];
    real phi_bottom[NQFACE][NK];
    real phi_top   [NQFACE][NK];
    real face_weight[NQFACE];
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

    real quadrature_node[2] = {-1.0/square_root(3.0), 1.0/square_root(3.0)};

    int nq = 0;

    for (int i = 0; i < NQFACE; ++i)
    {
        for (int j = 0; j < NQFACE; ++j)
            {
                update.phi   [nq][0] = 1.0;
                update.dphidx[nq][0] = 0.0;
                update.dphidy[nq][0] = 0.0;
                update.phi   [nq][1] = (square_root(3.0) * quadrature_node[i]) * (1.0);
                update.dphidx[nq][1] =  square_root(3.0);
                update.dphidy[nq][1] = 0.0;
                update.phi   [nq][2] = (1.0) * (square_root(3.0) * quadrature_node[j]);
                update.dphidx[nq][2] = 0.0;
                update.dphidy[nq][2] =  square_root(3.0);
                update.cell_weight[nq] = 1.0;

                nq += 1;
            }
    }   


    for (int i = 0; i < NQFACE; ++i)
    {
        update.phi_left[i][0]   = 1.0;
        update.phi_left[i][1]   = square_root(3.0) * -1.0;
        update.phi_left[i][2]   = square_root(3.0) * quadrature_node[i];
        update.phi_right[i][0]  = 1.0;
        update.phi_right[i][1]  = square_root(3.0) *  1.0;
        update.phi_right[i][2]  = square_root(3.0) * quadrature_node[i];
        update.phi_bottom[i][0] = 1.0;
        update.phi_bottom[i][1] = square_root(3.0) * quadrature_node[i];
        update.phi_bottom[i][2] = square_root(3.0) * -1.0;
        update.phi_top[i][0]    = 1.0;
        update.phi_top[i][1]    = square_root(3.0) * quadrature_node[i];
        update.phi_top[i][2]    = square_root(3.0) *  1.0;    
        update.face_weight[i]   = 1.0;
    }

    hipMalloc(&update.conserved, ni * nj * NK * 4 * sizeof(real));

    return update;
}

void update_struct_del(struct UpdateStruct update)
{
    hipFree(update.conserved);
}

void update_struct_set_conserved(struct UpdateStruct update, const real *conserved_host)
{
    int ni = update.ni;
    int nj = update.nj;
    int num_zones = ni * nj;

    hipMemcpy(
        update.conserved,
        conserved_host,
        num_zones * 2 * 4 * sizeof(real),
        hipMemcpyHostToDevice
    );

}

__global__ void update_struct_do_advance_cons(struct UpdateStruct update, real dt)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int ni = update.ni;
    int nj = update.nj;
    const real dx = (update.x1 - update.x0) / update.ni;
    const real dy = (update.y1 - update.y0) / update.nj;

    if (i < ni && j < nj)
    {
        int il = i - 1;
        int ir = i + 1;
        int jl = j - 1;
        int jr = j + 1;

        if (il == -1)
            il = 0;
        if (ir == update.ni)
            ir = update.ni - 1;
        if (jl == -1)
            jl = 0;
        if (jr == update.nj)
            jr = update.nj - 1;

        real       *czone   = &update.conserved[NK * 4 * (i  * nj +  j)];

        const real *cleft   = &update.conserved[NK * 4 * (il * nj + j)];
        const real *cright  = &update.conserved[NK * 4 * (ir * nj + j)];
        const real *cbottom = &update.conserved[NK * 4 * (i  * nj + jl)];
        const real *ctop    = &update.conserved[NK * 4 * (i  * nj + jr)];       

        real cons[4];
        real prim[4];
        real flux_x[4];
        real flux_y[4];

        real dw[4][NK];
    
        // "Cell" Term

        // loop over cell quadrature nodes             
        for (int nq = 0; nq < NQCELL; ++nq)
        {
            // calculate cons at the nodes
            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    cons[q] = czone[q * NK + l] * update.phi[nq][l];
                }
            }

            conserved_to_primitive(cons, prim);
            primitive_to_flux_vector(prim, flux_x, 0);
            primitive_to_flux_vector(prim, flux_y, 1);

            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    dw[q][l] +=  update.face_weight[nq] * (flux_x[q] * update.dphidx[l][nq] + flux_y[q] * update.dphidy[l][nq]);
                }
            }   
        }

        // Face terms

        real consp[4];
        real primp[4];
        real consm[4];
        real primm[4];

        //Left Face

        // loop over face quadrature nodes             
        for (int nq = 0; nq < NQFACE; ++nq)
        {
            // calculate cons at both sides of face nodes
            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    consp[q] = czone[q * NK + l] * update.phi_left[nq][l];
                    consm[q] = cleft[q * NK + l] * update.phi_right[nq][l];
                }
            }

            conserved_to_primitive(consp, primp);
            conserved_to_primitive(consm, primm);

            riemann_hlle(primm, primp, flux_x, 0);         

            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    dw[q][l] += flux_x[q] * update.phi_left[nq][l] * update.face_weight[nq];
                }
            }
        }

        //Right Face

        // loop over face quadrature nodes             
        for (int nq = 0; nq < NQFACE; ++nq)
        {
            // calculate cons at both sides of face nodes
            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    consm[q] = czone[q * NK + l] * update.phi_right[nq][l];
                    consp[q] = cright[q * NK + l] * update.phi_left[nq][l];
                }
            }

            conserved_to_primitive(consp, primp);
            conserved_to_primitive(consm, primm);

            riemann_hlle(primm, primp, flux_x, 0);         

            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    dw[q][l] += -1.0 * flux_x[q] * update.phi_right[nq][l] * update.face_weight[nq];
                }
            }
        }  

        //Bottom Face

        // loop over face quadrature nodes             
        for (int nq = 0; nq < NQFACE; ++nq)
        {
            // calculate cons at both sides of face nodes
            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    consp[q] = czone[q * NK + l] * update.phi_bottom[nq][l];
                    consm[q] = cbottom[q * NK + l] * update.phi_top[nq][l];
                }
            }

            conserved_to_primitive(consp, primp);
            conserved_to_primitive(consm, primm);

            riemann_hlle(primm, primp, flux_y, 1);         

            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    dw[q][l] += flux_y[q] * update.phi_left[nq][l] * update.face_weight[nq];
                }
            }
        }

        //Top Face

        // loop over face quadrature nodes             
        for (int nq = 0; nq < NQFACE; ++nq)
        {
            // calculate cons at both sides of face nodes
            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    consm[q] = czone[q * NK + l] * update.phi_top[nq][l];
                    consp[q] = ctop[q * NK + l] * update.phi_bottom[nq][l];
                }
            }

            conserved_to_primitive(consp, primp);
            conserved_to_primitive(consm, primm);

            riemann_hlle(primm, primp, flux_y, 1);         

            for (int q = 0; q < 4; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    dw[q][l] += -1.0 * flux_y[q] * update.phi_top[nq][l] * update.face_weight[nq];
                }
            }
        }


    }
}

int main()
{
    const int ni = 4096;
    const int nj = 4096;
    const int fold = 10;
    const real x0 = 0.0;
    const real x1 = 1.0;
    const real y0 = 0.0;
    const real y1 = 1.0;
    const real dx = (x1 - x0) / ni;
    const real dy = (y1 - y0) / nj;

    real *conserved = (real*) malloc(ni * nj * NK * 4 * sizeof(real));
    struct UpdateStruct update = update_struct_new(ni, nj, x0, x1, y0, y1);

    initial_conserved(conserved, ni, nj, x0, x1, y0, y1);
    update_struct_set_conserved(update, conserved);

    int iteration = 0;
    real time = 0.0;
    real dt = min2(dx, dy) * 0.05;

    int thread_per_dim_i = 8;
    int thread_per_dim_j = 8;
    int blocks_per_dim_i = (ni + thread_per_dim_i - 1) / thread_per_dim_i;
    int blocks_per_dim_j = (nj + thread_per_dim_j - 1) / thread_per_dim_j;
    dim3 group_size = dim3(blocks_per_dim_i, blocks_per_dim_j, 1);
    dim3 block_size = dim3(thread_per_dim_i, thread_per_dim_j, 1);

    while (time < 0.2)
    {
        clock_t start = clock();

        for (int i = 0; i < fold; ++i)
        {
            update_struct_do_advance_cons<<<group_size, block_size>>>(update, dt);
            time += dt;
            iteration += 1;
        }
        hipDeviceSynchronize();
        clock_t end = clock();

        real seconds = ((real) (end - start)) / CLOCKS_PER_SEC;
        real mzps = (ni * nj / 1e6) / seconds * fold;
        real mnps = mzps / 4;
        printf("[%d] t=%.3e Mzps=%.2f Mnps=%.2f\n", iteration, time, mzps, mnps);
    }

    update_struct_del(update);

    FILE* outfile = fopen("euler2d.dat", "w");
    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
	  //            real *prim = &primitive[4 * (i * nj + j)];
          //  real x = (i + 0.5) * dx;
          //  real y = (j + 0.5) * dy;
          //  fprintf(outfile, "%f %f %f %f %f %f\n", x, y, prim[0], prim[1], prim[2], prim[3]);
        }
    }
    fclose(outfile);

    //free(primitive);
    return 0;
}
