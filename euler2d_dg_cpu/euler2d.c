#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define ADIABATIC_GAMMA (7.0 / 5.0)
#define PI 3.14159265359
#define min2(a, b) (a) < (b) ? (a) : (b)
#define max2(a, b) (a) > (b) ? (a) : (b)
#define min3(a, b, c) min2(a, min2(b, c))
#define max3(a, b, c) max2(a, max2(b, c))
#define sign(x) copysign(1.0, x)
#define minabs(a, b, c) min3(fabs(a), fabs(b), fabs(c))
#define maxabs5(a, b, c, d, e) max2(max2(fabs(a), fabs(b)), max3(fabs(c), fabs(d), fabs(e)))

#define BETA 0.5
#define CFL 0.01

typedef double real;
#define square_root sqrt
#define power pow

#define NDIM 2
#define DG_ORDER 3
#define NFACE 3
#define NCELL 9
#define NK    6  // number of basis polynomials
#define NCONS 4  // number of conserved variables

#define SQRT_THREE square_root(3.0)
#define SQRT_FIVE  square_root(5.0)

struct Cells cell;

struct Cells
{
    struct Nodes
    {
        real xsi_x;         // normalized cell x coordinate [-1,1]
        real xsi_y;         // normalized cell y coordinate [-1,1]       
        real phi[NK];       // cell basis functions phi(xsi_x,xsi_y)
        real dphidx[NK];    // xsi_x derivatives of basis functions
        real dphidy[NK];    // xsi_y derivatives of basis functions
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

struct Cells set_cell(void)
{
    // reference: https://en.wikipedia.org/wiki/Gaussian_quadrature

    #if (DG_ORDER == 1)

        real xsi[NFACE] = {0.0};                                // Gaussian quadrature point
        real  gw[NFACE] = {2.0};                                // 1D Gaussian weight

    #elif (DG_ORDER == 2)

        real xsi[NFACE] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};      // Gaussian quadrature points
        real  gw[NFACE] = {1.0, 1.0};                           // 1D Gaussian weights

    #elif (DG_ORDER == 3)

        real xsi[NFACE] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)}; // Gaussian quadrature points
        real  gw[NFACE] = {5.0/9.0, 8.0/9.0, 5.0/9.0};          // 1D Gaussian weights

    #endif


    // cell nodes

    int nc = 0;

    for (int i = 0; i < NFACE; ++i)
    {
        for (int j = 0; j < NFACE; ++j)
        {   
            cell.node[nc].xsi_x     = xsi[i];
            cell.node[nc].xsi_y     = xsi[j];
            
            // nk = 0
            cell.node[nc].phi[0]    = p0(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidx[0] = p0_prime(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidy[0] = p0(xsi[i]) * p0_prime(xsi[j]);
            
            #if (DG_ORDER >= 2)
            
            // nk = 1 (y slope)
            cell.node[nc].phi[1]    = p0(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidx[1] = p0_prime(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidy[1] = p0(xsi[i]) * p1_prime(xsi[j]);
            
            // nk = 2 (x slope)
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
        // left face (xsi_x = -1.0)
        cell.node[nc].xsi_x     = -1.0;
        cell.node[nc].xsi_y     = xsi[n];

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
        cell.node[nc].xsi_x     = 1.0;
        cell.node[nc].xsi_y     = xsi[n];

        // nk = 0
        cell.faceri[n].phi[0]    = p0(1.0) * p0(xsi[n]);
        cell.faceri[n].dphidx[0] = p0_prime(1.0) * p0(xsi[n]);
        cell.faceri[n].dphidy[0] = p0(1.0) * p0_prime(xsi[n]);
        
        #if (DG_ORDER >= 2)
        
        // nk = 1
        cell.faceri[n].phi[1]    = p0(1.0) * p1(xsi[n]);
        cell.faceri[n].dphidx[1] = p0_prime(1.0) * p1(xsi[n]);
        cell.faceri[n].dphidy[1] = p0(1.0) * p1_prime(xsi[n]);

        // nk = 2
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
        cell.node[nc].xsi_x     = xsi[n];
        cell.node[nc].xsi_y     = -1.0;

        // nk = 0
        cell.facelj[n].phi[0]    = p0(xsi[n]) * p0(-1.0);
        cell.facelj[n].dphidx[0] = p0_prime(xsi[n]) * p0(-1.0);
        cell.facelj[n].dphidy[0] = p0(xsi[n]) * p0_prime(-1.0);
        
        #if (DG_ORDER >= 2)
        
        // nk = 1
        cell.facelj[n].phi[1]    = p0(xsi[n]) * p1(-1.0);
        cell.facelj[n].dphidx[1] = p0_prime(xsi[n]) * p1(-1.0);
        cell.facelj[n].dphidy[1] = p0(xsi[n]) * p1_prime(-1.0);

        // nk = 2
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
        cell.node[nc].xsi_x     = xsi[n];
        cell.node[nc].xsi_y     = 1.0;

        // nk = 0
        cell.facerj[n].phi[0]    = p0(xsi[n]) * p0(1.0);
        cell.facerj[n].dphidx[0] = p0_prime(xsi[n]) * p0(1.0);
        cell.facerj[n].dphidy[0] = p0(xsi[n]) * p0_prime(1.0);

        #if (DG_ORDER >= 2)

        // nk = 1
        cell.facerj[n].phi[1]    = p0(xsi[n]) * p1(1.0);
        cell.facerj[n].dphidx[1] = p0_prime(xsi[n]) * p1(1.0);
        cell.facerj[n].dphidy[1] = p0(xsi[n]) * p1_prime(1.0);

        // nk = 2
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

real minmod(real a, real b, real c)
{
        real x1 = fabs(sign(a) + sign(b)) * (sign(a) + sign(c));
        real x2 = minabs(a, b, c);
        real x  = 0.25 * x1 * x2;

        return x;
}

real minmodTVB(real w1, real w0l, real w0, real w0r, real dl)
{
    real a = w1 * SQRT_THREE;
    real b = (w0 - w0l) * BETA;
    real c = (w0r - w0) * BETA;

    const real M = 0.0; //Cockburn & Shu, JCP 141, 199 (1998) eq. 3.7 suggest M~50.0
    //const real Mtilde = 0.5; //Schaal+
    if (fabs(a) <= M * dl * dl)
    //if (fabs(a) <= Mtilde * dl)
    {        
        return w1;
    }
    else
    {
        real x1 = fabs(sign(a) + sign(b)) * (sign(a) + sign(c));
        real x2 = minabs(a, b, c);
        real x = (0.25 / SQRT_THREE) * x1 * x2;

        return x;
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

void riemann_llf(const real *pl, const real *pr, real *flux, int direction)
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

    for (int q = 0; q < NCONS; ++q)
    {
        flux[q] = 0.5 * (fl[q] + fr[q]) - 0.5 * ap * (ur[q] - ul[q]); // Check this
    }
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

    for (int q = 0; q < NCONS; ++q)
    {
        flux[q] = (fl[q] * ap - fr[q] * am - (ul[q] - ur[q]) * ap * am) / (ap - am);
    }
}

void riemann_hllc(const real *pl, const real *pr, real *flux, int direction)
{
    enum { d, px, py, e }; // Conserved
    enum { rho, vx, vy, p }; // Primitive

    real ul[NCONS];
    real ur[NCONS];
    real ulstar[NCONS];
    real urstar[NCONS];
    real fl[NCONS];
    real fr[NCONS];
    real al[2];
    real ar[2];

    const real vnl = primitive_to_velocity_component(pl, direction);
    const real vnr = primitive_to_velocity_component(pr, direction);

    primitive_to_conserved(pl, ul);
    primitive_to_conserved(pr, ur);
    primitive_to_flux_vector(pl, fl, direction);
    primitive_to_flux_vector(pr, fr, direction);
    primitive_to_outer_wavespeeds(pl, al, direction);
    primitive_to_outer_wavespeeds(pr, ar, direction);

    const real am = min3(0.0, al[0], ar[0]);
    const real ap = max3(0.0, al[1], ar[1]);

    real lc = (
        + (pr[p] - pr[rho] * vnr * (ap - vnr))
        - (pl[p] - pl[rho] * vnl * (am - vnl))) / (pl[rho] * (am - vnl) - pr[rho] * (ap - vnr));

    real ffl = pl[rho] * (am - vnl) / (am - lc);
    real ffr = pr[rho] * (ap - vnr) / (ap - lc);
    
    ulstar[d] = ffl;
    ulstar[e] = ffl * (ul[e] / pl[rho] + (lc - vnl) * (lc + pl[p] / (pl[rho] * (am - vnl))));
    ulstar[px] = ffl * ((lc - vnl) * (direction==0) + pl[vx]);
    ulstar[py] = ffl * ((lc - vnl) * (direction==1) + pl[vy]);

    urstar[d] = ffr;
    urstar[e] = ffr * (ur[e] / pr[rho] + (lc - vnl) * (lc + pr[p] / (pr[rho] * (ap - vnl))));
    urstar[px] = ffr * ((lc - vnl) * (direction==0) + pl[vx]);
    urstar[py] = ffr * ((lc - vnl) * (direction==1) + pl[vy]);

    const real s = 0.0; //stationary face s = x / t

    if      ( s <= am )       for (int i=0; i<NCONS; ++i) flux[i] = fl[i];
    else if ( am<s && s<=lc ) for (int i=0; i<NCONS; ++i) flux[i] = fl[i] + am * (ulstar[i] - ul[i]);
    else if ( lc<s && s<=ap ) for (int i=0; i<NCONS; ++i) flux[i] = fr[i] + ap * (urstar[i] - ur[i]);
    else if ( ap<s          ) for (int i=0; i<NCONS; ++i) flux[i] = fr[i];

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

            real prim[NCONS];
            real cons[NCONS];
            real rho, vx, vy, pressure;

            for (int l = 0; l < NK; ++l)
            {
                for (int q = 0; q < NCONS; ++q)
                {
                    wij[q * NK + l] = 0.0;
                }
            }

            // loop over cell quadrature points
            for (int qp = 0; qp < NCELL; ++qp){

                // global coordinates of the quadrature points
                real xq = x + (cell.node[qp].xsi_x / 2.0) * dx;
                real yq = y + (cell.node[qp].xsi_y / 2.0) * dy;

                // real r2 = power(xq - xmid, 2) + power(yq - ymid, 2);

                // Schaal shock tube problem 5.2, run to t=0.228
                if (xq < xmid)     
                {
                    rho = 1.0;
                    pressure = 1.0;
                }
                else
                {
                    rho = 0.125;
                    pressure = 0.1;
                }
                vx = 0.0;
                vy = 0.0;

                prim[0] = rho;
                prim[1] = vx;
                prim[2] = vy;
                prim[3] = pressure;

                primitive_to_conserved(prim, cons);

                for (int l = 0; l < NK; ++l)
                {
                    for (int q = 0; q < NCONS; ++q)
                    {
                        wij[q * NK + l] += 0.25 * cons[q] * cell.node[qp].phi[l] * cell.node[qp].gw;
                    }
                }
            }

            /*for (int q = 0; q < NCONS; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    if(j==4) printf("%d %d %d %f\n",i,q,l,wij[q * NK + l]);
                }
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
    real *dweights;
    real *trouble;
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

    update.weights  = (real*) malloc(ni * nj * NCONS * NK * sizeof(real));
    update.dweights = (real*) malloc(ni * nj * NCONS * NK * sizeof(real));
    update.trouble  = (real*) malloc(ni * nj * sizeof(real));

    return update;
}

void update_struct_del(struct UpdateStruct update)
{
    free(update.weights);
    free(update.dweights);
    free(update.trouble);
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

void update_struct_get_trouble(struct UpdateStruct update, real *trouble_host)
{
    int num_zones = update.ni * update.nj;
    memcpy(trouble_host,
        update.trouble,
        num_zones * sizeof(real)
    );
}

void compute_delta_weights(struct UpdateStruct update)
{
    int ni = update.ni;
    int nj = update.nj;

    real cons[NCONS];
    real prim[NCONS];
    real consm[NCONS];
    real consp[NCONS];
    real primm[NCONS];
    real primp[NCONS];

    real w1ij[NCONS];
    real dw0l[NCONS];
    real dw0r[NCONS];

    real flux_x[NCONS];
    real flux_y[NCONS];
    
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
            
            /* */ real *wij = &update.weights[NCONS * NK * (i  * nj + j )];
            /* */ real *dwij=&update.dweights[NCONS * NK * (i  * nj + j )];

            const real *wli = &update.weights[NCONS * NK * (il * nj + j )];
            const real *wri = &update.weights[NCONS * NK * (ir * nj + j )];
            const real *wlj = &update.weights[NCONS * NK * (i  * nj + jl)];
            const real *wrj = &update.weights[NCONS * NK * (i  * nj + jr)];
   
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
                        consm[q] += wli[NK * q + l] * cell.faceri[n].phi[l]; // right face of zone i-1 
                        consp[q] += wij[NK * q + l] * cell.faceli[n].phi[l]; // left face of zone i                     
                    }
                }

                conserved_to_primitive(consm, primm);
                conserved_to_primitive(consp, primp);
                
                riemann_hlle(primm, primp, flux_x, 0);
                
                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux_x[q] * cell.faceli[n].phi[l] * cell.faceli[n].gw;  
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
            
                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] -= flux_x[q] * cell.faceri[n].phi[l] * cell.faceri[n].gw;
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

            // Volume term
            #if (DG_ORDER >= 2)

            for (int n = 0; n < NCELL; ++n)
            {
                for (int q = 0; q < NCONS; ++q)
                {
                    cons[q] = 0.0;

                    for (int l = 0; l < NK; ++l)
                    {
                        cons[q] += wij[NK * q + l] * cell.node[n].phi[l];
                    }
                }

                conserved_to_primitive(cons, prim);

                primitive_to_flux_vector(prim, flux_x, 0);
                primitive_to_flux_vector(prim, flux_y, 1);

                for (int q = 0; q < NCONS; ++q)
                {
                    for (int l = 0; l < NK; ++l)
                    {
                        dwij[NK * q + l] += flux_x[q] * cell.node[n].dphidx[l] * cell.node[n].gw;
                        dwij[NK * q + l] += flux_y[q] * cell.node[n].dphidy[l] * cell.node[n].gw; 
                        //printf("q = %d l = %d %f %f\n", q, l, flux_x[q] * cell.node[n].dphidx[l] * cell.node[n].gw,flux_y[q] * cell.node[n].dphidy[l] * cell.node[n].gw);  
                    }
                }

            }
            #endif
        }
    }
}

void add_delta_weights(struct UpdateStruct update, real dt)
{
    int ni = update.ni;
    int nj = update.nj;
    const real dx = (update.x1 - update.x0) / update.ni;
    const real dy = (update.y1 - update.y0) / update.nj;

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
            
            real *wij = &update.weights[NCONS * NK * (i  * nj + j )];
            real *dwij= &update.dweights[NCONS * NK * (i  * nj + j )];

            for (int q = 0; q < NCONS; ++q)
            {
                for (int l = 0; l < NK; ++l)
                {
                    wij[ NK * q + l] += 0.5 * dwij[NK * q + l] * dt / dx; //assumes dy=dx
                }
            }

        }
    }
}

void mark_troubled_cells(struct UpdateStruct update)
{
    int ni = update.ni;
    int nj = update.nj;

    const real dx = (update.x1 - update.x0) / update.ni;
    const real dy = (update.y1 - update.y0) / update.nj;
    const real dvol = 4; // in xsi coordinate [-1,1] x [-1,1]

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
            
            /* */ real *wij = &update.weights[NCONS * NK * (i  * nj + j )];

            const real *wli = &update.weights[NCONS * NK * (il * nj + j )];
            const real *wri = &update.weights[NCONS * NK * (ir * nj + j )];
            const real *wlj = &update.weights[NCONS * NK * (i  * nj + jl)];
            const real *wrj = &update.weights[NCONS * NK * (i  * nj + jr)];

            /* */ real  *tij = &update.trouble[i * nj + j];

            // Troubled Cell Indicator G. Fu & C.-W. Shu (JCP, 347, 305 (2017))

            // eq. 2.3 compute max bar p and double bar p for rho and E

            // q=0 for rho

            real maxpj = maxabs5(wij[0 * NK + 0], wli[0 * NK + 0], wri[0 * NK + 0], wlj[0 * NK + 0], wrj[0 * NK + 0]);

            // first order
            real a = 4.0 * wli[0 * NK + 0];
            real b = 4.0 * wri[0 * NK + 0];
            real c = 4.0 * wlj[0 * NK + 0];
            real d = 4.0 * wrj[0 * NK + 0];

            #if (DG_ORDER >= 2)

            a += 8.0 * SQRT_THREE * wli[0 * NK + 2];
            b -= 8.0 * SQRT_THREE * wri[0 * NK + 2];
            c += 8.0 * SQRT_THREE * wlj[0 * NK + 1];
            d -= 8.0 * SQRT_THREE * wrj[0 * NK + 1];
            
            #endif

            #if (DG_ORDER >= 3)

            a += 24.0 * SQRT_FIVE * wli[0 * NK + 5];
            b += 24.0 * SQRT_FIVE * wri[0 * NK + 5];
            c += 24.0 * SQRT_FIVE * wlj[0 * NK + 3];
            d += 24.0 * SQRT_FIVE * wrj[0 * NK + 3];
 
            #endif

            real pbb_li = fabs(wij[0 * NK + 0] - a / dvol);
            real pbb_ri = fabs(wij[0 * NK + 0] - b / dvol);
            real pbb_lj = fabs(wij[0 * NK + 0] - c / dvol);
            real pbb_rj = fabs(wij[0 * NK + 0] - d / dvol);

            real trouble_rho = (pbb_li + pbb_ri + pbb_lj + pbb_rj) / maxpj;

            // q=3 for energy

            maxpj = maxabs5(wij[3 * NK + 0], wli[3 * NK + 0], wri[3 * NK + 0], wlj[3 * NK + 0], wrj[3 * NK + 0]);

            // first order
            a = 4.0 * wli[3 * NK + 0];
            b = 4.0 * wri[3 * NK + 0];
            c = 4.0 * wlj[3 * NK + 0];
            d = 4.0 * wrj[3 * NK + 0];
  
            #if (DG_ORDER >= 2)

            a += 8.0 * SQRT_THREE * wli[3 * NK + 2];
            b -= 8.0 * SQRT_THREE * wri[3 * NK + 2];
            c += 8.0 * SQRT_THREE * wlj[3 * NK + 1];
            d -= 8.0 * SQRT_THREE * wrj[3 * NK + 1];
             
            #endif

            #if (DG_ORDER >= 3)

            a += 24.0 * SQRT_FIVE * wli[3 * NK + 5];
            b += 24.0 * SQRT_FIVE * wri[3 * NK + 5];
            c += 24.0 * SQRT_FIVE * wlj[3 * NK + 3];
            d += 24.0 * SQRT_FIVE * wrj[3 * NK + 3];

            #endif

            pbb_li = fabs(wij[3 * NK + 0] - a / dvol);
            pbb_ri = fabs(wij[3 * NK + 0] - b / dvol);
            pbb_lj = fabs(wij[3 * NK + 0] - c / dvol);
            pbb_rj = fabs(wij[3 * NK + 0] - d / dvol);
 
            real trouble_e = (pbb_li + pbb_ri + pbb_lj + pbb_rj) / maxpj;

            *tij = max2( trouble_rho, trouble_e );
        }
    }
}

void limit_conserved_weights(struct UpdateStruct update)
{
    real cons[NCONS];
    real prim[NCONS];
    
    int ni = update.ni;
    int nj = update.nj;

    const real dx = (update.x1 - update.x0) / update.ni;
    const real dy = (update.y1 - update.y0) / update.nj;

    real wtilde[NK * NCONS];

    //const real ck = 0.6;       
    //const real ck = 0.03;        //k=1; see Fu & Shu Sec. 3
    const real ck = 0.03; //k=2; see Fu & Shu Sec. 3
    //const real ck = 0.03 * 4 ; //k=3; see Fu & Shu Sec. 3
    
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

            /* */ real *wij = &update.weights[NCONS * NK * (i  * nj + j )];

            const real *wli = &update.weights[NCONS * NK * (il * nj + j )];
            const real *wri = &update.weights[NCONS * NK * (ir * nj + j )];
            const real *wlj = &update.weights[NCONS * NK * (i  * nj + jl)];
            const real *wrj = &update.weights[NCONS * NK * (i  * nj + jr)];

            const real  *tij = &update.trouble[i * nj + j];

            if ( *tij > 0.4 )
            {
                for (int q = 0; q < NCONS; ++q)
                {
                    // real minmodTVB(real w1, real w0l, real w0, real w0r, real dl)
                    
                    // x slopes
                    wtilde[ NK * q + 2] = minmodTVB(wij[ NK * q + 2], wli[ NK * q + 0], wij[ NK * q + 0], wri[ NK * q + 0], dx);
                    
                    // y slopes 
                    wtilde[ NK * q + 1] = minmodTVB(wij[ NK * q + 1], wlj[ NK * q + 0], wij[ NK * q + 0], wrj[ NK * q + 0], dy);

                    
                    if ( (wtilde[ NK * q + 2] != wij[ NK * q + 2]) ||
                         (wtilde[ NK * q + 1] != wij[ NK * q + 1]) ) 
                    {
                        //if (j==0 && q==0) printf("cons 1: %d %f %f %f\n\n", i, wij[ NK * q + 2], wtilde[ NK * q + 2], wij[ NK * q + 5]);                        


                        wij[ NK * q + 2] = wtilde[ NK * q + 2];
                        wij[ NK * q + 1] = wtilde[ NK * q + 1];
                            
                        for (int l = 3; l < NK; ++l)
                        {
                            wij[ NK * q + l] = 0.0;
                        } 
                        //if (j==0 && q==0) printf("cons 2: %d %f %f %f\n\n", i, wij[ NK * q + 2], wtilde[ NK * q + 2], wij[ NK * q + 5]);                        
                       
                    }
                }
            }
        }
    }
}

void limit_characteristic_weights(struct UpdateStruct update)
{    
    int ni = update.ni;
    int nj = update.nj;

    const real dx = (update.x1 - update.x0) / update.ni;
    const real dy = (update.y1 - update.y0) / update.nj;

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

                  real *wij = &update.weights[NCONS * NK * (i  * nj + j )];

            const real *wli = &update.weights[NCONS * NK * (il * nj + j )];
            const real *wri = &update.weights[NCONS * NK * (ir * nj + j )];
            const real *wlj = &update.weights[NCONS * NK * (i  * nj + jl)];
            const real *wrj = &update.weights[NCONS * NK * (i  * nj + jr)];

            const real  *tij = &update.trouble[i * nj + j];

            real wtilde[NK * NCONS];

            real cons[NCONS];
            real prim[NCONS];

            real w0[NCONS];
            real w0l[NCONS];
            real w0r[NCONS];
            real w0b[NCONS];
            real w0t[NCONS];

            // slopes of conserved variables
            real w1[NCONS]; // y slopes
            real w2[NCONS]; // x slopes

            // slopes of characteristic variables
            real c1[NCONS]; // y slopes
            real c2[NCONS]; // x slopes

            // limited "tilde" slopes
            real w1t[NCONS];
            real w2t[NCONS]; 
            real c1t[NCONS];
            real c2t[NCONS]; 

            // characteristic version of difference of mean values (l=0) to neighbor zones
            real cl[NCONS]; // left
            real cr[NCONS]; // right   
            real cb[NCONS]; // bottom
            real ct[NCONS]; // top 

            if ( *tij > 0.03 )
            {
                for (int q = 0; q < NCONS; ++q)
                {
                    // mean values (l=0) of conserved variables in the cell and nearest neighbor cells
                    w0[q]  = wij[q * NK + 0];
                    w0l[q] = wli[q * NK + 0]; // left 
                    w0r[q] = wri[q * NK + 0]; // right
                    w0b[q] = wlj[q * NK + 0]; // bottom
                    w0t[q] = wrj[q * NK + 0]; // top

                    // slopes (l=1, l=2) of conserved variables in the cell
                    w1[q] =  wij[q * NK + 1]; // y slopes
                    w2[q] =  wij[q * NK + 2]; // x slopes
                }
            
                conserved_to_primitive(w0, prim);

                const real cs2 = primitive_to_sound_speed_squared(prim);
                real cs  = square_root(cs2);
                const real g1 = ADIABATIC_GAMMA - 1.0;
                real vx = prim[1];
                real vy = prim[2];
                real k = 0.5 * (vx * vx + vy * vy);
                real h = (cs2 / g1) + k;
                real phi = g1 * k;
                real beta = 1.0 / (2.0 * cs2);

                real lx[4][4]= {
                      {beta*(phi+cs*vx),  -beta*(g1*vx+cs),  -beta*g1*vy,       beta*g1},
                      {(1.0-2.0*beta*phi), 2.0*beta*g1*vx,   2.0*beta*g1*vy,    -2.0*beta*g1},
                      {beta*(phi-cs*vx),  -beta*(g1*vx-cs),  -beta*g1*vy,       beta*g1},
                      {vy,                      0.0,            -1.0,            0.0}};

                real ly[4][4] = {
                      {beta*(phi+cs*vy),  -beta*g1*vx,  -beta*(g1*vy+cs),       beta*g1},
                      {(1.0-2.0*beta*phi), 2.0*beta*g1*vx,   2.0*beta*g1*vy,    -2.0*beta*g1},
                      {beta*(phi-cs*vy),  -beta*g1*vx,  -beta*(g1*vy-cs),       beta*g1},
                      {-vx,                      1.0,            0.0,            0.0}};
                          
                real rx[4][4] = {
                      { 1.0,        1.0,        1.0,        0.0},
                      {(vx - cs),   vx,     (vx + cs),      0.0},
                      {  vy,        vy,         vy,        -1.0},
                      {(h - cs*vx), k,      (h + cs*vx),    -vy}};
                          
                real ry[4][4] = {
                      { 1.0,        1.0,        1.0,        0.0},
                      {vx,   vx,     vx,      1.0},
                      {  vy-cs,        vy,         vy+cs,        0.0},
                      {(h - cs*vy), k,      (h + cs*vy),    vx}};
/*
                real result[NCONS][NCONS];
                // Initializing elements of matrix mult to 0.
                for (int qi = 0; qi < NCONS; ++qi) {
                   for (int qj = 0; qj < NCONS; ++qj) {
                      result[qi][qj] = 0;
                   }
                }
             
                // Multiplying first and second matrices and storing it in result
                for (int qi = 0; qi < NCONS; ++qi) {
                   for (int qj = 0; qj < NCONS; ++qj) {
                      for (int qk = 0; qk < NCONS; ++qk) {
                         result[qi][qj] += ry[qi][qk] * ly[qk][qj];
                      }
                   }
                }

                printf("\nOutput Matrix:\n");
                for (int qi = 0; qi < NCONS; ++qi) 
                {
                    for (int qj = 0; qj < NCONS; ++qj) 
                    {
                        printf("%f  ", ry[qi][qj]);
                        if (qj == (NCONS - 1) ) printf("\n");
                    }
                }*/

                // convert to characteristic variables
                for (int qi = 0; qi < NCONS; ++qi)
                {
                    c2[qi] = 0.0;
                    c1[qi] = 0.0;
                    cl[qi] = 0.0;
                    cr[qi] = 0.0;
                    cb[qi] = 0.0;
                    ct[qi] = 0.0;

                    for (int qj = 0; qj < NCONS; ++qj)
                    {
                        c2[qi] += lx[qi][qj] *   w2[qj]; // x slopes
                        cl[qi] += lx[qi][qj] *  (w0[qj]  - w0l[qj]); // left difference
                        cr[qi] += lx[qi][qj] * (w0r[qj]  -  w0[qj]); // right difference

                        c1[qi] += ly[qi][qj] *   w1[qj]; // y slopes
                        cb[qi] += ly[qi][qj] *  (w0[qj]  - w0b[qj]); // bottom difference
                        ct[qi] += ly[qi][qj] * (w0t[qj]  -  w0[qj]); // top difference
                    }
                }

                // limit characteristic slopes (for l=1, l=2)
                for (int q = 0; q < NCONS; ++q)
                {
                    c1t[q] = minmod(SQRT_THREE * c1[q], BETA * cb[q], BETA * ct[q]) / SQRT_THREE;
                    c2t[q] = minmod(SQRT_THREE * c2[q], BETA * cl[q], BETA * cr[q]) / SQRT_THREE;
                    //if (j==0)printf("char: %f %f %f %f \n", SQRT_THREE * c2[q], BETA * cl[q], BETA * cr[q], SQRT_THREE * c2t[q]);
                }
                
                // compute limited conservative slopes (for l=1, l=2)
                for (int qi = 0; qi < NCONS; ++qi)
                {
                    w2t[qi] = 0.0;
                    w1t[qi] = 0.0;

                    for (int qj = 0; qj < NCONS; ++qj)
                    {
                        w1t[qi] += ry[qi][qj] * c1t[qj]; // y slopes
                        w2t[qi] += rx[qi][qj] * c2t[qj]; // x slopes
                    }
                }
                                
                for (int q = 0; q < NCONS; ++q)
                {
                    if ( (c2t[q] != c2[q]) || ( c1t[q] != c1[q]) )
                    {              
                        //if (j==0 && q==0)printf("char: %d %f %f %f %f \n\n",i, c2[q], c2t[q], w2[q], w2t[q]);

                        wij[NK * q + 2] = w2t[q];
                        wij[NK * q + 1] = w1t[q];
                            
                        for (int l = 3; l < NK; ++l)
                        {
                            wij[NK * q + l] = 0.0;
                        }                        
                    }
                }
            }
        }
    }
}


int main()
{
    const int ni = 64;
    const int nj = 64;
    const int fold = 1;
    const real x0 = 0.0;
    const real x1 = 1.0;
    const real y0 = 0.0;
    const real y1 = 1.0;
    const real dx = (x1 - x0) / ni;
    const real dy = (y1 - y0) / nj;

    real cons[NCONS];
    real prim[NCONS];

    cell = set_cell();

    real *weights_host = (real*) malloc(ni * nj * NCONS * NK * sizeof(real));
    real *trouble_host = (real*) malloc(ni * nj * sizeof(real));

    struct UpdateStruct update = update_struct_new(ni, nj, x0, x1, y0, y1);

    initial_weights(weights_host, ni, nj, x0, x1, y0, y1);
    update_struct_set_weights(update, weights_host);
    mark_troubled_cells(update);
    limit_characteristic_weights(update);

    int iteration = 0;
    real time = 0.0;
    real dt = dx * CFL;

    while (time < 0.228)
    {
        clock_t start = clock();

        for (int i = 0; i < fold; ++i)
        {
            mark_troubled_cells(update);
            //limit_conserved_weights(update);
            limit_characteristic_weights(update);
            compute_delta_weights(update);
            add_delta_weights(update, dt);
            
            time += dt;
            iteration += 1;
        }
        clock_t end = clock();

        real seconds = ((real) (end - start)) / CLOCKS_PER_SEC;
        real mcps = (ni * nj / 1e6) / seconds * fold;
        printf("[%d] t=%.3e Mcps=%.2f Mpps=%.2f Mnps=%.2f\n", iteration, time, mcps, NK * mcps, NCELL * mcps);
    }

    update_struct_get_weights(update, weights_host);
    update_struct_get_trouble(update, trouble_host);

    update_struct_del(update);

    FILE* outfile = fopen("euler2d.dat", "w");

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
            real *wij = &weights_host[NCONS * NK * (i * nj + j)];
            real *tij = &trouble_host[i * nj + j];

            real x = (i + 0.5) * dx;
            real y = (j + 0.5) * dy;

            for (int q = 0; q < NCONS; ++q)
            {
                cons[q] = wij[NK * q + 0]; 
            }

            conserved_to_primitive(cons, prim);
            //fprintf(outfile, "%f %f %f %f %f %f %f\n", x, y, cons[0], cons[1], cons[2], cons[3], *tij);
            fprintf(outfile, "%f %f %f %f %f %f %f\n", x, y, prim[0], prim[1], prim[2], prim[3], *tij);
            //fprintf(outfile, "%f %f %f %f %f %f\n", x, y, wij[0], wij[2], wij[5], *tij);
        }
    }
    fclose(outfile);
    free(weights_host);
    free(trouble_host);
    return 0;
}

