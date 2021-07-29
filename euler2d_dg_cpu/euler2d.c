#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define ADIABATIC_GAMMA (5.0 / 3.0)
#define PI 3.14159265359
#define min2(a, b) (a) < (b) ? (a) : (b)
#define max2(a, b) (a) > (b) ? (a) : (b)
#define min3(a, b, c) min2(a, min2(b, c))
#define max3(a, b, c) max2(a, max2(b, c))
#define sign(x) copysign(1.0, x)
#define minabs(a, b, c) min3(fabs(a), fabs(b), fabs(c))
#define maxabs(a, b, c, d, e) max2(max2(fabs(a), fabs(b)), max3(fabs(c), fabs(d), fabs(e)))

#define BETA 1.0
#define CFL 0.03

typedef double real;
#define square_root sqrt
#define power pow

#define NDIM 2
#define DG_ORDER 2
#define NFACE 2
#define NCELL 4
#define NK    3  // number of basis polynomials
#define NCONS 4  // number of conserved variables

#define SQRT_THREE square_root(3.0)

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

real p4(const real xsi)
{
    return 3.0 / 8.0 * (35.0 * xsi * xsi * xsi * xsi - 30.0 * xsi * xsi + 3.0);
}

real p4_prime(const real xsi)
{
    return 3.0 / 8.0 * (140.0 * xsi * xsi * xsi - 60.0 * xsi);
}

struct Cells set_cell(void)
{
    // reference: https://en.wikipedia.org/wiki/Gaussian_quadrature

    #if (DG_ORDER == 1)

        real xsi[NFACE] = {0.0};                                // Gaussian quadrature point
        real  gw[NFACE] = {2.0};                                // 1D Gaussian weight

        real xsi_gl[NFACE+2] = {-1.0, 0.0, 1.0};                // Gauss-Lobatto quadrature points
        real  glw[NFACE+2] = {1.0/3.0, 4.0/3.0, 1.0/3.0};       // 1D Gauss-Lobatto weights

    #elif (DG_ORDER == 2)

        real xsi[NFACE] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};      // Gaussian quadrature points
        real  gw[NFACE] = {1.0, 1.0};                           // 1D Gaussian weights

        real xsi_gl[NFACE+2] = {-1.0, -sqrt(1.0/5.0), sqrt(1.0/5.0), 1.0}; // Gauss-Lobatto quadrature points
        real  glw[NFACE+2] = {1.0/6.0, 5.0/6.0, 5.0/6.0, 1.0/6.0};         // 1D Gauss-Lobatto weights
    
    #elif (DG_ORDER == 3)

        real xsi[NFACE] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)}; // Gaussian quadrature points
        real  gw[NFACE] = {5.0/9.0, 8.0/9.0, 5.0/9.0};          // 1D Gaussian weights
    
        real xsi_gl[NFACE+2] = {-1.0,-sqrt(3.0/7.0),0.0,sqrt(3.0/7.0),1.0};                    // Gauss-Lobatto quadrature points
        real  glw[NFACE+2] = {1.0/10.0, 49.0/90.0, 32.0/45.0, 49.0/90.0, 1.0/10.0};             // 1D Gauss-Lobatto weights

    #elif (DG_ORDER == 4)

        real xsi[NFACE] = {-sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0)), -sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)), sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)), sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0))}; // Gaussian quadrature points
        real  gw[NFACE] = {(18.0 - sqrt(30.0))/36.0, (18.0 + sqrt(30.0))/36.0, (18.0 + sqrt(30.0))/36.0, (18.0 - sqrt(30.0))/36.0};  // 1D Gaussian weights    
    
        real xsi_gl[NFACE+2] = {-1.0, -sqrt(1.0/3.0+2.0*sqrt(7.0)/21.0), -sqrt(1.0/3.0-2.0*sqrt(7.0)/21.0), sqrt(1.0/3.0-2.0*sqrt(7.0)/21.0), sqrt(1.0/3.0+2.0*sqrt(7.0)/21.0),1.0}; // Gauss-Lobatto quadrature points
        real  glw[NFACE+2] = {1.0/15.0, (14.0-sqrt(7.0))/30.0, (14.0+sqrt(7.0))/30.0, (14.0+sqrt(7.0))/30.0, (14.0-sqrt(7.0))/30.0, 1.0/15.0};   // 1D Gauss-Lobatto weights
    
    #elif (DG_ORDER == 5)

        real xsi[NFACE] = {-1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0)), -1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)), 0.0, 1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)), 1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0))}; // Gaussian quadrature points
        real  gw[NFACE] = {(322.0-13.0*sqrt(70.0))/900.0, (322.0+13.0*sqrt(70.0))/900.0, 128.0/225.0, (322.0+13.0*sqrt(70.0))/900.0, (322.0-13.0*sqrt(70.0))/900.0};          // 1D Gaussian weights    
    
        real xsi_gl[NFACE+2] = {-1.0, -sqrt(5.0/11.0+2.0/11.0*sqrt(5.0/3.0)), -sqrt(5.0/11.0-2.0/11.0*sqrt(5.0/3.0)), 0.0, sqrt(5.0/11.0-2.0/11.0*sqrt(5.0/3.0)), sqrt(5.0/11.0+2.0/11.0*sqrt(5.0/3.0)), 1.0};  // Gauss-Lobatto quadrature points
        real  glw[NFACE+2] = {1.0/21.0, (124.0-7.0*sqrt(15.0))/350.0, (124.0+7.0*sqrt(15.0))/350.0, 256.0/525.0, (124.0+7.0*sqrt(15.0))/350.0, (124.0-7.0*sqrt(15.0))/350.0, 1.0/21.0}; // 1D Gauss-Lobatto weights

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

            #if (DG_ORDER >= 4)

            // nk = 6
            cell.node[nc].phi[6]    = p0(xsi[i]) * p3(xsi[j]);
            cell.node[nc].dphidx[6] = p0_prime(xsi[i]) * p3(xsi[j]);
            cell.node[nc].dphidy[6] = p0(xsi[i]) * p3_prime(xsi[j]);
            
            // nk = 7
            cell.node[nc].phi[7]    = p1(xsi[i]) * p2(xsi[j]);
            cell.node[nc].dphidx[7] = p1_prime(xsi[i]) * p2(xsi[j]);
            cell.node[nc].dphidy[7] = p1(xsi[i]) * p2_prime(xsi[j]);
            
            // nk = 8
            cell.node[nc].phi[8]    = p2(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidx[8] = p2_prime(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidy[8] = p2(xsi[i]) * p1_prime(xsi[j]);

            // nk = 9
            cell.node[nc].phi[9]    = p3(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidx[9] = p3_prime(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidy[9] = p3(xsi[i]) * p0_prime(xsi[j]);
            
            #endif

            #if (DG_ORDER >= 5)

            // nk = 10
            cell.node[nc].phi[10]    = p0(xsi[i]) * p4(xsi[j]);
            cell.node[nc].dphidx[10] = p0_prime(xsi[i]) * p4(xsi[j]);
            cell.node[nc].dphidy[10] = p0(xsi[i]) * p4_prime(xsi[j]);

            // nk = 11
            cell.node[nc].phi[11]    = p1(xsi[i]) * p3(xsi[j]);
            cell.node[nc].dphidx[11] = p1_prime(xsi[i]) * p3(xsi[j]);
            cell.node[nc].dphidy[11] = p1(xsi[i]) * p3_prime(xsi[j]);            

            // nk = 12
            cell.node[nc].phi[12]    = p2(xsi[i]) * p2(xsi[j]);
            cell.node[nc].dphidx[12] = p2_prime(xsi[i]) * p2(xsi[j]);
            cell.node[nc].dphidy[12] = p2(xsi[i]) * p2_prime(xsi[j]);

            // nk = 13
            cell.node[nc].phi[13]    = p3(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidx[13] = p3_prime(xsi[i]) * p1(xsi[j]);
            cell.node[nc].dphidy[13] = p3(xsi[i]) * p1_prime(xsi[j]);

            // nk = 14
            cell.node[nc].phi[14]    = p4(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidx[14] = p4_prime(xsi[i]) * p0(xsi[j]);
            cell.node[nc].dphidy[14] = p4(xsi[i]) * p0_prime(xsi[j]);

            #endif
            
            cell.node[nc].gw = gw[i] * gw[j];   // 2D Gaussian weight

            //printf("%d %d %d %f %f %f\n",i,j,nc,cell.node[nc].phi[0], cell.node[nc].phi[1],cell.node[nc].phi[2]);
            //printf("%d %d %d %f %f %f\n",i,j,nc,cell.node[nc].dphidx[0], cell.node[nc].dphidx[1],cell.node[nc].dphidx[2]);
            //printf("%d %d %d %f %f %f\n\n",i,j,nc,cell.node[nc].dphidy[0], cell.node[nc].dphidy[1],cell.node[nc].dphidy[2]);

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

        #if (DG_ORDER >= 4)

        // nk = 6
        cell.node[nc].phi[6]    = p0(-1.0) * p3(xsi[n]);
        cell.node[nc].dphidx[6] = p0_prime(-1.0) * p3(xsi[n]);
        cell.node[nc].dphidy[6] = p0(-1.0) * p3_prime(xsi[n]);
            
        // nk = 7
        cell.node[nc].phi[7]    = p1(-1.0) * p2(xsi[n]);
        cell.node[nc].dphidx[7] = p1_prime(-1.0) * p2(xsi[n]);
        cell.node[nc].dphidy[7] = p1(-1.0) * p2_prime(xsi[n]);
            
        // nk = 8
        cell.node[nc].phi[8]    = p2(-1.0) * p1(xsi[n]);
        cell.node[nc].dphidx[8] = p2_prime(-1.0) * p1(xsi[n]);
        cell.node[nc].dphidy[8] = p2(-1.0) * p1_prime(xsi[n]);

        // nk = 9
        cell.node[nc].phi[9]    = p3(-1.0) * p0(xsi[n]);
        cell.node[nc].dphidx[9] = p3_prime(-1.0) * p0(xsi[n]);
        cell.node[nc].dphidy[9] = p3(-1.0) * p0_prime(xsi[n]);
        
        #endif

        #if (DG_ORDER >= 5)

        // nk = 10
        cell.node[nc].phi[10]    = p0(-1.0) * p4(xsi[n]);
        cell.node[nc].dphidx[10] = p0_prime(-1.0) * p4(xsi[n]);
        cell.node[nc].dphidy[10] = p0(-1.0) * p4_prime(xsi[n]);

        // nk = 11
        cell.node[nc].phi[11]    = p1(-1.0) * p3(xsi[n]);
        cell.node[nc].dphidx[11] = p1_prime(-1.0) * p3(xsi[n]);
        cell.node[nc].dphidy[11] = p1(-1.0) * p3_prime(xsi[n]);            

        // nk = 12
        cell.node[nc].phi[12]    = p2(-1.0) * p2(xsi[n]);
        cell.node[nc].dphidx[12] = p2_prime(-1.0) * p2(xsi[n]);
        cell.node[nc].dphidy[12] = p2(-1.0) * p2_prime(xsi[n]);

        // nk = 13
        cell.node[nc].phi[13]    = p3(-1.0) * p1(xsi[n]);
        cell.node[nc].dphidx[13] = p3_prime(-1.0) * p1(xsi[n]);
        cell.node[nc].dphidy[13] = p3(-1.0) * p1_prime(xsi[n]);

        // nk = 14
        cell.node[nc].phi[14]    = p4(-1.0) * p0(xsi[n]);
        cell.node[nc].dphidx[14] = p4_prime(-1.0) * p0(xsi[n]);
        cell.node[nc].dphidy[14] = p4(-1.0) * p0_prime(xsi[n]);

        #endif

        cell.faceli[n].gw = gw[n] * (-1.0); // 1D Gaussian weight * nhat

        // right face (x = +1)

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

        #if (DG_ORDER >= 4)

        // nk = 6
        cell.node[nc].phi[6]    = p0(1.0) * p3(xsi[n]);
        cell.node[nc].dphidx[6] = p0_prime(1.0) * p3(xsi[n]);
        cell.node[nc].dphidy[6] = p0(1.0) * p3_prime(xsi[n]);
            
        // nk = 7
        cell.node[nc].phi[7]    = p1(1.0) * p2(xsi[n]);
        cell.node[nc].dphidx[7] = p1_prime(1.0) * p2(xsi[n]);
        cell.node[nc].dphidy[7] = p1(1.0) * p2_prime(xsi[n]);
            
        // nk = 8
        cell.node[nc].phi[8]    = p2(1.0) * p1(xsi[n]);
        cell.node[nc].dphidx[8] = p2_prime(1.0) * p1(xsi[n]);
        cell.node[nc].dphidy[8] = p2(1.0) * p1_prime(xsi[n]);

        // nk = 9
        cell.node[nc].phi[9]    = p3(1.0) * p0(xsi[n]);
        cell.node[nc].dphidx[9] = p3_prime(1.0) * p0(xsi[n]);
        cell.node[nc].dphidy[9] = p3(1.0) * p0_prime(xsi[n]);
        
        #endif

        #if (DG_ORDER >= 5)

        // nk = 10
        cell.node[nc].phi[10]    = p0(1.0) * p4(xsi[n]);
        cell.node[nc].dphidx[10] = p0_prime(1.0) * p4(xsi[n]);
        cell.node[nc].dphidy[10] = p0(1.0) * p4_prime(xsi[n]);

        // nk = 11
        cell.node[nc].phi[11]    = p1(1.0) * p3(xsi[n]);
        cell.node[nc].dphidx[11] = p1_prime(1.0) * p3(xsi[n]);
        cell.node[nc].dphidy[11] = p1(1.0) * p3_prime(xsi[n]);            

        // nk = 12
        cell.node[nc].phi[12]    = p2(1.0) * p2(xsi[n]);
        cell.node[nc].dphidx[12] = p2_prime(1.0) * p2(xsi[n]);
        cell.node[nc].dphidy[12] = p2(1.0) * p2_prime(xsi[n]);

        // nk = 13
        cell.node[nc].phi[13]    = p3(1.0) * p1(xsi[n]);
        cell.node[nc].dphidx[13] = p3_prime(1.0) * p1(xsi[n]);
        cell.node[nc].dphidy[13] = p3(1.0) * p1_prime(xsi[n]);

        // nk = 14
        cell.node[nc].phi[14]    = p4(1.0) * p0(xsi[n]);
        cell.node[nc].dphidx[14] = p4_prime(1.0) * p0(xsi[n]);
        cell.node[nc].dphidy[14] = p4(1.0) * p0_prime(xsi[n]);

        #endif

        cell.faceri[n].gw = gw[n] * (1.0); // 1D Gaussian weight * nhat

        // bottom face (y = -1)

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

        #if (DG_ORDER >= 4)

        // nk = 6
        cell.node[nc].phi[6]    = p0(xsi[n]) * p3(-1.0);
        cell.node[nc].dphidx[6] = p0_prime(xsi[n]) * p3(-1.0);
        cell.node[nc].dphidy[6] = p0(xsi[n]) * p3_prime(-1.0);
            
        // nk = 7
        cell.node[nc].phi[7]    = p1(xsi[n]) * p2(-1.0);
        cell.node[nc].dphidx[7] = p1_prime(xsi[n]) * p2(-1.0);
        cell.node[nc].dphidy[7] = p1(xsi[n]) * p2_prime(-1.0);
            
        // nk = 8
        cell.node[nc].phi[8]    = p2(xsi[n]) * p1(-1.0);
        cell.node[nc].dphidx[8] = p2_prime(xsi[n]) * p1(-1.0);
        cell.node[nc].dphidy[8] = p2(xsi[n]) * p1_prime(-1.0);

        // nk = 9
        cell.node[nc].phi[9]    = p3(xsi[n]) * p0(-1.0);
        cell.node[nc].dphidx[9] = p3_prime(xsi[n]) * p0(-1.0);
        cell.node[nc].dphidy[9] = p3(xsi[n]) * p0_prime(-1.0);
            
        #endif

        #if (DG_ORDER >= 5)

        // nk = 10
        cell.node[nc].phi[10]    = p0(xsi[n]) * p4(-1.0);
        cell.node[nc].dphidx[10] = p0_prime(xsi[n]) * p4(-1.0);
        cell.node[nc].dphidy[10] = p0(xsi[n]) * p4_prime(-1.0);

        // nk = 11
        cell.node[nc].phi[11]    = p1(xsi[n]) * p3(-1.0);
        cell.node[nc].dphidx[11] = p1_prime(xsi[n]) * p3(-1.0);
        cell.node[nc].dphidy[11] = p1(xsi[n]) * p3_prime(-1.0);            

        // nk = 12
        cell.node[nc].phi[12]    = p2(xsi[n]) * p2(-1.0);
        cell.node[nc].dphidx[12] = p2_prime(xsi[n]) * p2(-1.0);
        cell.node[nc].dphidy[12] = p2(xsi[n]) * p2_prime(-1.0);

        // nk = 13
        cell.node[nc].phi[13]    = p3(xsi[n]) * p1(-1.0);
        cell.node[nc].dphidx[13] = p3_prime(xsi[n]) * p1(-1.0);
        cell.node[nc].dphidy[13] = p3(xsi[n]) * p1_prime(-1.0);

        // nk = 14
        cell.node[nc].phi[14]    = p4(xsi[n]) * p0(-1.0);
        cell.node[nc].dphidx[14] = p4_prime(xsi[n]) * p0(-1.0);
        cell.node[nc].dphidy[14] = p4(xsi[n]) * p0_prime(-1.0);

        #endif

        cell.facelj[n].gw = gw[n] * (-1.0); // 1D Gaussian weight * nhat

        // top face (y = +1)

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

        #if (DG_ORDER >= 4)

        // nk = 6
        cell.node[nc].phi[6]    = p0(xsi[n]) * p3(1.0);
        cell.node[nc].dphidx[6] = p0_prime(xsi[n]) * p3(1.0);
        cell.node[nc].dphidy[6] = p0(xsi[n]) * p3_prime(1.0);
            
        // nk = 7
        cell.node[nc].phi[7]    = p1(xsi[n]) * p2(1.0);
        cell.node[nc].dphidx[7] = p1_prime(xsi[n]) * p2(1.0);
        cell.node[nc].dphidy[7] = p1(xsi[n]) * p2_prime(1.0);
            
        // nk = 8
        cell.node[nc].phi[8]    = p2(xsi[n]) * p1(1.0);
        cell.node[nc].dphidx[8] = p2_prime(xsi[n]) * p1(1.0);
        cell.node[nc].dphidy[8] = p2(xsi[n]) * p1_prime(1.0);

        // nk = 9
        cell.node[nc].phi[9]    = p3(xsi[n]) * p0(1.0);
        cell.node[nc].dphidx[9] = p3_prime(xsi[n]) * p0(1.0);
        cell.node[nc].dphidy[9] = p3(xsi[n]) * p0_prime(1.0);
            
        #endif

        #if (DG_ORDER >= 5)

        // nk = 10
        cell.node[nc].phi[10]    = p0(xsi[n]) * p4(1.0);
        cell.node[nc].dphidx[10] = p0_prime(xsi[n]) * p4(1.0);
        cell.node[nc].dphidy[10] = p0(xsi[n]) * p4_prime(1.0);

        // nk = 11
        cell.node[nc].phi[11]    = p1(xsi[n]) * p3(1.0);
        cell.node[nc].dphidx[11] = p1_prime(xsi[n]) * p3(1.0);
        cell.node[nc].dphidy[11] = p1(xsi[n]) * p3_prime(1.0);            

        // nk = 12
        cell.node[nc].phi[12]    = p2(xsi[n]) * p2(1.0);
        cell.node[nc].dphidx[12] = p2_prime(xsi[n]) * p2(1.0);
        cell.node[nc].dphidy[12] = p2(xsi[n]) * p2_prime(1.0);

        // nk = 13
        cell.node[nc].phi[13]    = p3(xsi[n]) * p1(1.0);
        cell.node[nc].dphidx[13] = p3_prime(xsi[n]) * p1(1.0);
        cell.node[nc].dphidy[13] = p3(xsi[n]) * p1_prime(1.0);

        // nk = 14
        cell.node[nc].phi[14]    = p4(xsi[n]) * p0(1.0);
        cell.node[nc].dphidx[14] = p4_prime(xsi[n]) * p0(1.0);
        cell.node[nc].dphidy[14] = p4(xsi[n]) * p0_prime(1.0);

        #endif

        cell.facerj[n].gw = gw[n] * (1.0); // 1D Gaussian weight * nhat
    }
    return cell;
}

real minmod(real w1, real w0l, real w0, real w0r, real dx)
{
    real a = w1 / SQRT_THREE;
    real b = (w0 - w0l) * BETA;
    real c = (w0r - w0) * BETA;

    const real M = 0.0; //Cockburn & Shu, JCP 141, 199 (1998) eq. 3.7 suggest M~50.0

    if (fabs(a) <= M * dx * dx)
    {        
        return w1;
    }
    else
    {
        return 0.25 * SQRT_THREE * fabs(sign(a) + sign(b)) * (sign(a) + sign(c)) * minabs(a, b, c);
    }
}

/*

void zero_weights_to_characteristic_x(const real *cons, real *charx)
{
    real prim[NCONS];

    //conserved_to_primitive(cons,prim);

    real rho = 1.5;
    real vx = 0.5;
    real vy = 0.25;
    real pressure = 0.75;

    prim[0] = rho;
    prim[1] = vx;
    prim[2] = vy;
    prim[3] = pressure;

    const real g1 = ADIABATIC_GAMMA - 1.0;

    real cs2 = primitive_to_sound_speed_squared(prim);

    real cs  = square_root(cs2);
    
    real k = 0.5 * (vx * vx + vy * vy);
    real h = (cs2 / g1) + k;
    real phi = g1 * k;
    real beta = 1.0 / (2.0 * cs2);

    // Schaal+ Appendix E

    real lx[4][4] = {
          {beta*(phi+cs*vx),  -beta*(g1*vx+cs),  -beta*g1*vy,       beta*g1},
          {(1.0-2.0*beta*phi), 2.0*beta*g1*vx,   2.0*beta*g1*vy,    -2.0*beta*g1},
          {beta*(phi-cs*vx),  -beta*(g1*vx-cs),  -beta*g1*vy,       beta*g1},
          {vy,                      0.0,            -1.0,            0.0}
                    };

    real ly[4][4] = {
          {beta*(phi+cs*vy),  -beta*g1*vx,  -beta*(g1*vy+cs),       beta*g1},
          {(1.0-2.0*beta*phi), 2.0*beta*g1*vx,   2.0*beta*g1*vy,    -2.0*beta*g1},
          {beta*(phi-cs*vy),  -beta*g1*vx,  -beta*(g1*vy-cs),       beta*g1},
          {-vx,                      1.0,            0.0,            0.0}
                    };

    real rx[4][4] = {
          { 1.0,        1.0,        1.0,        0.0},
          {(vx - cs),   vx,     (vx + cs),      0.0},
          {  vy,        vy,         vy,        -1.0},
          {(h - cs*vx), k,      (h + cs*vx),    -vy}
                    };

    real ry[4][4] = {
          { 1.0,        1.0,        1.0,        0.0},
          {vx,   vx,     vx,      1.0},
          {  vy-cs,        vy,         vy+cs,        0.0},
          {(h - cs*vy), k,      (h + cs*vy),    vx}
                    };          
}
*/

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

/*void characterstic_slope_limiter(const real *wij, real *w0l, real *w0r, real *w0b, real *w0t)
{
    real w0[NCONS];
    real w1[NCONS];
    real w2[NCONS];

    // slopes of conserved variables
    real a[NCONS];
    real b[NCONS];    
    real c[NCONS];
    real d[NCONS];    
    real e[NCONS];

    // slopes of characteristic variables
    real ca[NCONS];
    real cb[NCONS];    
    real cc[NCONS];
    real cd[NCONS];    
    real ce[NCONS];
    
    //real prim[NCONS];

    for (int q = 0; q < NCONS; ++q)
    {
        w0[q] = wij[q * NK];
        w1[q] = wij[q * NK + 1]; // y slopes
        w2[q] = wij[q * NK + 2]; // x slopes

        a[q] =  w2[q] / SQRT_THREE;
        b[q] = (w0[q] - w0l[q]) * BETA;
        c[q] = (w0r[q] - w0[q]) * BETA;

        d[q] =  w1[q] / SQRT_THREE;
        e[q] = (w0[q] - w0b[q]) * BETA;
        f[q] = (w0t[q] - w0[q]) * BETA;
    }
    // Use 0th values of the weights, w0, as the conserved variables
    // to calculate primitive variabbles

    conserved_to_primitive(w0, prim);

    const real cs2 = primitive_to_sound_speed_squared(prim);
    real cs  = square_root(cs2);
    const real g1 = ADIABATIC_GAMMA - 1.0;

    real k = 0.5 * (vx * vx + vy * vy);
    real h = (cs2 / g1) + k;
    real phi = g1 * k;
    real beta = 1.0 / (2.0 * cs2);

    // Schaal+ Appendix E

    real lx[4][4] = {
          {beta*(phi+cs*vx),  -beta*(g1*vx+cs),  -beta*g1*vy,       beta*g1},
          {(1.0-2.0*beta*phi), 2.0*beta*g1*vx,   2.0*beta*g1*vy,    -2.0*beta*g1},
          {beta*(phi-cs*vx),  -beta*(g1*vx-cs),  -beta*g1*vy,       beta*g1},
          {vy,                      0.0,            -1.0,            0.0}
                    };

    // compute slopes of characterstic variables
    for (int i = 0, i < NCONS; ++i)
    {
        ca[i] = 0.0;
        cb[i] = 0.0;
        cc[i] = 0.0;
        cd[i] = 0.0;
        ce[i] = 0.0;
        cf[i] = 0.0;

        for (int j = 0, j < NCONS; ++j)
        {
            ca[i] += lx[i][j] * a[j]; // check this
            cb[i] += lx[i][j] * b[j];
            cc[i] += lx[i][j] * c[j];
            cd[i] += lx[i][j] * d[j];
            ce[i] += lx[i][j] * e[j];
            cf[i] += lx[i][j] * f[j];
        }
    }

                    // x slopes
                wtilde[ NK * q + 2] = minmodTVB(wij[ NK * q + 2], wij[ NK * q + 0], wli[ NK * q + 0], wri[ NK * q + 0], dx);
                // y slopes 
                wtilde[ NK * q + 1] = minmodTVB(wij[ NK * q + 1], wij[ NK * q + 0], wlj[ NK * q + 0], wrj[ NK * q + 0], dy);
                
                if (wtilde[ NK * q + 2] != wij[ NK * q + 2])
                {
                    wij[ NK * q + 2] = wtilde[ NK * q + 2];
                    for (int l = 3; l < NK; ++l)
                        {
                            wij[ NK * q + l] = 0.0;
                        }
                }
                if (wtilde[ NK * q + 1] != wij[ NK * q + 1])
                {
                    wij[ NK * q + 1] = wtilde[ NK * q + 1];
                    for (int l = 3; l < NK; ++l)
                        {
                            wij[ NK * q + l] = 0.0;
                        }
                }

}
*/
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

            //real amag = 1e-2;

            //wij[0 * NK] = 1.0 + amag * 1.0 * sin(2.0*PI*x);
            //wij[1 * NK] = 1.0 + amag * (-1.0) * sin(2.0*PI*x);
            //wij[2 * NK] = amag *  (1.0) * sin(2.0*PI*x);
            //wij[3 * NK] = (1.0 / ADIABATIC_GAMMA) / (ADIABATIC_GAMMA -1.0) + amag * 1.5 * sin(2.0*PI*x);
            //if (j==0) printf("%f %e %e\n",x,sin(2.0*PI*x),wij[1 * NK]);
            
            /*real prim[NCONS];
            real cons[NCONS];
            real flux[NCONS];
            real rho, vx, vy, pressure; 

            // For |y| > 0.25, we set Vx = -0.5 and ρ = 1, for |y| ≤ 0.25, Vx = 0.5 and ρ = 2. 
/*
            if ( y < 0.75 * (y1-y0) && y > 0.25*(y1-y0) )
            {
                vx = 0.5;
                rho = 2.0;
            }
            else
            {
                vx = -0.5;
                rho = 1.0;
            }
            vy = 0.02*sin(2.0*PI*x);
            pressure = 2.5;
            
            //rho      = 2.0 + 0.5*sin(2.0*PI*x);
            //vx       = 3.0;
            //vy       = 0.0;
            //pressure = 2.0;

            prim[0] = rho;
            prim[1] = vx;
            prim[2] = vy;
            prim[3] = pressure;

            primitive_to_conserved(prim,cons);
            //primitive_to_flux_vector(prim,flux,0);

            //printf("cons: %f %f %f %f\n", cons[0],cons[1],cons[2],cons[3]);
            //printf("flux %f %f %f %f\n\n", flux[0],flux[1],flux[2],flux[3]);

            wij[0 * NK] = cons[0];
            wij[1 * NK] = cons[1];
            wij[2 * NK] = cons[2];
            wij[3 * NK] = cons[3]; 
*/
            //if (square_root(r2) < 0.125)

            real prim[NCONS];
            real cons[NCONS];
            real rho, vx, vy, pressure;

            if (x < xmid)    
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

            wij[0 * NK] = cons[0];
            wij[1 * NK] = cons[1];
            wij[2 * NK] = cons[2];
            wij[3 * NK] = cons[3]; 
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
    real *dweights;
    real *trouble; //
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
            /*
            //Periodic BC
            if (il == -1)
                il = ni-1;

            if (ir == ni)
                ir = 0;

            if (jl == -1)
                jl = nj-1;

            if (jr == nj)
                jr = 0;
            
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
                
                riemann_hllc(primm, primp, flux_x, 0);
                
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
                
                riemann_hllc(primm, primp, flux_x, 0);
            
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
                
                riemann_hllc(primm, primp, flux_y, 1);

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
                
                riemann_hllc(primm, primp, flux_y, 1);

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

void limit_weights(struct UpdateStruct update)
{
    int ni = update.ni;
    int nj = update.nj;

    const real dx = (update.x1 - update.x0) / update.ni;
    const real dy = (update.y1 - update.y0) / update.nj;
    const real dvol = 4;

    real cons[NCONS];
    real prim[NCONS];

    const real ck = 1e30;       //k=1; see Fu & Shu Sec. 3
    //const real ck = 0.03;        //k=1; see Fu & Shu Sec. 3
    //const real ck = 0.03 * 2 ; //k=2; see Fu & Shu Sec. 3
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
            /*
            //Periodic BC
            if (il == -1)
                il = ni-1;

            if (ir == ni)
                ir = 0;

            if (jl == -1)
                jl = nj-1;

            if (jr == nj)
                jr = 0;
            
            /* */ real *wij = &update.weights[NCONS * NK * (i  * nj + j )];

            const real *wli = &update.weights[NCONS * NK * (il * nj + j )];
            const real *wri = &update.weights[NCONS * NK * (ir * nj + j )];
            const real *wlj = &update.weights[NCONS * NK * (i  * nj + jl)];
            const real *wrj = &update.weights[NCONS * NK * (i  * nj + jr)];

            /* */ real  *tij = &update.trouble[i * nj + j];

            // Troubled Cell Indicator G. Fu & C.-W. Shu (JCP, 347, 305 (2017))

            // Try Second order first

            // eq. 2.3 compute max bar p and double bar p for rho and E
            
            // integral of p over it's own cell

            // l=0 for rho

            real maxpj = maxabs(wij[0], wli[0], wri[0], wlj[0], wrj[0]);

            real pbb_li = fabs(wij[0] - (4.0 * wli[0] + 8.0 * wli[2]) / dvol);
            real pbb_ri = fabs(wij[0] - (4.0 * wri[0] - 8.0 * wri[2]) / dvol);
            real pbb_lj = fabs(wij[0] - (4.0 * wlj[0] + 8.0 * wlj[1]) / dvol);
            real pbb_rj = fabs(wij[0] - (4.0 * wrj[0] - 8.0 * wrj[1]) / dvol);

            real trouble_rho = (pbb_li + pbb_ri + pbb_lj + pbb_rj) / maxpj;

            //printf("%d %d %f %f %f %f\n",i,j,(4.0 * wli[0] + 8.0 * wli[2]),wij[0],wli[0],wli[2]);
            //printf("%d %d %f %f %f %f %f\n",i,j,pbb_li, pbb_ri, pbb_lj, pbb_rj, maxpj);

            // l=3 for energy

            /**/ maxpj = maxabs(wij[3*NK+0], wli[3*NK+0], wri[3*NK+0], wlj[3*NK+0], wrj[3*NK+0]);

            /**/ pbb_li = fabs(wij[3*NK+0] - (4.0 * wli[3*NK+0] + 8.0 * wli[3*NK+2]) / dvol);
            /**/ pbb_ri = fabs(wij[3*NK+0] - (4.0 * wri[3*NK+0] - 8.0 * wri[3*NK+2]) / dvol);
            /**/ pbb_lj = fabs(wij[3*NK+0] - (4.0 * wlj[3*NK+0] + 8.0 * wlj[3*NK+1]) / dvol);
            /**/ pbb_rj = fabs(wij[3*NK+0] - (4.0 * wrj[3*NK+0] - 8.0 * wrj[3*NK+1]) / dvol);

            real trouble_e = (pbb_li + pbb_ri + pbb_lj + pbb_rj) / maxpj;

            //printf("trouble: %d %d %f %f\n",i,j,trouble_rho,trouble_e);

            *tij = trouble_e;

            if ( trouble_rho > ck )
            {
            // Do characteristic slope limiting
            // slope limiting 

                for (int q = 0; q < NCONS; ++q)
                {
                    #if (DG_ORDER == 2)

                    // x slopes
                    wij[ NK * q + 2] = minmod(wij[ NK * q + 2], wli[ NK * q + 0], wij[ NK * q + 0], wri[ NK * q + 0], dx);

                    // y slopes 
                    wij[ NK * q + 1] = minmod(wij[ NK * q + 1], wlj[ NK * q + 0], wij[ NK * q + 0], wrj[ NK * q + 0], dy);

                    #endif 
                }
            }
        }
    }
}

int main()
{
    const int ni = 100;
    const int nj = 100;
    const int fold = 10;
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

    int iteration = 0;
    real time = 0.0;
    real dt = dx * CFL;

    while (time < 0.1)
    {
        clock_t start = clock();

        for (int i = 0; i < fold; ++i)
        {
            compute_delta_weights(update);
            add_delta_weights(update, dt);
            limit_weights(update);

            time += dt;
            iteration += 1;
        }
        clock_t end = clock();

        real seconds = ((real) (end - start)) / CLOCKS_PER_SEC;
        real mcps = (ni * nj / 1e6) / seconds * fold;
        printf("[%d] t=%.3e Mcps=%.2f Mnps=%.2f\n", iteration, time, mcps, NCELL * mcps);
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

            fprintf(outfile, "%f %f %f %f %f %f %f\n", x, y, prim[0], prim[1], prim[2], prim[3], *tij);
        }
    }
    fclose(outfile);
    free(weights_host);
    free(trouble_host);
    return 0;
}

/*
    real w0[NCONS];
    real w1[NCONS];
    real w2[NCONS];

    // slopes of conserved variables
    real a[NCONS];
    real b[NCONS];    
    real c[NCONS];
    real d[NCONS];    
    real e[NCONS];

    // slopes of characteristic variables
    real ca[NCONS];
    real cb[NCONS];    
    real cc[NCONS];
    real cd[NCONS];    
    real ce[NCONS];
    
    //real prim[NCONS];

    for (int q = 0; q < NCONS; ++q)
    {
        w0[q] = wij[q * NK];
        w1[q] = wij[q * NK + 1]; // y slopes
        w2[q] = wij[q * NK + 2]; // x slopes

        a[q] =  w2[q] / SQRT_THREE;
        b[q] = (w0[q] - w0l[q]) * BETA;
        c[q] = (w0r[q] - w0[q]) * BETA;

        d[q] =  w1[q] / SQRT_THREE;
        e[q] = (w0[q] - w0b[q]) * BETA;
        f[q] = (w0t[q] - w0[q]) * BETA;
    }
    // Use 0th values of the weights, w0, as the conserved variables
    // to calculate primitive variabbles

    conserved_to_primitive(w0, prim);

    const real cs2 = primitive_to_sound_speed_squared(prim);
    real cs  = square_root(cs2);
    const real g1 = ADIABATIC_GAMMA - 1.0;

    real k = 0.5 * (vx * vx + vy * vy);
    real h = (cs2 / g1) + k;
    real phi = g1 * k;
    real beta = 1.0 / (2.0 * cs2);

    // Schaal+ Appendix E

    real lx[4][4] = {
          {beta*(phi+cs*vx),  -beta*(g1*vx+cs),  -beta*g1*vy,       beta*g1},
          {(1.0-2.0*beta*phi), 2.0*beta*g1*vx,   2.0*beta*g1*vy,    -2.0*beta*g1},
          {beta*(phi-cs*vx),  -beta*(g1*vx-cs),  -beta*g1*vy,       beta*g1},
          {vy,                      0.0,            -1.0,            0.0}
                    };

    // compute slopes of characterstic variables
    for (int i = 0, i < NCONS; ++i)
    {
        ca[i] = 0.0;
        cb[i] = 0.0;
        cc[i] = 0.0;
        cd[i] = 0.0;
        ce[i] = 0.0;
        cf[i] = 0.0;

        for (int j = 0, j < NCONS; ++j)
        {
            ca[i] += lx[i][j] * a[j]; // check this
            cb[i] += lx[i][j] * b[j];
            cc[i] += lx[i][j] * c[j];
            cd[i] += lx[i][j] * d[j];
            ce[i] += lx[i][j] * e[j];
            cf[i] += lx[i][j] * f[j];
        }
    }
    // Schaal+ Appendix E

    real lx[4][4] = {
          {beta*(phi+cs*vx),  -beta*(g1*vx+cs),  -beta*g1*vy,       beta*g1},
          {(1.0-2.0*beta*phi), 2.0*beta*g1*vx,   2.0*beta*g1*vy,    -2.0*beta*g1},
          {beta*(phi-cs*vx),  -beta*(g1*vx-cs),  -beta*g1*vy,       beta*g1},
          {vy,                      0.0,            -1.0,            0.0}
                    };

    real ly[4][4] = {
          {beta*(phi+cs*vy),  -beta*g1*vx,  -beta*(g1*vy+cs),       beta*g1},
          {(1.0-2.0*beta*phi), 2.0*beta*g1*vx,   2.0*beta*g1*vy,    -2.0*beta*g1},
          {beta*(phi-cs*vy),  -beta*g1*vx,  -beta*(g1*vy-cs),       beta*g1},
          {-vx,                      1.0,            0.0,            0.0}
                    };

    real rx[4][4] = {
          { 1.0,        1.0,        1.0,        0.0},
          {(vx - cs),   vx,     (vx + cs),      0.0},
          {  vy,        vy,         vy,        -1.0},
          {(h - cs*vx), k,      (h + cs*vx),    -vy}
                    };

    real ry[4][4] = {
          { 1.0,        1.0,        1.0,        0.0},
          {vx,   vx,     vx,      1.0},
          {  vy-cs,        vy,         vy+cs,        0.0},
          {(h - cs*vy), k,      (h + cs*vy),    vx}
                    };  
            // compute LeftEigenX and LeftEigenY matrices from wij[ NK * q + 0]

            // compute slopes of characteristic variables:
            //  c1[q] = LeftEigenX(wij[ NK * q + 1])
            //  c2[q] = LeftEigenX(wij[ NK * q + 2])
            //  left[q] = LeftEigenX(wij[NK * q + 0] - wli[NK * q + 0])
            //  right[q] = LeftEigenX(wri[NK * q + 0] - wij[NK * q + 0])
            //  bottom[q] = LeftEigenY(wij[NK * q + 0] - wlj[NK * q + 0])
            //  right[q] = LeftEigenY(wrj[NK * q + 0] - wij[NK * q + 0])

            // use minmodTVB to compute limited characteristic slopes
            // c1tilde[q] and c2tilde[q] 

            // if c1tilde[q] != c1[q] OR c2tilde[q] != c2[q]
            // w1tilde[q] = RightEigenX(c1tilde[q])
            // w2tilde[q] = RightEigenX(c2tilde[q])
            // set
            // wij[ NK * q + 1] = w1tilde[q]
            // wij[ NK * q + 2] = w2tilde[q]
            // wij[ NK * q + (l>2)] = 0.0

            // slope limiting 
            for (int q = 0; q < NCONS; ++q)
            {
                #if (DG_ORDER == 2)

                // x slopes
                wij[ NK * q + 2] = minmodTVB(wij[ NK * q + 2], wli[ NK * q + 0], wij[ NK * q + 0], wri[ NK * q + 0], dx);

                // y slopes 
                wij[ NK * q + 1] = minmodTVB(wij[ NK * q + 1], wlj[ NK * q + 0], wij[ NK * q + 0], wrj[ NK * q + 0], dy);

                #elif (DG_ORDER >= 2)

                // x slopes
                wtilde[ NK * q + 2] = minmodTVB(wij[ NK * q + 2], wij[ NK * q + 0], wli[ NK * q + 0], wri[ NK * q + 0], dx);
                // y slopes 
                wtilde[ NK * q + 1] = minmodTVB(wij[ NK * q + 1], wij[ NK * q + 0], wlj[ NK * q + 0], wrj[ NK * q + 0], dy);
                
                if (wtilde[ NK * q + 2] != wij[ NK * q + 2])
                {
                    wij[ NK * q + 2] = wtilde[ NK * q + 2];
                    for (int l = 3; l < NK; ++l)
                        {
                            wij[ NK * q + l] = 0.0;
                        }
                }
                if (wtilde[ NK * q + 1] != wij[ NK * q + 1])
                {
                    wij[ NK * q + 1] = wtilde[ NK * q + 1];
                    for (int l = 3; l < NK; ++l)
                        {
                            wij[ NK * q + l] = 0.0;
                        }
                }
                #endif
            }
        }
    } 
    */
