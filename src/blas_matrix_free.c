/**
 * @file blas_matrix_free.c
 * @brief BLAS lib
 * @author Ziteng Wang
 * @version 1.0
 * @date 2013-03-03
 */
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "multigrid_functs.h"


/**
 * @brief Compute Residue Vector 1D
 *
 * @param b right hand
 * @param u solution
 * @param r residue vector
 * @param k current level of grids
 * @param level indicator of each level
 */
void compute_r_1d(REAL *b,
                  REAL *u,
                  REAL *r,
                  INT   k,
                  INT *level)
{
    INT i,n;

    n = level[k+1]-level[k];

    // r[level[k]] = b[level[k]]-2*u[level[k]]+u[level[k]+1];
    for (i = 1; i < n-1; i++) {
        r[level[k]+i] = b[level[k]+i]-2*u[level[k]+i]+u[level[k]+i+1]+u[level[k]+i-1];
    }
    r[level[k]+n-1] = b[level[k]+n-1]-2*u[level[k]+n-1]+u[level[k]+n-2];
}

/**
 * @brief Compute Residue Vector 2D
 *
 * @param u right hand
 * @param b solution
 * @param r residue vector
 * @param k current level of grids
 * @param level indicator of each level
 * @param nxk number of cells in x direction in level k
 * @param nyk number of cells in y direction in level k
 */
void compute_r_2d(REAL *u,
                  REAL *b,
                  REAL *r,
                  INT k,
                  INT *level,
                  INT *nxk,
                  INT *nyk)
{
    INT i,j;

    //r[level[k]] = b[level[k]]-4*u[level[k]]+u[level[k]+1]+u[level[k]+nxk[k]];
    //for (j = 1; j < nxk[k]-1; j++) {
    //    r[level[k]+j] = b[level[k]+j]-4*u[level[k]+j]+u[level[k]+j+1]+u[level[k]+j-1]+u[level[k]+j+nxk[k]];       
    //}
    //r[level[k]+nxk[k]-1] = b[level[k]+nxk[k]-1]-4*u[level[k]+nxk[k]-1]+u[level[k]+nxk[k]-2]+
    //                       u[level[k]+2*nxk[k]-1];
    for (i = 1; i < nyk[k]-1; i++) {
//        r[level[k]+i*nxk[k]] = b[level[k]+i*nxk[k]]-4*u[level[k]+i*nxk[k]]+u[level[k]+i*nxk[k]+1]+
//                               u[level[k]+(i+1)*nxk[k]]+u[level[k]+(i-1)*nxk[k]];
        for (j = 1; j < nxk[k]-1; j++) {
            r[level[k]+i*nxk[k]+j] = b[level[k]+i*nxk[k]+j]-4*u[level[k]+i*nxk[k]+j]+u[level[k]+i*nxk[k]+j+1]+
                                     u[level[k]+i*nxk[k]+j-1]+u[level[k]+(i+1)*nxk[k]+j]+
                                     u[level[k]+(i-1)*nxk[k]+j];
        }
//        r[level[k]+(i+1)*nxk[k]-1] = b[level[k]+(i+1)*nxk[k]-1]-4*u[level[k]+(i+1)*nxk[k]-1]+
//                                     u[level[k]+(i+1)*nxk[k]-2]+u[level[k]+(i+2)*nxk[k]-1]+
//                                     u[level[k]+(i)*nxk[k]-1];
    }
    //r[level[k]+(nyk[k]-1)*nxk[k]] = b[level[k]+(nyk[k]-1)*nxk[k]]-4*u[level[k]+(nyk[k]-1)*nxk[k]]+
    //                                u[level[k]+(nyk[k]-1)*nxk[k]+1]+u[level[k]+(nyk[k]-2)*nxk[k]];
    //for (j = 1; j < nxk[k]-1; j++) {
    //    r[level[k]+(nyk[k]-1)*nxk[k]+j] = b[level[k]+(nyk[k]-1)*nxk[k]+j]-4*u[level[k]+(nyk[k]-1)*nxk[k]+j]+
    //                                      u[level[k]+(nyk[k]-1)*nxk[k]+j+1]+u[level[k]+(nyk[k]-1)*nxk[k]+j-1]+
    //                                      u[level[k]+(nyk[k]-2)*nxk[k]+j];
    //}
    //r[level[k]+nyk[k]*nxk[k]-1] = b[level[k]+nyk[k]*nxk[k]-1]-4*u[level[k]+nyk[k]*nxk[k]-1]+
    //                              u[level[k]+nyk[k]*nxk[k]-2]+u[level[k]+(nyk[k]-1)*nxk[k]-1];
}    
            

/**
 * @brief Compute Residue Vector 3D
 *
 * @param u right hand vector
 * @param b solution
 * @param r residue vector
 * @param k current level of the multigrids
 * @param level indicator of each level k
 * @param nxk number of cells in x direction in level k
 * @param nyk number of cells in y direction in level k
 * @param nzk number of cells in z direction in level k
 */
void compute_r_3d(REAL *u,
                  REAL *b,
                  REAL *r,
                  INT k,
                  INT *level,
                  INT *nxk,
                  INT *nyk,
                  INT *nzk)
{
    int i,j,h;       
	int i0,j0,j1,j2,k0,k1,k2,k3,k4,k5,k6;
    const int levelk = level[k];
	const int nxkk = nxk[k];
	const int nykk = nyk[k];
	const int nzkk = nzk[k];

    // middle part of the cubic
    for (i = 1; i < nzkk-1; i++) {
		i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j++) {
		    j0 = i0+j*nxkk;
			j1 = i0+(j+1)*nxkk;
			j2 = i0+(j-1)*nxkk;
			for (h = 1; h < nxkk-1; h++) {
/*                r[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h] = b[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]-
                                                         6*u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
                                                         u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1]+
                                                         u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1]+
                                                         u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h]+
                                                         u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h]+
                                                         u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h]+
                                                         u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]+h];
 */
				k0 = j0+h;
				k1 = j0+h-1;
				k2 = j0+h+1;
				k3 = j1+h;
				k4 = j2+h;
				k5 = k0+nxkk*nykk;
				k6 = k0-nxkk*nykk;
				r[k0]=b[k0]-6*u[k0]+u[k1]+u[k2]+u[k3]+u[k4]+u[k5]+u[k6];
			}
        }
	}
}    


/**
 * @brief Compute L2 Norm
 *
 * @param r residue vector
 * @param level indicator of each level k
 * @param k current level of the multigrid
 *
 * @return L2 Norm
 */
REAL computenorm(REAL *r,
                 INT *level,
                 INT k)
{  
    INT i,n;
    REAL squarnorm;

    squarnorm = 0.0;
    n = level[k+1]-level[k];
    for (i = 1; i < n; i++) {
		//printf("%f\n",u[i]);
        squarnorm = squarnorm + r[level[k]+i]*r[level[k]+i];    
    }
    squarnorm = sqrt(squarnorm);
    return squarnorm;
} 

/**
 * @brief x = y
 *
 * @param x vector x
 * @param y vector y
 * @param level indicator of each level of multigrids
 * @param k current level of multigrids
 */
void xequaly(REAL *x,
             REAL *y,
             INT *level,
             INT k)
{
    INT i,n;

    n = level[k+1] - level[k];

    for (i = 0; i < n; i++) {
        x[level[k]+i] = y[level[k]+i];
    }
}

/**
 * @brief x = y + c*z
 *
 * @param x vector x
 * @param y vector y
 * @param c coefficient c
 * @param z vector z
 * @param level indicator of each level of multigrids
 * @param k current level of multigrids
 */
void xequalypcz(REAL *x,
                REAL *y,
                REAL c,
                REAL *z,
                INT *level,
                INT k)
{
    INT i, n;

    n = level[k+1]-level[k];

    for (i = 0; i < n; i++) {
        x[level[k]+i] = y[level[k]+i] + c*z[level[k]+i];            
    }
}

/**
 * @brief Compute quadratic form of A 1D
 *
 * @param u vector u
 * @param level indicator of each level k of multigrids
 * @param k current level of multigrids
 * @param nxk number of cells in x direction
 *
 * @return quadratic form of A 
 */
REAL compute_quad_norm_1d(REAL *u,
                          INT *level, 
                          INT k,
                          INT *nxk)
{
    INT i,n;
    REAL norm;

    n = level[k+1] - level[k];
    norm = 0.0;

    norm += 2*u[level[k]]*u[level[k]] - u[level[k]]*u[level[k]+1];
    for (i = 1; i < n-1; i++) {
        norm += 2*u[level[k]+i]*u[level[k]+i] - u[level[k]+i]*u[level[k]+i-1] - u[level[k]+i+1]*u[level[k]+i];            
    }
    norm += 2*u[level[k]+n-1]*u[level[k]+n-1] - u[level[k]+n-1]*u[level[k]+n-2];
    
    return norm;
}

/**
 * @brief Compute Quadratic form of A 2D
 *
 * @param u vector u
 * @param level indicator of each level of multigrids
 * @param k current level of multigrids
 * @param nxk number of cells in x direction in level k
 * @param nyk number of cells in y direction in level k
 *
 * @return norm
 */
REAL compute_quad_norm_2d(REAL *u,
                          INT *level,
                          INT k,
                          INT *nxk,
                          INT *nyk)
{
    INT i,j;
    REAL norm;

    norm = 0.0;

    // first line
    norm += u[level[k]]*u[level[k]]*4-u[level[k]]*u[level[k]+1]-u[level[k]]*u[level[k]+nxk[k]];
    for (i = 1; i < nxk[k]-1; i++) {
        norm += u[level[k]+i]*u[level[k]+i]*4-u[level[k]+i-1]*u[level[k]+i]-
                u[level[k]+i]*u[level[k]+i+1]-u[level[k]+i]*u[level[k]+nxk[k]+i];
    }
    norm += u[level[k]+nxk[k]-1]*u[level[k]+nxk[k]-1]*4-u[level[k]+nxk[k]-1]*u[level[k]+nxk[k]-2]-
            u[level[k]+nxk[k]-1]*u[level[k]+2*nxk[k]-1];
    // middle part
    for (i = 1; i < nyk[k]-1; i++) {
        norm += u[level[k]+i*nxk[k]]*u[level[k]+i*nxk[k]]*4-u[level[k]+i*nxk[k]]*u[level[k]+i*nxk[k]+1]-
                u[level[k]+i*nxk[k]]*u[level[k]+(i+1)*nxk[k]]-
                u[level[k]+i*nxk[k]]*u[level[k]+(i-1)*nxk[k]];
        for (j = 1; j < nxk[k]-1; j++) {
            norm += u[level[k]+i*nxk[k]+j]*u[level[k]+i*nxk[k]+j]*4-
                    u[level[k]+i*nxk[k]+j]*u[level[k]+i*nxk[k]+j+1]-
                    u[level[k]+i*nxk[k]+j]*u[level[k]+i*nxk[k]+j-1]-
                    u[level[k]+i*nxk[k]+j]*u[level[k]+(i+1)*nxk[k]+j]-
                    u[level[k]+i*nxk[k]+j]*u[level[k]+(i-1)*nxk[k]+j];
        }
        norm += u[level[k]+(i+1)*nxk[k]-1]*u[level[k]+(i+1)*nxk[k]-1]*4-
                u[level[k]+(i+1)*nxk[k]-1]*u[level[k]+(i+1)*nxk[k]-2]-
                u[level[k]+(i+1)*nxk[k]-1]*u[level[k]+(i+2)*nxk[k]-1]-
                u[level[k]+(i+1)*nxk[k]-1]*u[level[k]+(i)*nxk[k]-1];
    }
    // last line
    norm += u[level[k]+(nyk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]-1)*nxk[k]]*4-
            u[level[k]+(nyk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]-1)*nxk[k]+1]-
            u[level[k]+(nyk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]-2)*nxk[k]];
    for (i = 1; i < nxk[k]-1; i++) {
        norm += u[level[k]+(nyk[k]-1)*nxk[k]+i]*u[level[k]+(nyk[k]-1)*nxk[k]+i]*4-
                u[level[k]+(nyk[k]-1)*nxk[k]+i]*u[level[k]+(nyk[k]-1)*nxk[k]+i+1]-
                u[level[k]+(nyk[k]-1)*nxk[k]+i]*u[level[k]+(nyk[k]-1)*nxk[k]+i-1]-
                u[level[k]+(nyk[k]-1)*nxk[k]+i]*u[level[k]+(nyk[k]-2)*nxk[k]+i];
    }
    norm += u[level[k]+nxk[k]*nyk[k]-1]*u[level[k]+nxk[k]*nyk[k]-1]*4-
            u[level[k]+nxk[k]*nyk[k]-1]*u[level[k]+nxk[k]*nyk[k]-2]-
            u[level[k]+nxk[k]*nyk[k]-1]*u[level[k]+nxk[k]*nyk[k]-1-nxk[k]];

    return norm;
}

/**
 * @brief Compute Quadratic form of A 3D
 *
 * @param u vector u
 * @param level indicator of each level of multigrids
 * @param k current level of multigrids
 * @param nxk number of cells in x direction in level k
 * @param nyk number of cells in y direction in level k
 * @param nzk number of cells in z direction in level k
 *
 * @return norm
 */
REAL compute_quad_norm_3d(REAL *u,
                          INT *level,
                          INT k,
                          INT *nxk,
                          INT *nyk,
                          INT *nzk)
{
    INT i,j,h;
    REAL norm;

    norm = 0.0;
    // bottom plane
    // first line
    norm += u[level[k]]*u[level[k]]*6-u[level[k]]*u[level[k]+1]-u[level[k]]*u[level[k]+nxk[k]]-
            u[level[k]]*u[level[k]+nxk[k]*nyk[k]];
    for (h = 1; h < nxk[k]-1; h++) {
        norm += u[level[k]+h]*u[level[k]+h]*6-u[level[k]+h]*u[level[k]+h+1]-
                u[level[k]+h]*u[level[k]+h-1]-u[level[k]+h]*u[level[k]+h+nxk[k]]-
                u[level[k]+h]*u[level[k]+h+nxk[k]*nyk[k]];
    }
    norm += u[level[k]+nxk[k]-1]*u[level[k]+nxk[k]-1]*6-u[level[k]+nxk[k]-1]*u[level[k]+nxk[k]-2]-
            u[level[k]+nxk[k]-1]*u[level[k]+2*nxk[k]-1]-
            u[level[k]+nxk[k]-1]*u[level[k]+nxk[k]-1+nxk[k]*nyk[k]];
    // middle parts
    for (j = 1; j < nyk[k]-1; j++) {
        norm += u[level[k]+j*nxk[k]]*u[level[k]+j*nxk[k]]*6-u[level[k]+j*nxk[k]]*u[level[k]+j*nxk[k]+1]-
                u[level[k]+j*nxk[k]]*u[level[k]+(j-1)*nxk[k]]-
                u[level[k]+j*nxk[k]]*u[level[k]+(j+1)*nxk[k]]-
                u[level[k]+j*nxk[k]]*u[level[k]+j*nxk[k]+nxk[k]*nyk[k]];
        for (h = 1; h < nxk[k]-1; h++) {
            norm += u[level[k]+j*nxk[k]+h]*u[level[k]+j*nxk[k]+h]*6-
                    u[level[k]+j*nxk[k]+h]*u[level[k]+j*nxk[k]+h+1]-
                    u[level[k]+j*nxk[k]+h]*u[level[k]+j*nxk[k]+h-1]-
                    u[level[k]+j*nxk[k]+h]*u[level[k]+(j+1)*nxk[k]+h]-
                    u[level[k]+j*nxk[k]+h]*u[level[k]+(j-1)*nxk[k]+h]-
                    u[level[k]+j*nxk[k]+h]*u[level[k]+j*nxk[k]+h+nxk[k]*nyk[k]];
        }
        norm += u[level[k]+(j+1)*nxk[k]-1]*u[level[k]+(j+1)*nxk[k]-1]*6-
                u[level[k]+(j+1)*nxk[k]-1]*u[level[k]+(j+1)*nxk[k]-2]-
                u[level[k]+(j+1)*nxk[k]-1]*u[level[k]+(j+2)*nxk[k]-1]-
                u[level[k]+(j+1)*nxk[k]-1]*u[level[k]+(j)*nxk[k]-1]-
                u[level[k]+(j+1)*nxk[k]-1]*u[level[k]+(j+1)*nxk[k]-1+nxk[k]*nyk[k]];
    }
    // last line
    norm += u[level[k]+(nyk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]-1)*nxk[k]]*6-
            u[level[k]+(nyk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]-1)*nxk[k]+1]-
            u[level[k]+(nyk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]-2)*nxk[k]]-
            u[level[k]+(nyk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]-1)*nxk[k]+nxk[k]*nyk[k]];
    for (h = 1; h < nxk[k]-1; h++) {
        norm += u[level[k]+(nyk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]-1)*nxk[k]+h]*6-
                u[level[k]+(nyk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]-1)*nxk[k]+h+1]-
                u[level[k]+(nyk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]-1)*nxk[k]+h-1]-
                u[level[k]+(nyk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]-2)*nxk[k]+h]-
                u[level[k]+(nyk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]-1)*nxk[k]+h+nxk[k]*nyk[k]];
    }
    norm += u[level[k]+nxk[k]*nyk[k]-1]*u[level[k]+nxk[k]*nyk[k]-1]*6-
            u[level[k]+nxk[k]*nyk[k]-1]*u[level[k]+nxk[k]*nyk[k]-2]-
            u[level[k]+nxk[k]*nyk[k]-1]*u[level[k]+nxk[k]*nyk[k]-1-nxk[k]]-
            u[level[k]+nxk[k]*nyk[k]-1]*u[level[k]+nxk[k]*nyk[k]-1+nxk[k]*nyk[k]];

    // middle part of the x-y-z cubic
    for (i = 1; i < nzk[k]-1; i++) {
        // first line
        norm += u[level[k]+i*nxk[k]*nyk[k]]*u[level[k]+i*nxk[k]*nyk[k]]*6-
                u[level[k]+i*nxk[k]*nyk[k]]*u[level[k]+i*nxk[k]*nyk[k]+1]-
                u[level[k]+i*nxk[k]*nyk[k]]*u[level[k]+i*nxk[k]*nyk[k]+nxk[k]]-
                u[level[k]+i*nxk[k]*nyk[k]]*u[level[k]+(i+1)*nxk[k]*nyk[k]]-
                u[level[k]+i*nxk[k]*nyk[k]]*u[level[k]+(i-1)*nxk[k]*nyk[k]];
        for (h = 1; h < nxk[k]-1; h++) {
            norm += u[level[k]+i*nxk[k]*nyk[k]+h]*u[level[k]+i*nxk[k]*nyk[k]+h]*6-
                    u[level[k]+i*nxk[k]*nyk[k]+h]*u[level[k]+i*nxk[k]*nyk[k]+h+1]-
                    u[level[k]+i*nxk[k]*nyk[k]+h]*u[level[k]+i*nxk[k]*nyk[k]+h-1]-
                    u[level[k]+i*nxk[k]*nyk[k]+h]*u[level[k]+i*nxk[k]*nyk[k]+h+nxk[k]]-
                    u[level[k]+i*nxk[k]*nyk[k]+h]*u[level[k]+(i+1)*nxk[k]*nyk[k]+h]-
                    u[level[k]+i*nxk[k]*nyk[k]+h]*u[level[k]+(i-1)*nxk[k]*nyk[k]+h];
        }
        norm += u[level[k]+i*nxk[k]*nyk[k]+nxk[k]-1]*u[level[k]+i*nxk[k]*nyk[k]+nxk[k]-1]*6-
                u[level[k]+i*nxk[k]*nyk[k]+nxk[k]-1]*u[level[k]+i*nxk[k]*nyk[k]+nxk[k]-2]-
                u[level[k]+i*nxk[k]*nyk[k]+nxk[k]-1]*u[level[k]+i*nxk[k]*nyk[k]+2*nxk[k]-1]-
                u[level[k]+i*nxk[k]*nyk[k]+nxk[k]-1]*u[level[k]+(i+1)*nxk[k]*nyk[k]+nxk[k]-1]-
                u[level[k]+i*nxk[k]*nyk[k]+nxk[k]-1]*u[level[k]+(i-1)*nxk[k]*nyk[k]+nxk[k]-1];
        // middle lines
        for (h = 1; j < nyk[k]-1; j++) {
            norm += u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]]*u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]]*6-
                    u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]]*u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+1]-
                    u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]]*u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]]-
                    u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]]*u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]]-
                    u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]]*u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]]-
                    u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]]*u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]];
            for (h = 1; h < nxk[k]-1; h++) {
                norm += u[level[k]+(i*nyk[k]+j)*nxk[k]+h]*u[level[k]+(i*nyk[k]+j)*nxk[k]+h]*6-
                        u[level[k]+(i*nyk[k]+j)*nxk[k]+h]*u[level[k]+(i*nyk[k]+j)*nxk[k]+h+1]-
                        u[level[k]+(i*nyk[k]+j)*nxk[k]+h]*u[level[k]+(i*nyk[k]+j)*nxk[k]+h-1]-
                        u[level[k]+(i*nyk[k]+j)*nxk[k]+h]*u[level[k]+(i*nyk[k]+j-1)*nxk[k]+h]-
                        u[level[k]+(i*nyk[k]+j)*nxk[k]+h]*u[level[k]+(i*nyk[k]+j+1)*nxk[k]+h]-
                        u[level[k]+(i*nyk[k]+j)*nxk[k]+h]*u[level[k]+((i+1)*nyk[k]+j)*nxk[k]+h]-
                        u[level[k]+(i*nyk[k]+j)*nxk[k]+h]*u[level[k]+((i-1)*nyk[k]+j)*nxk[k]+h];
            }
            norm += u[level[k]+(i*nxk[k]+j+1)*nxk[k]-1]*u[level[k]+(i*nxk[k]+j+1)*nxk[k]-1]*6-
                    u[level[k]+(i*nxk[k]+j+1)*nxk[k]-1]*u[level[k]+(i*nxk[k]+j+1)*nxk[k]-2]-
                    u[level[k]+(i*nxk[k]+j+1)*nxk[k]-1]*u[level[k]+(i*nxk[k]+j)*nxk[k]-1]-
                    u[level[k]+(i*nxk[k]+j+1)*nxk[k]-1]*u[level[k]+(i*nxk[k]+j+2)*nxk[k]-1]-
                    u[level[k]+(i*nxk[k]+j+1)*nxk[k]-1]*u[level[k]+((i+1)*nxk[k]+j+1)*nxk[k]-1]-
                    u[level[k]+(i*nxk[k]+j+1)*nxk[k]-1]*u[level[k]+((i-1)*nxk[k]+j+1)*nxk[k]-1];
        }
        // last line
        norm += u[level[k]+((i+1)*nyk[k]-1)*nxk[k]]*u[level[k]+((i+1)*nyk[k]-1)*nxk[k]]*6-
                u[level[k]+((i+1)*nyk[k]-1)*nxk[k]]*u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+1]-
                u[level[k]+((i+1)*nyk[k]-1)*nxk[k]]*u[level[k]+((i+1)*nyk[k]-2)*nxk[k]]-
                u[level[k]+((i+1)*nyk[k]-1)*nxk[k]]*u[level[k]+((i)*nyk[k]-1)*nxk[k]]-
                u[level[k]+((i+1)*nyk[k]-1)*nxk[k]]*u[level[k]+((i+2)*nyk[k]-1)*nxk[k]];
        for (h = 1; h < nxk[k]-1; h++) {
            norm += u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+h]*u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+h]*6-
                    u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+h]*u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+h+1]-
                    u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+h]*u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+h-1]-
                    u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+h]*u[level[k]+((i+1)*nyk[k]-2)*nxk[k]+h]-
                    u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+h]*u[level[k]+((i)*nyk[k]-1)*nxk[k]+h]-
                    u[level[k]+((i+1)*nyk[k]-1)*nxk[k]+h]*u[level[k]+((i+2)*nyk[k]-1)*nxk[k]+h];
        }
        norm += u[level[k]+(i+1)*nxk[k]*nyk[k]-1]*u[level[k]+(i+1)*nxk[k]*nyk[k]-1]*6-
                u[level[k]+(i+1)*nxk[k]*nyk[k]-1]*u[level[k]+(i+1)*nxk[k]*nyk[k]-2]-
                u[level[k]+(i+1)*nxk[k]*nyk[k]-1]*u[level[k]+(i+1)*nxk[k]*nyk[k]-1-nxk[k]]-
                u[level[k]+(i+1)*nxk[k]*nyk[k]-1]*u[level[k]+(i+2)*nxk[k]*nyk[k]-1]-
                u[level[k]+(i+1)*nxk[k]*nyk[k]-1]*u[level[k]+(i)*nxk[k]*nyk[k]-1];
    }

    // top plane
    // first line
    norm += u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]]*u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]]*6-
            u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]]*u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+1]-
            u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]]*u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+nxk[k]]-
            u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]]*u[level[k]+(nzk[k]-2)*nxk[k]*nyk[k]];
    for (h = 1; h < nxk[k]-1; h++) {
        norm += u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+h]*u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+h]*6-
                u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+h]*u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+h-1]-
                u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+h]*u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+h+1]-
                u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+h]*u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+h+nxk[k]]-
                u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+h]*u[level[k]+(nzk[k]-2)*nxk[k]*nyk[k]+h];
    }
    norm += 6*u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+nxk[k]-1]*
              u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+nxk[k]-1]-
              u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+nxk[k]-1]*
              u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+nxk[k]-2]-
              u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+nxk[k]-1]*
              u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+2*nxk[k]-1]-
              u[level[k]+(nzk[k]-1)*nxk[k]*nyk[k]+nxk[k]-1]*
              u[level[k]+(nzk[k]-2)*nxk[k]*nyk[k]+nxk[k]-1];
    // middle lines
    for (j = 1; j < nyk[k]-1; j++) {
        norm += u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]]*u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]]*6-
                u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]]*u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+1]-
                u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]]*u[level[k]+((nzk[k]-1)*nyk[k]+j+1)*nxk[k]]-
                u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]]*u[level[k]+((nzk[k]-1)*nyk[k]+j-1)*nxk[k]]-
                u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]]*u[level[k]+((nzk[k]-2)*nyk[k]+j)*nxk[k]];
        for (h = 1; h < nxk[k]-1; h++) {
            norm += 6*u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+h]*
                      u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+h]-
                      u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+h]*
                      u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+h+1]-
                      u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+h]*
                      u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+h-1]-
                      u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+h]*
                      u[level[k]+((nzk[k]-1)*nyk[k]+j+1)*nxk[k]+h]-
                      u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+h]*
                      u[level[k]+((nzk[k]-1)*nyk[k]+j-1)*nxk[k]+h]-
                      u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]+h]*
                      u[level[k]+((nzk[k]-2)*nyk[k]+j)*nxk[k]+h];
        }
        norm += 6*u[level[k]+((nzk[k]-1)*nyk[k]+j+1)*nxk[k]-1]*
                  u[level[k]+((nzk[k]-1)*nyk[k]+j+1)*nxk[k]-1]-
                  u[level[k]+((nzk[k]-1)*nyk[k]+j+1)*nxk[k]-1]*
                  u[level[k]+((nzk[k]-1)*nyk[k]+j+1)*nxk[k]-2]-
                  u[level[k]+((nzk[k]-1)*nyk[k]+j+1)*nxk[k]-1]*
                  u[level[k]+((nzk[k]-1)*nyk[k]+j+2)*nxk[k]-1]-
                  u[level[k]+((nzk[k]-1)*nyk[k]+j+1)*nxk[k]-1]*
                  u[level[k]+((nzk[k]-1)*nyk[k]+j)*nxk[k]-1]-
                  u[level[k]+((nzk[k]-1)*nyk[k]+j+1)*nxk[k]-1]*
                  u[level[k]+((nzk[k]-2)*nyk[k]+j+1)*nxk[k]-1];
    }
    // last line
    norm += u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]]*6-
            u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]+1]-
            u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]*nzk[k]-2)*nxk[k]]-
            u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]]*u[level[k]+(nyk[k]*nzk[k]-1-nyk[k])*nxk[k]];
    for (h = 1; h < nxk[k]-1; h++) {
        norm += u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]+h]*6-
                u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]+h+1]-
                u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]+h-1]-
                u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]*nzk[k]-2)*nxk[k]+h]-
                u[level[k]+(nyk[k]*nzk[k]-1)*nxk[k]+h]*u[level[k]+(nyk[k]*nzk[k]-1-nyk[k])*nxk[k]+h];
    }
    norm += u[level[k]+nyk[k]*nzk[k]*nxk[k]-1]*u[level[k]+nyk[k]*nzk[k]*nxk[k]-1]*6-
            u[level[k]+nyk[k]*nzk[k]*nxk[k]-1]*u[level[k]+nyk[k]*nzk[k]*nxk[k]-2]-
            u[level[k]+nyk[k]*nzk[k]*nxk[k]-1]*u[level[k]+nyk[k]*nzk[k]*nxk[k]-1-nxk[k]]-
            u[level[k]+nyk[k]*nzk[k]*nxk[k]-1]*u[level[k]+nyk[k]*(nzk[k]-1)*nxk[k]-1];

    return norm;
}

/**
 * @brief x = Ay
 *
 * @param x vector x
 * @param y vector y
 * @param level  indicator of each level of multigrids
 * @param k current level of multigrids
 */
void xequalay_1d(REAL *x,
                 REAL *y,
                 INT *level,
                 INT k)
{
    INT i,n;
    REAL *btemp;
    btemp = (REAL *)malloc(level[k+1]*sizeof(REAL));
    for (i = 0; i < level[k+1]; i++) {
        btemp[i] = 0.0;
    }
    // compute (-x)
    compute_r_1d(btemp, y, x, k, level);
    n = level[k+1]-level[k];
    for (i = 0; i < n; i++) {
        x[level[k]+i] = (-1)*x[level[k]+i];
    }
    free(btemp);
}


/**
 * @brief x = Ay 2D
 *
 * @param x vector x
 * @param y vector y
 * @param level  indicator of each level of multigrids
 * @param k current level of multigrid
 * @param nxk number of cells in x direction in level k
 * @param nyk number of cells in y direction in level k
 */
void xequalay_2d(REAL *x,
                 REAL *y,
                 INT *level,
                 INT k,
                 INT *nxk,
                 INT *nyk)
{
    INT i,n;
    REAL *btemp;
    btemp = (REAL *)malloc(level[k+1]*sizeof(REAL));
    for (i = 0; i < level[k+1]; i++) {
        btemp[i] = 0.0;
    }
    // compute (-x)
    compute_r_2d(y, btemp, x, k, level, nxk, nyk);
    n = level[k+1]-level[k];
    for (i = 0; i < n; i++) {
        x[level[k]+i] = (-1)*x[level[k]+i];
    }
    free(btemp);
}

/**
 * @brief x = Ay 3D
 *
 * @param x vector x
 * @param y vector y
 * @param level indicator of each level of multigrids
 * @param k current level of multigrids
 * @param nxk number of cells in x direction in level k
 * @param nyk number of cells in y direction in level k
 * @param nzk number of cells in z direction in level k
 */
void xequalay_3d(REAL *x,
                 REAL *y,
                 INT *level,
                 INT k,
                 INT *nxk,
                 INT *nyk,
                 INT *nzk)
{
    INT i,n;
    REAL *btemp;
    btemp = (REAL *)malloc(level[k+1]*sizeof(REAL));
    for (i = 0; i < level[k+1]; i++) {
        btemp[i] = 0.0;
    }
    // compute (-x)
    compute_r_3d(y, btemp, x, k, level, nxk, nyk, nzk);
    n = level[k+1]-level[k];
    for (i = 0; i < n; i++) {
        x[level[k]+i] = (-1)*x[level[k]+i];
    }
    //free(btemp);
}


/**
 * @brief <x,y>
 *
 * @param x vector x
 * @param y vector y
 * @param level indicator of each level of multigrids
 * @param k current level of multigrids
 *
 * @return 
 */
REAL innerproductxy(REAL *x, 
					REAL *y, 
					INT *level, 
					INT k)
{
	INT i,n;
	REAL innerproduct;

	innerproduct = 0.0;
	n = level[k+1] - level[k];

	for	(i = 0; i < n; i++){
		innerproduct = innerproduct + x[level[k]+i]*y[level[k]+i];
	}

	return innerproduct;
}


/**
 * @brief set initial F
 *
 * @param x number x
 *
 * @return F
 */
REAL F(REAL x)
{
  REAL r = exp(x) * (3.0 + x) * x; 
  return r;
}
  
