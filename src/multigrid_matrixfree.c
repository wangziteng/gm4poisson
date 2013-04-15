/**
 * @file multigrid_matrixfree.c
 * @brief All the algorithm of v-cycle, pcg and full multigrid
 * @author Ziteng Wang
 * @version 1.0
 * @date 2013-04-02
 */
#include <stdio.h>
#include <stdlib.h>
#include "malloc.h"
#include <math.h>
#include <omp.h>

#include "multigrid_functs.h"

/**
 * @brief a v-cycle process of multigrid for 1d poisson problem
 *
 * @param u solution vector
 * @param b right hand vector
 * @param level a vector containing information of all grids. For exmaple, level[0] is the first position of level 0.
 * @param maxlevel number of levels
 */
void multigriditeration1d(REAL *u,
                          REAL *b,
                          INT *level,
                          INT maxlevel)
{
    INT n,i,j,k;
    REAL *r;

	// OpenMP settings
	//INT nthreads = omp_get_thread_num();
	int nthreads = 1;

	r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));

    // forward sweep 
    for (k = 0; k < maxlevel-1; k++) {
        // initial
		n = level[k+1] - level[k];
	#pragma omp parallel for private(i) num_threads(nthreads)
		for (i = 0; i < n; i++){
			r[i] = 0.0;
		}
		if (k>0){
		#pragma omp parallel for private(i) num_threads(nthreads)
			for (i = 0; i < (level[k+1]-level[k]);i++){
				u[level[k]+i] = 0.0;
			}
		}
        // pre-smoothing, G-S as smoother
        for (i = 0; i < 3; i++) {
//		#pragma omp parallel for private(i) num_threads(nthreads)
            for (j = 1; j < n; j=j+1) {
                u[level[k]+j] = (b[level[k]+j]+u[level[k]-1+j]+u[level[k]+1+j])/2;          
            }
			//#pragma omp parallel for private(i) num_threads(nthreads)
   //         for (j = 2; j < n; j=j+2) {
   //             u[level[k]+j] = (b[level[k]+j]+u[level[k]-1+j]+u[level[k]+1+j])/2;          
   //         }
   		}
		compute_r_1d(b, u, r, k, level);

        // restriction on coarser grids
		n = level[k+2]-level[k+1];
		b[level[k+1]] = 0.0;
//	#pragma omp parallel for private(i) num_threads(nthreads)
        for (j = 1; j < n; j++) {
			b[level[k+1]+j] = (2*r[level[k]+2*j]+r[level[k]+2*j-1]+r[level[k]+2*j+1]);  
        }   
    } 
    
    // coarsest grid
    if (k==maxlevel-1) {
        i = level[maxlevel-1];
        u[i+1] = b[i+1]/2;
    }

    // back sweep
    for (k = maxlevel-1; k > 0; k--) {
        n = level[k+1] - level[k];
        // interpolation to finer grids
//	#pragma omp parallel for private(i) num_threads(nthreads)
        for (i = 0; i < n-1; i++) {
			u[level[k-1]+2*i] += u[level[k]+i];
            u[level[k-1]+2*i+1] += (u[level[k]+i]+u[level[k]+i+1])/2;
        }
        u[level[k-1]+2*n-1] += u[level[k]+n-1]/2;

        // post-smoothing, G-S as smoother
        n = level[k] - level[k-1];
        for (i = 0; i < 3; i++) {
//		#pragma omp parallel for private(i) num_threads(nthreads)
            for (j = 1; j < n; j=j+1) {
                u[level[k-1]+j] = (b[level[k-1]+j]+u[level[k-1]-1+j]+u[level[k-1]+1+j])/2;    
            }
		//#pragma omp parallel for private(i) num_threads(nthreads)
  //          for (j = 2; j < n; j=j+2) {
  //              u[level[k-1]+j] = (b[level[k-1]+j]+u[level[k-1]-1+j]+u[level[k-1]+1+j])/2;    
  //          }

        }
    }
    free(r);
}

/**
 * @brief a v-cycle process of multigrid for 1d poisson problem
 *
 * @param u solution vector
 * @param b right hand vector
 * @param level a vector containing information of all grids. For exmaple, level[0] is the first position of level 0.
 * @param maxlevel number of all levels
 * @param nx number of cells on x direction
 * @param ny number of cells on y direction
 */
void multigriditeration2d(REAL *u,
                          REAL *b,
                          INT *level,
                          INT maxlevel,
                          INT nx,
                          INT ny)
{
    INT i,k,i1;
    REAL *r;
    INT *nxk, *nyk;
	int nthreads = 4;

	nxk = (INT *)malloc(maxlevel*sizeof(INT));
	nyk = (INT *)malloc(maxlevel*sizeof(INT));
	r = (REAL *)malloc(3*nx*ny*sizeof(REAL));
    nxk[0] = nx+1; nyk[0] = ny+1;

	// forward sweep
    for (k = 0; k < maxlevel-1; k++) {
        nxk[k+1] = ((int) (nxk[k]-1)/2)+1; nyk[k+1] = ((int) (nyk[k]-1)/2)+1;
		level[k+1] = level[k]+(nxk[k])*(nyk[k]);
        
        // pre-smoothing, GS as smoother
		// initial some vectors
	#pragma omp parallel for private(i) num_threads(nthreads)
		for (i = 0; i < (level[k+1]-level[k]); i++){
			r[level[k]+i] = 0.0;
		}
		if (k>0) {
		#pragma omp parallel for private(i) num_threads(nthreads)
			for (i = 0; i < (level[k+1]-level[k]); i++){
				u[level[k]+i] = 0.0;
			}
		}

		// Gauss-Seidel 2 colors
		for (i = 0; i < 3; i++){
			gsiteration_2color_2d(u, b, level, k, maxlevel, nxk, nyk);
		}           
		compute_r_2d(u, b, r, k, level, nxk, nyk);
        // restriction on coarser grids
        coarsergrid7pointrestriction2d(b, r, level, k, nxk, nyk);
	}

    // coarsest grid
    if(k==maxlevel-1){
        u[level[k]+4] = b[level[k]+4]/4;
		level[k+1] = level[k]+nxk[k]*nyk[k];
    }

    // back sweep
        
    // interpolation on finer grids
    for (k = maxlevel-1; k > 0; k--) {
		finergridinterpolation2d(u, level, k, nxk, nyk);
        // post-smoothing
		k = k-1;
		for (i1 = 0; i1<3; i1++){
			// Gauss-Seidel 2 colors
			gsiteration_2color_2d(u, b, level, k, maxlevel, nxk, nyk);	
		}
		compute_r_2d(u, b, r, k, level, nxk, nyk);
		k = k+1;
    }
    free(r);
}

/**
 * @brief a v-cycle process of multigrid for 1d poisson problem
 *
 * @param u solution vector
 * @param b right hand vector
 * @param level a vector containing information of all grids. For exmaple, level[0] is the first position of level 0.
 * @param maxlevel number of levels
 * @param nx number of cells on x direction
 * @param ny number of cells on y direction
 * @param nz number of cells on z direction
 */
void multigriditeration3d(REAL *u,
						  REAL *b,
                          INT *level,
                          INT maxlevel,
                          INT nx,
                          INT ny,
                          INT nz)
{
    INT n,i,k;
    REAL *r;
    INT *nxk,*nyk,*nzk,*presmoothtime,*prosmoothtime;
	// OpenMP settings
	//INT nthreads = omp_get_thread_num();
//	int nthreads = 2;

    nxk = (INT *)malloc((maxlevel+1)*sizeof(INT));
	nyk = (INT *)malloc((maxlevel+1)*sizeof(INT));
	nzk = (INT *)malloc((maxlevel+1)*sizeof(INT));
	presmoothtime = (INT *)malloc((maxlevel+1)*sizeof(INT));
	prosmoothtime = (INT *)malloc((maxlevel+1)*sizeof(INT));
	r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	nxk[0] = nx+1; nyk[0] = ny+1; nzk[0] = nz+1;
    n = level[1] - level[0];

	// set times of post and pre smoothing
	for (k = 0; k < maxlevel; k++){
		presmoothtime[k] = 3;
	}
	presmoothtime[0] = 1;
	presmoothtime[1] = 2;
	presmoothtime[2] = 3;
	for (k = 0; k < maxlevel; k++){
		prosmoothtime[k] = 3;
	}
	prosmoothtime[0] = 1;
	prosmoothtime[1] = 2;
	prosmoothtime[2] = 3;

	
    // forward sweep

	// finest level
	nxk[1] = (int) (nxk[0]+1)/2;
	nyk[1] = (int) (nyk[0]+1)/2;
	nzk[1] = (int) (nzk[0]+1)/2;

	for	(i = 0; i < presmoothtime[0]; i++){
		gsiteration_2color_3d(u, b, level, 0, maxlevel, nxk, nyk, nzk);
	}
	compute_r_3d(u, b, r, 0, level, nxk, nyk, nzk);
	// restriction on coarser grid
	coarsergrid7pointrestriction3d(b, r, level, 0, nxk, nyk, nzk);

	// coarser grids
    for (k = 1; k < maxlevel-1; k++) {

		nxk[k+1] = (int) (nxk[k]+1)/2;
		nyk[k+1] = (int) (nyk[k]+1)/2;
		nzk[k+1] = (int) (nzk[k]+1)/2;

		// initial vectors
		for (i = 0; i < level[k+1]-level[k]; i++){
			u[level[k]+i] = 0.0;
		}

        // pre-smoothing 
		for	(i = 0; i < presmoothtime[k]; i++){
			gsiteration_2color_3d(u, b, level, k, maxlevel, nxk, nyk, nzk);
//			gsiteration3d(u, b, level, k, maxlevel, nxk, nyk, nzk);
		}
		compute_r_3d(u, b, r, k, level, nxk, nyk, nzk);
		// restriction on coarser grid
		coarsergrid7pointrestriction3d(b, r, level, k, nxk, nyk, nzk);
	}
	free(r);

	// coarsest level
	u[level[maxlevel-1]+13] = b[level[maxlevel-1]+13]/6;

	// back sweep
	for (k = maxlevel-1; k > 0; k--) {
		// interpolation from coarser grid		
		finergridinterpolation3d(u, level, k, nxk, nyk, nzk);
		// post smoothing
		k = k-1;
		for (i = 0; i < prosmoothtime[k]; i++){
			gsiteration_2color_3d(u, b, level, k, maxlevel, nxk, nyk, nzk);
		}
		k = k+1;
	}
	free(nxk);
	free(nyk);
	free(nzk);
}

/**
 * @brief Full multigrid for 1d poisson equation
 *
 * @param u solution vector
 * @param b right hand vector
 * @param level a vector containing information of all grids. For exmaple, level[0] is the first position of level 0
 * @param maxlevel number of all levels
 * @param nx number of cells in x direction
 * @param rtol relative tolerance
 */
void fullmultigrid_1d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
                      INT nx,
					  REAL rtol)
{
    INT n,i,j,k;
    REAL *r;
    INT *nxk;
	REAL norm_r;

	// initial
    nxk = (INT *)malloc((maxlevel+1)*sizeof(INT));
	r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	nxk[0] = nx; 
    n = level[1] - level[0];

    for (k = 0; k < maxlevel-1; k++) {
		nxk[k+1] = (int)(nxk[k]+1)/2;

		// initial some vectors
		for (i = 0; i < (level[k+1]-level[k]); i++){
			r[level[k]+i] = 0.0;
		}
		if (k>0) {
			for (i = 0; i < (level[k+1]-level[k]); i++){
				u[level[k]+i] = 0.0;
			}
		}

		// GS iteration
        for (i = 0; i < 3; i++) {
            for (j = 1; j < n; j++) {
                u[level[k]+j] = (b[level[k]+j]+u[level[k]-1+j]+u[level[k]+1+j])/2;          
            }
   		}
		compute_r_1d(b, u, r, k, level);
		norm_r = computenorm(r, level, k);
		printf("pre-residue of level[%d] = %f\n",k,norm_r);
        // restriction on coarser grids
		n = level[k+2]-level[k+1];
        for (j = 0; j < n; j++) {
            if (j==0) {
                b[level[k+1]] = 0;
            }
		    else{
				b[level[k+1]+j] = (2*r[level[k]+2*j]+r[level[k]+2*j-1]+r[level[k]+2*j+1])/4*4;  
		    }
        }   
    }

    // coarsest grid
    if(k==maxlevel-1){
        i = level[maxlevel-1];
        u[i+1] = b[i+1]/2;
    }

    // FULL multigrid
    while (k>0) {
        n = level[k+1] - level[k];
        // interpolation to finer grids
        for (i = 0; i < n; i++) {
			u[level[k-1]+2*i] += u[level[k]+i];
            if(i<n-1){
                u[level[k-1]+2*i+1] += (u[level[k]+i]+u[level[k]+i+1])/2;
            }
            else {
                u[level[k-1]+2*i+1] += u[level[k]+i]/2;
            }
        }
		compute_r_1d(u, b, r, k-1, level);
        norm_r = computenorm(r, level, k-1);
		printf("full-residue of level[%d] = %f\n",k-1,norm_r);
		k = k-1;
		// GS iteration
        for (i = 0; i < 3; i++) {
            for (j = 1; j < n; j++) {
                u[level[k]+j] = (b[level[k]+j]+u[level[k]-1+j]+u[level[k]+1+j])/2;          
            }
   		}
		compute_r_1d(u, b, r, k, level);
        norm_r = computenorm(r, level, k);
		printf("residue of level[%d] = %f\n",k,norm_r);
        while (norm_r>rtol) {
            vcycle_at_levelk_1d(u, b, r, level, k, maxlevel);
            compute_r_1d(u, b, r, k, level);
            norm_r = computenorm(r, level, k);
        }
    }
    free(r);
}


static vcycle_at_levelk_1d(REAL *u, 
	 				       REAL *b, 
						   REAL *r, 
						   INT *level, 
						   INT k,
						   INT maxlevel)
{
	INT i,j,h,n;
	REAL norm_r0,norm_r;

	// forward sweep
    for (h = k; h < maxlevel-1; h++) {
        
        // pre-smoothing, GS as smoother
		// initial some vectors
		for (i = 0; i < (level[h+1]-level[h]); i++){
			r[level[h]+i] = 0.0;
		}
		if (h>k) {
			for (i = 0; i < (level[h+1]-level[h]); i++){
				u[level[h]+i] = 0.0;
			}
		}

		// GS iteration
        for (i = 0; i < 3; i++) {
            for (j = 1; j < n; j++) {
                u[level[h]+j] = (b[level[h]+j]+u[level[h]-1+j]+u[level[h]+1+j])/2;          
            }
   		}
		compute_r_1d(b, u, r, h, level);
		norm_r = computenorm(r, level, h);
		printf("pre-residue of level[%d] = %f\n",h,norm_r);
        // restriction on coarser grids
		n = level[h+2]-level[h+1];
        for (j = 0; j < n; j++) {
            if (j==0) {
                b[level[h+1]] = 0;
            }
		    else{
				b[level[h+1]+j] = (2*r[level[h]+2*j]+r[level[h]+2*j-1]+r[level[h]+2*j+1])/4*4;  
		    }
        }   
	}

    // coarsest grid
    if(h==maxlevel-1){
        i = level[maxlevel-1];
        u[i+1] = b[i+1]/2;
    }

    // back sweep
        
    // interpolation on finer grids
    for (h = maxlevel-1; h > k; h--) {
        n = level[h+1] - level[h];
        // interpolation to finer grids
        for (i = 0; i < n; i++) {
			u[level[h-1]+2*i] += u[level[h]+i];
            if(i<n-1){
                u[level[h-1]+2*i+1] += (u[level[h]+i]+u[level[h]+i+1])/2;
            }
            else {
                u[level[h-1]+2*i+1] += u[level[h]+i]/2;
            }
        }
		compute_r_1d(u, b, r, h-1, level);
        norm_r = computenorm(r, level, h-1);
		printf("full-residue of level[%d] = %f\n",h-1,norm_r);
		// post-smoothing
		h = h-1;
		// GS iteration
        for (i = 0; i < 3; i++) {
            for (j = 1; j < n; j++) {
                u[level[h]+j] = (b[level[h]+j]+u[level[h]-1+j]+u[level[h]+1+j])/2;          
            }
   		}
		compute_r_1d(u, b, r, h, level);
		norm_r0 = computenorm(r, level, h);
		printf("norm_r0 = %f\n",norm_r0);
		h = h+1;
    }
}

/**
 * @brief Full multigrid for 2d poisson problem
 *
 * @param u solution vector
 * @param b right hand vector
 * @param level a vector containing information of all grids. For exmaple, level[0] is the first position of level 0
 * @param maxlevel number of total levels
 * @param nx number of cells in x direction
 * @param ny number of cells in y direction
 * @param rtol relative tolerance
 */
void fullmultigrid_2d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
                      INT nx,
                      INT ny,
					  REAL rtol)
{
    INT n,i,k;
    REAL *r;
    INT *nxk,*nyk;
	REAL norm_r;

	// initial
    nxk = (INT *)malloc(maxlevel*sizeof(INT));
	nyk = (INT *)malloc(maxlevel*sizeof(INT));
	r = (REAL *)malloc(3*nx*ny*sizeof(REAL));
	nxk[0] = nx+1; nyk[0] = ny+1;
    n = level[1] - level[0];

    for (k = 0; k < maxlevel-1; k++) {
		nxk[k+1] = (int)(nxk[k]+1)/2;
		nyk[k+1] = (int)(nyk[k]+1)/2;
		level[k+1] = level[k]+nxk[k]*nyk[k];
        // GS iteration

		// initial some vectors
		for (i = 0; i < (level[k+1]-level[k]); i++){
			r[level[k]+i] = 0.0;
		}
		if (k>0) {
			for (i = 0; i < (level[k+1]-level[k]); i++){
				u[level[k]+i] = 0.0;
			}
		}
		// Gauss-Seidel 2 colors
		for (i = 0; i < 1; i++){
			gsiteration_2color_2d(u, b, level, k, maxlevel, nxk, nyk);
		}           
		compute_r_2d(u, b, r, k, level, nxk, nyk);
		norm_r = computenorm(r, level, k);
		printf("pre-residue of level[%d] = %f\n",k,norm_r);
        // restriction on coarser grids
        coarsergrid7pointrestriction2d(b, r, level, k, nxk, nyk);
    }

    // coarsest grid
    if(k==maxlevel-1){
        u[level[k]+4] = b[level[k]+4]/4;
		level[k+1] = level[k]+nxk[k]*nyk[k];
    }

    // FULL multigrid
    while (k>0) {
        // interpolation from coarser grid
        finergridinterpolation2d(u, level, k, nxk, nyk);
		compute_r_2d(u, b, r, k-1, level, nxk, nyk);
		for (i=0;i<(level[k]-level[k-1]);i++){
			printf("u[%d]=%lf\n",level[k-1]+i,u[level[k-1]+i]);
		}
        norm_r = computenorm(r, level, k-1);
		printf("full-residue of level[%d] = %f\n",k-1,norm_r);
		k = k-1;
        for (i = 0; i < 3; i++) {
            gsiteration_2color_2d(u, b, level, k, maxlevel, nxk, nyk);
        }
        compute_r_2d(u, b, r, k, level, nxk, nyk);
        norm_r = computenorm(r, level, k);
		printf("residue of level[%d] = %f\n",k,norm_r);
        while (norm_r>rtol) {
            vcycle_at_levelk_2d(u, b, r, level, k, maxlevel, nxk, nyk);
            compute_r_2d(u, b, r, k, level, nxk, nyk);
            norm_r = computenorm(r, level, k);
        }
    }
    free(r);
}


static vcycle_at_levelk_2d(REAL *u, 
		  			       REAL *b, 
						   REAL *r, 
						   INT *level, 
						   INT k, 
						   INT maxlevel, 
						   INT *nxk, 
						   INT *nyk)
{
	INT i,h,i1;
	REAL norm_r0;

	// forward sweep
    for (h = k; h < maxlevel-1; h++) {
        
        // pre-smoothing, GS as smoother
		// initial some vectors
		for (i = 0; i < (level[h+1]-level[h]); i++){
			r[level[h]+i] = 0.0;
		}
		if (h>k) {
			for (i = 0; i < (level[h+1]-level[h]); i++){
				u[level[h]+i] = 0.0;
			}
		}

		// Gauss-Seidel 2 colors
		for (i = 0; i < 1; i++){
			gsiteration_2color_2d(u, b, level, h, maxlevel, nxk, nyk);
		}           
		compute_r_2d(u, b, r, h, level, nxk, nyk);
	    norm_r0 = computenorm(r, level, h);
		printf("norm_r0 = %f\n",norm_r0);
        // restriction on coarser grids
        coarsergrid7pointrestriction2d(b, r, level, h, nxk, nyk);
	}

    // coarsest grid
    if(h==maxlevel-1){
        u[level[h]+4] = b[level[h]+4]/4;
    }

    // back sweep
        
    // interpolation on finer grids
    for (h = maxlevel-1; h > k; h--) {
		finergridinterpolation2d(u, level, h, nxk, nyk);
        // post-smoothing
		h = h-1;
		for (i1 = 0; i1<1; i1++){
			// Gauss-Seidel 2 colors
			gsiteration_2color_2d(u, b, level, h, maxlevel, nxk, nyk);	
		}
		compute_r_2d(u, b, r, h, level, nxk, nyk);
		norm_r0 = computenorm(r, level, h);
		printf("norm_r0 = %f\n",norm_r0);
		h = h+1;
    }
}

/**
 * @brief Full multigrid for 3d poisson problem
 *
 * @param u solution vector
 * @param b right hand vector
 * @param level a vector containing information of all grids. For exmaple, level[0] is the first position of level 0 
 * @param maxlevel number of total levels
 * @param nx number of cells in x direction
 * @param ny number of cells in y direction
 * @param nz number of cells in z direction
 * @param rtol relative tolerance
 */
void fullmultigrid_3d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
                      INT nx,
                      INT ny,
                      INT nz,
					  REAL rtol)
{
    INT n,i,k;
    REAL *r;
    INT *nxk,*nyk,*nzk;
	REAL res, norm_r;

	// initial
    nxk = (INT *)malloc(maxlevel*sizeof(INT));
	nyk = (INT *)malloc(maxlevel*sizeof(INT));
	nzk = (INT *)malloc(maxlevel*sizeof(INT));
	r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	nxk[0] = nx; nyk[0] = ny; nzk[0] = nz;
    n = level[1] - level[0];

    for (k = 0; k < maxlevel-1; k++) {
		nxk[k+1] = (int)(nxk[k]+1)/2;
		nyk[k+1] = (int)(nyk[k]+1)/2;
		nzk[k+1] = (int)(nzk[k]+1)/2;
        // GS iteration
        for (i = 0; i < 3; i++) {
            gsiteration_2color_3d(u, b, level, k, maxlevel, nxk, nyk, nzk);
        }
		compute_r_3d(u, b, r, k, level, nxk, nyk, nzk);
		res = computenorm(r,level,k);
		printf("residue of level[%d] = %f\n",k,res);
        // restriction on coarser grid
        coarsergrid7pointrestriction3d(b, r, level, k, nxk, nyk, nzk);
    }

    // coarsest grid
    u[level[maxlevel]-1] = b[level[maxlevel]-1]/6;

    // FULL multigrid
    while (k>0) {
        // interpolation from coarser grid
        finergridinterpolation3d(u, level, k, nxk, nyk, nzk);
		compute_r_3d(u, b, r, k-1, level, nxk, nyk, nzk);
        norm_r = computenorm(r, level, k-1);
		printf("residue of level[%d] = %f\n",k-1,norm_r);
		k = k-1;
        for (i = 0; i < 3; i++) {
            gsiteration_2color_3d(u, b, level, k, maxlevel, nxk, nyk, nzk);
        }
        compute_r_3d(u, b, r, k, level, nxk, nyk, nzk);
        norm_r = computenorm(r, level, k);
		printf("residue of level[%d] = %f\n",k,norm_r);
        while (norm_r>rtol) {
            vcycle_at_levelk3d(u, b, r, level, k, maxlevel, nxk, nyk, nzk);
            compute_r_3d(u, b, r, k, level, nxk, nyk, nzk);
            norm_r = computenorm(r, level, k);
        }
    }
    free(r);
}


static vcycle_at_levelk3d(REAL *u, 
					      REAL *b, 
					      REAL *r, 
					      INT *level, 
					      INT k, 
					      INT maxlevel, 
					      INT *nxk, 
					      INT *nyk, 
					      INT *nzk)
{
	INT i,j,h;
	REAL res;

    for (h = k; h < maxlevel-1; h++) {

        // pre-smoothing 
		for	(i = 0; i < 3; i++){
			gsiteration_2color_3d(u, b, level, h, maxlevel, nxk, nyk, nzk);

			compute_r_3d(u, b, r, h, level, nxk, nyk, nzk);
			res = computenorm(r,level,h);
			printf("residue of level[%d] = %f\n",h,res);
		}

		// restriction on coarser grid
		coarsergrid7pointrestriction3d(b, r, level, h, nxk, nyk, nzk);
	}

	// coarsest level
	u[level[maxlevel-1]] = b[level[maxlevel-1]]/6;

	// back sweep
	for (h = maxlevel-1; h > k; h--) {
		// interpolation from coarser grid		
		finergridinterpolation3d(u, level, h, nxk, nyk, nzk);

		// post smoothing
		h = h-1;
		for (i = 0; i < 3; i++){
			gsiteration_2color_3d(u, b, level, h, maxlevel, nxk, nyk, nzk);

			compute_r_3d(u, b, r, h, level, nxk, nyk, nzk);
			res = computenorm(r,level,h);
			printf("residue of level[%d] = %f\n",h,res);
		}
		h = h+1;
	}
}
void pcg_1d(REAL *u,
            REAL *b,
            INT *level,
			INT maxlevel,
            INT nx,
            REAL rtol,
			INT maxiteration)
{
    INT i,k,done;
	INT nxk[1];
    REAL *p, *r, *z, *q;
    REAL rh0, rh1, rh2, alfa, beta, normb, normr, resid;

    p = (REAL *)malloc(level[1]*sizeof(REAL));
    r = (REAL *)malloc(3*nx*sizeof(REAL));
	z = (REAL *)malloc(3*nx*sizeof(REAL));
	q = (REAL *)malloc(level[1]*sizeof(REAL));
	done = 0;
    k = 0;
	nxk[0] = nx+1; 

    // initial residue and other vector
	for (i = 0; i < 3*nx; i++){
		z[i] = 0.0;
		r[i] = 0.0;
	}
	for (i = 0; i < level[1]; i++){
		p[i] = 0.0;
		q[i] = 0.0;
	}

	compute_r_1d(u, b, r, 0, level);
	normr = computenorm(r, level, 0);
	normb = computenorm(b, level, 0);
	if (normb==0.0) {
		normb==1.0;
	}
	if ((resid = normr / normb) <= rtol) {
		return;
	}

    multigriditeration1d(z, r, level, maxlevel);
	//xequaly(z, r, level, 0);	
    rh0 = innerproductxy(r, z, level, 0);
    xequaly(p, z, level, 0);

    while (!done && k < maxiteration) {
		// init z
        for (i = 0; i < level[1]; i++) z[i] = 0;

        // calculating alpha
        xequalay_1d(q, p, level, 0);
		rh2 = innerproductxy(q, p, level, 0);
		alfa = rh0/rh2;

        // update vector u, r
		xequalypcz(u, u, alfa, p, level, 0);
		xequalypcz(r, r, (-alfa), q, level, 0);
		normr = computenorm(r, level, 0);
		if ((resid = normr / normb) <= rtol) {
			done = 1;     
		}

        // update z and beta
        multigriditeration1d(z, r, level, maxlevel);
		//xequaly(z, r, level, 0);	
        rh1 = innerproductxy(r, z, level, 0);
        beta = rh1 / rh0;

        // update p
        xequalypcz(p, z, beta, p, level, 0);

		rh0 = rh1;
		k++;
	}
	printf("iteration number = %d\n",k);
    free(r);
    free(q);
    free(p);
	free(z);
}

void pcg_2d(REAL *u,
            REAL *b,
            INT *level,
			INT maxlevel,
            INT nx,
            INT ny,
            REAL rtol,
			INT maxiteration)
{
    INT i,k,done;
	INT nxk[1], nyk[1];
    REAL *p, *r, *z, *q;
    REAL rh0, rh1, rh2, alfa, beta, normb, normr, resid;

    p = (REAL *)malloc(level[1]*sizeof(REAL));
    r = (REAL *)malloc(3*nx*ny*sizeof(REAL));
	z = (REAL *)malloc(3*nx*ny*sizeof(REAL));
	q = (REAL *)malloc(level[1]*sizeof(REAL));
	done = 0;
    k = 0;
	nxk[0] = nx+1; nyk[0] = ny+1; 

    // initial residue and other vector
	for (i = 0; i < 3*nx*ny; i++){
		z[i] = 0.0;
		r[i] = 0.0;
	}
	for (i = 0; i < level[1]; i++){
		p[i] = 0.0;
		q[i] = 0.0;
	}

	compute_r_2d(u, b, r, 0, level, nxk, nyk);
	normr = computenorm(r, level, 0);
	normb = computenorm(b, level, 0);
	if (normb==0.0) {
		normb==1.0;
	}
	if ((resid = normr / normb) <= rtol) {
		return;
	}

    multigriditeration2d(z, r, level, maxlevel, nx, ny);
	//xequaly(z, r, level, 0);	
    rh0 = innerproductxy(r, z, level, 0);
    xequaly(p, z, level, 0);

    while (!done && k < maxiteration) {
		// init z
        for (i = 0; i < level[1]; i++) z[i] = 0;

        // calculating alpha
        xequalay_2d(q, p, level, 0, nxk, nyk);
		rh2 = innerproductxy(q, p, level, 0);
		alfa = rh0/rh2;

        // update vector u, r
		xequalypcz(u, u, alfa, p, level, 0);
		xequalypcz(r, r, (-alfa), q, level, 0);
		normr = computenorm(r, level, 0);
		if ((resid = normr / normb) <= rtol) {
			done = 1;     
		}

        // update z and beta
        multigriditeration2d(z, r, level, maxlevel, nx, ny);
		//xequaly(z, r, level, 0);	
        rh1 = innerproductxy(r, z, level, 0);
        beta = rh1 / rh0;

        // update p
        xequalypcz(p, z, beta, p, level, 0);

		rh0 = rh1;
		k++;
	}
	printf("iteration number = %d\n",k);
    free(r);
    free(q);
    free(p);
	free(z);
}

/**
 * @brief A multigrid-preconditioned CG process for 3d Poisson problem
 *
 * @param u solution vector
 * @param b right hand vector
 * @param level a vector containing information of all grids. For exmaple, level[0] is the first position of level 0
 * @param maxlevel number of total levels 
 * @param nx number of cells in x direction
 * @param ny number of cells in y direction
 * @param nz number of cells in z direction
 * @param rtol relative tolerance
 * @param maxiteration maximum CG iteration number
 */
void pcg_3d(REAL *u,
            REAL *b,
            INT *level,
			INT maxlevel,
            INT nx,
            INT ny,
            INT nz, 
            REAL rtol,
			INT maxiteration)
{
    INT i,k,done;
	INT nxk[1], nyk[1], nzk[1];
    REAL *p, *r, *z, *q;
    REAL rh0, rh1, rh2, alfa, beta, normb, normr, resid;

    p = (REAL *)malloc(level[1]*sizeof(REAL));
    r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	z = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	q = (REAL *)malloc(level[1]*sizeof(REAL));
	done = 0;
    k = 0;
	nxk[0] = nx+1; nyk[0] = ny+1; nzk[0] = nz+1;

    // initial residue and other vector
	for (i = 0; i < level[1]; i++){
		z[i] = 0.0;
		r[i] = 0.0;
	}
	for (i = 0; i < level[1]; i++){
		p[i] = 0.0;
		q[i] = 0.0;
	}
	compute_r_3d(u, b, r, 0, level, nxk, nyk, nzk);
	normr = computenorm(r, level, 0);
	normb = computenorm(b, level, 0);
	if (normb==0.0) {
		normb==1.0;
	}
	if ((resid = normr / normb) <= rtol) {
		return;
	}

	multigriditeration3d(z, r, level, maxlevel, nx, ny, nz);
	//xequaly(z, r, level, 0);	
    rh0 = innerproductxy(r, z, level, 0);
    xequaly(p, z, level, 0);

    while (!done && k < maxiteration) {
		// init z
        for (i = 0; i < level[1]; i++) z[i] = 0;

        // calculating alpha
        xequalay_3d(q, p, level, 0, nxk, nyk, nzk);
		rh2 = innerproductxy(q, p, level, 0);
		alfa = rh0/rh2;

        // update vector u, r
		xequalypcz(u, u, alfa, p, level, 0);
		xequalypcz(r, r, (-alfa), q, level, 0);
		normr = computenorm(r, level, 0);
		if ((resid = normr / normb) <= rtol) {
			done = 1;     
		}

        // update z and beta
        multigriditeration3d(z, r, level, maxlevel, nx, ny, nz);
		//xequaly(z, r, level, 0);	
        rh1 = innerproductxy(r, z, level, 0);
        beta = rh1 / rh0;

        // update p
        xequalypcz(p, z, beta, p, level, 0);

		rh0 = rh1;
		k++;
	}
	printf("iteration number = %d\n",k);
    free(r);
    free(q);
    free(p);
	free(z);
}

/**
 * @brief 2 Color G-S iteration for 2d poisson problem
 *
 * @param u solution vector
 * @param b right hand vector
 * @param level a vector containing information of all grids. For exmaple, level[0] is the first position of level 0
 * @param k current level
 * @param maxlevel number of total levels
 * @param nxk number of cells in x direction of level k
 * @param nyk number of cells in y direction of level k
 */
void gsiteration_2color_2d(REAL *u,
						   REAL *b,
						   INT *level,
                           INT k,
                           INT maxlevel,
                           INT *nxk, 
                           INT *nyk)
{
	INT h,i;
	int nthreads = 1;
	// red
    for (h = 1; h < nyk[k]-1; h = h+2) {
		#pragma omp parallel for private(i) num_threads(nthreads)
        for (i = 1; i < nxk[k]-1; i=i+2) {
			u[level[k]+nxk[k]*h+i] = (b[level[k]+nxk[k]*h+i]+u[level[k]+nxk[k]*h+i+1]+u[level[k]+nxk[k]*h+i-1]+
									  u[level[k]+nxk[k]*(h+1)+i]+u[level[k]+nxk[k]*(h-1)+i])/4;
        }
	}
    for (h = 2; h < nyk[k]-1; h = h+2) {
		#pragma omp parallel for private(i) num_threads(nthreads)
        for (i = 2; i < nxk[k]-1; i=i+2) {
            u[level[k]+nxk[k]*h+i] = (b[level[k]+nxk[k]*h+i]+u[level[k]+nxk[k]*h+i+1]+u[level[k]+nxk[k]*h+i-1]+	
									  u[level[k]+nxk[k]*(h+1)+i]+u[level[k]+nxk[k]*(h-1)+i])/4;
        }
	}
	// black
    for (h = 1; h < nyk[k]-1; h = h+2) {
		#pragma omp parallel for private(i) num_threads(nthreads)
        for (i = 2; i < nxk[k]-1; i=i+2) {
			u[level[k]+nxk[k]*h+i] = (b[level[k]+nxk[k]*h+i]+u[level[k]+nxk[k]*h+i+1]+u[level[k]+nxk[k]*h+i-1]+	
									  u[level[k]+nxk[k]*(h+1)+i]+u[level[k]+nxk[k]*(h-1)+i])/4;
        }
	}
    for (h = 2; h < nyk[k]-1; h = h+2) {
		#pragma omp parallel for private(i) num_threads(nthreads)
        for (i = 1; i < nxk[k]-1; i=i+2) {
			u[level[k]+nxk[k]*h+i] = (b[level[k]+nxk[k]*h+i]+u[level[k]+nxk[k]*h+i+1]+u[level[k]+nxk[k]*h+i-1]+
									  u[level[k]+nxk[k]*(h+1)+i]+u[level[k]+nxk[k]*(h-1)+i])/4;
        }
	}
}


/**
 * @brief 2 Color G-S iteration for 2d poisson problem
 *
 * @param u solution vector
 * @param b right hand vector
 * @param level a vector containing information of all grids. For exmaple, level[0] is the first position of level 0
 * @param k current level
 * @param maxlevel number of total levels
 * @param nxk number of cells in x direction of level k
 * @param nyk number of cells in y direction of level k
 * @param nzk number of cells in z direction of level k
 */
void gsiteration_2color_3d(REAL *u,
                           REAL *b,
                           INT *level,
                           INT k,
                           INT maxlevel,
                           INT *nxk, 
                           INT *nyk, 
                           INT *nzk)
{
	INT i,j,h;
	INT i0,j0,j1,j2,k0,k1,k2,k3,k4,k5,k6;
	const int levelk = level[k];
	const int nxkk = nxk[k];
	const int nykk = nyk[k];
	const int nzkk = nzk[k];
	const int nxyk = nxkk*nykk;

	//i0 = level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h;
	//ir = level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1;
	//il = level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1;
	//ifront = level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1;
	//ia = level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h;
	//iu = level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h;
	//id = level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h;
	

	// red point of 2*i,2*j,2*h
    for (i = 2; i < nzkk-1; i = i+2) {
		i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
			j0 = i0+j*nxkk;
			j1 = i0+(j+1)*nxkk;
			j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
				k0 = j0+h;
				//k1 = j0+h-1;
				//k2 = j0+h+1;
				k3 = j1+h;
				k4 = j2+h;
				k5 = k0+nxyk;
				k6 = k0-nxyk;
                //u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h] =
                //    (b[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
                //     u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1]+
                //     u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1]+
                //     u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h]+
                //     u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h]+
                //     u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h]+
                //     u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]+h])/6;
				u[k0]=(b[k0]+u[k0+1]+u[k0-1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
			}
        }
    }
    // red point of 2*i,2*j+1,2*h+1
	for (i = 2; i < nzkk-1; i = i+2) {
		i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
			j0 = i0+j*nxkk;
			j1 = i0+(j+1)*nxkk;
			j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
				k0 = j0+h;
				//k1 = j0+h-1;
				//k2 = j0+h+1;
				k3 = j1+h;
				k4 = j2+h;
				k5 = k0+nxyk;
				k6 = k0-nxyk;
				u[k0]=(b[k0]+u[k0+1]+u[k0-1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
			}
        }
    }
    //for (i = 2; i < nzk[k]-1; i = i+2) {
    //    for (j = 1; j < nyk[k]-1; j = j+2) {
    //        for (h = 1; h < nxk[k]-1; h = h+2) {
    //            u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h] =
    //            (b[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //             u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1]+
    //             u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1]+
    //             u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h]+
    //             u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h]+
    //             u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //             u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]+h])/6;
    //        }
    //    }
    //}
    // points of (2*i+1,2*j+1,2h)
	for (i = 1; i < nzkk-1; i = i+2) {
		i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
			j0 = i0+j*nxkk;
			j1 = i0+(j+1)*nxkk;
			j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
				k0 = j0+h;
				//k1 = j0+h-1;
				//k2 = j0+h+1;
				k3 = j1+h;
				k4 = j2+h;
				k5 = k0+nxyk;
				k6 = k0-nxyk;
				u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
			}
        }
    }
    //for (i = 1; i < nzk[k]-1; i = i+2) {
    //    for (j = 1; j < nyk[k]-1; j = j+2) {
    //        for (h = 2; h < nxk[k]-1; h = h+2) {
    //            u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h] =
    //                (b[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //                 u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1]+
    //                 u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1]+
    //                 u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h]+
    //                 u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h]+
    //                 u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //                 u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]+h])/6;
    //        }
    //    }
    //}
    // points of 2*i+1, 2*j, 2*h+1
    //for (i = 1; i < nzk[k]-1; i = i+2) {
    //    for (j = 2; j < nyk[k]-1; j = j+2) {
    //        for (h = 1; h < nxk[k]-1; h = h+2) {
    //            u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h] =
    //                (b[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //                 u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1]+
    //                 u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1]+
    //                 u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h]+
    //                 u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h]+
    //                 u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //                 u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]+h])/6;
    //        }
    //    }
    //}
	for (i = 1; i < nzkk-1; i = i+2) {
		i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
			j0 = i0+j*nxkk;
			j1 = i0+(j+1)*nxkk;
			j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
				k0 = j0+h;
				//k1 = j0+h-1;
				//k2 = j0+h+1;
				k3 = j1+h;
				k4 = j2+h;
				k5 = k0+nxyk;
				k6 = k0-nxyk;
				u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
			}
        }
    }
    
	// Black points 
    // 2*i,2*j,2*h+1
    //for (i = 2; i < nzk[k]-1; i = i+2) {
    //    for (j = 2; j < nyk[k]-1; j = j+2) {
    //        for (h = 1; h < nxk[k]-1; h = h+2) {
    //            u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h] =
    //            (b[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //             u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1]+
    //             u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1]+
    //             u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h]+
    //             u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h]+
    //             u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //             u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]+h])/6;
    //        }
    //    }
    //}
	for (i = 2; i < nzkk-1; i = i+2) {
		i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
			j0 = i0+j*nxkk;
			j1 = i0+(j+1)*nxkk;
			j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
				k0 = j0+h;
				//k1 = j0+h-1;
				//k2 = j0+h+1;
				k3 = j1+h;
				k4 = j2+h;
				k5 = k0+nxyk;
				k6 = k0-nxyk;
				u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
			}
        }
    }
    // 2*i,2*j+1,2*h
/*    for (i = 2; i < nzk[k]-1; i = i+2) {
        for (j = 1; j < nyk[k]-1; j = j+2) {
            for (h = 2; h < nxk[k]-1; h = h+2) {
                u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h] =
                (b[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
                 u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1]+
                 u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1]+
                 u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h]+
                 u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h]+
                 u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h]+
                 u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]+h])/6;
            }
        }
    }
 */ 
	for (i = 2; i < nzkk-1; i = i+2) {
		i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
			j0 = i0+j*nxkk;
			j1 = i0+(j+1)*nxkk;
			j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
				k0 = j0+h;
				//k1 = j0+h-1;
				//k2 = j0+h+1;
				k3 = j1+h;
				k4 = j2+h;
				k5 = k0+nxyk;
				k6 = k0-nxyk;
				u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
			}
        }
    }	
	// 2*i+1,2*j,2*h
    //for (i = 1; i < nzk[k]-1; i = i+2) {
    //    for (j = 2; j < nyk[k]-1; j = j+2) {
    //        for (h = 2; h < nxk[k]-1; h = h+2) {
    //            u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h] =
    //            (b[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //             u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1]+
    //             u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1]+
    //             u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h]+
    //             u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h]+
    //             u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //             u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]+h])/6;
    //        }
    //    }
    //}
	for (i = 1; i < nzkk-1; i = i+2) {
		i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
			j0 = i0+j*nxkk;
			j1 = i0+(j+1)*nxkk;
			j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
				k0 = j0+h;
				//k1 = j0+h-1;
				//k2 = j0+h+1;
				k3 = j1+h;
				k4 = j2+h;
				k5 = k0+nxyk;
				k6 = k0-nxyk;
				u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
			}
        }
    }
	for (i = 1; i < nzkk-1; i = i+2) {
		i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
			j0 = i0+j*nxkk;
			j1 = i0+(j+1)*nxkk;
			j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
				k0 = j0+h;
				//k1 = j0+h-1;
				//k2 = j0+h+1;
				k3 = j1+h;
				k4 = j2+h;
				k5 = k0+nxyk;
				k6 = k0-nxyk;
				u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
			}
        }
    }
    // 2*i+1,2*j+1,2*h+1
//    for (i = 1; i < nzk[k]-1; i = i+2) {
//        for (j = 1; j < nyk[k]-1; j = j+2) {
//            for (h = 1; h < nxk[k]-1; h = h+2) {
//                u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h] =
//                (b[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
//                 u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1]+
//                 u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h-1]+
//                 u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h]+
//                 u[level[k]+i*nxk[k]*nyk[k]+(j-1)*nxk[k]+h]+
//                 u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h]+
//                 u[level[k]+(i-1)*nxk[k]*nyk[k]+j*nxk[k]+h])/6;
//            }
//        }
//    }
}

static coarsergrid7pointrestriction2d(REAL *b, 
									REAL *r, 
									INT *level, 
									INT k, 
									INT *nxk, 
									INT *nyk)
{
	INT i,j;
	int nthreads = 1;

	for (i = 1; i < nyk[k+1]-1; i++){
//		#pragma omp parallel for private(j) num_threads(nthreads)
		for (j = 1; j < nxk[k+1]-1; j++){
			b[level[k+1]+i*nxk[k+1]+j] = (r[level[k]+2*i*nxk[k]+2*j]*2+
										  r[level[k]+2*i*nxk[k]+2*j+1]+
										  r[level[k]+2*i*nxk[k]+2*j-1]+
										  r[level[k]+(2*i+1)*nxk[k]+2*j]+
										  r[level[k]+(2*i-1)*nxk[k]+2*j]+
										  r[level[k]+(2*i+1)*nxk[k]+2*j+1]+
										  r[level[k]+(2*i-1)*nxk[k]+2*j-1])/2;
		}
		b[level[k+1]+(i+1)*nxk[k+1]-1] = (r[level[k]+(2*i+1)*nxk[k]-2]+
										  r[level[k]+(2*i)*nxk[k]-2])/2;
	}
//	#pragma omp parallel for private(j) num_threads(nthreads)
	for (j = 1; j < nxk[k+1]-1; j++){
		b[level[k+1]+(nyk[k+1]-1)*nxk[k+1]+j] = (r[level[k]+(nyk[k]-1)*nxk[k]+2*j-1]+
												 r[level[k]+(nyk[k]-2)*nxk[k]+2*j-1])/2;
	}
	b[level[k+1]+nxk[k+1]*nyk[k+1]-1] = r[level[k]+(nyk[k]-1)*nxk[k]-2]/2;
}

static coarsergrid7pointrestriction3d(REAL *b, 
									REAL *r, 
									INT *level, 
									INT k, 
									INT *nxk, 
									INT *nyk,
									INT *nzk)
{
	INT i,j,h;
	int i0,j0,j1,j2,j3,k0,k1,k2,k3,k4,k5,k6;
	int i01,j01,j11,j21,j31,k01,k11,k21,k31;
    const int levelk = level[k];
	const int levelk1 = level[k+1];
	const int nxkk = nxk[k];
	const int nxkk1 = nxk[k+1];
	const int nykk = nyk[k];
	const int nykk1 = nyk[k+1];
	const int nzkk = nzk[k];
	const int nzkk1 = nzk[k+1];

	int nxyk = nxkk*nykk;
	int nxyk1 = nxkk1*nykk1;

	for (i = 1; i < nzkk1-1; i++){
		i0 = levelk+2*i*nxyk;
		i01 = levelk1+i*nxyk1; 
		for (j = 1; j < nykk1-1; j++){
			j0 = i0+2*j*nxkk;
			j01 = i01+j*nxkk1;
			for (h = 1; h < nxkk1-1; h++){
				k0 = j0+2*h;
				k1 = k0-nxkk;
				k2 = k0+nxkk;
				k3 = k0-nxyk;
				k4 = k0-nxyk-nxkk;
				k5 = k0+nxyk;
				k6 = k0+nxyk+nxkk;
				b[j01+h] = (r[k0]*2+r[k0-1]+r[k0+1]+r[k1]+r[k1-1]+r[k2]+r[k2+1]+r[k3]
							+r[k3-1]+r[k4]+r[k4-1]+r[k5]+r[k5+1]+r[k6]+r[k6+1])/4;						
							

				//b[j01+h] = (r[k0]*2+r[k0-1]+r[k0+1]+r[k0-nxkk]+r[k0-nxkk-1]+
				//			r[k0+nxkk]+r[k0+nxkk+1]+r[k0-nxyk]+r[k0-nxyk-1]+
				//			r[k0-nxyk-nxkk]+r[k0-nxyk-nxkk-1]+r[k0+nxyk]+
				//			r[k0+nxyk+1]+r[k0+nxyk+nxkk]+r[k0+nxyk+nxkk+1])/4;


				//b[level[k+1]+i*nxk[k+1]*nyk[k+1]+j*nxk[k+1]+h] = 
				//	(r[level[k]+2*i*nxk[k]*nyk[k]+2*j*nxk[k]+2*h]*2+
				//	 r[level[k]+2*i*nxk[k]*nyk[k]+2*j*nxk[k]+2*h-1]+
				//	 r[level[k]+2*i*nxk[k]*nyk[k]+2*j*nxk[k]+2*h+1]+
				//	 r[level[k]+2*i*nxk[k]*nyk[k]+(2*j-1)*nxk[k]+2*h]+
				//	 r[level[k]+2*i*nxk[k]*nyk[k]+(2*j-1)*nxk[k]+2*h-1]+
				//	 r[level[k]+2*i*nxk[k]*nyk[k]+(2*j+1)*nxk[k]+2*h]+
				//	 r[level[k]+2*i*nxk[k]*nyk[k]+(2*j+1)*nxk[k]+2*h+1]+
				//	 r[level[k]+(2*i-1)*nxk[k]*nyk[k]+2*j*nxk[k]+2*h]+
				//	 r[level[k]+(2*i-1)*nxk[k]*nyk[k]+2*j*nxk[k]+2*h-1]+
				//	 r[level[k]+(2*i-1)*nxk[k]*nyk[k]+(2*j-1)*nxk[k]+2*h]+
				//	 r[level[k]+(2*i-1)*nxk[k]*nyk[k]+(2*j-1)*nxk[k]+2*h-1]+
				//	 r[level[k]+(2*i+1)*nxk[k]*nyk[k]+2*j*nxk[k]+2*h]+
				//	 r[level[k]+(2*i+1)*nxk[k]*nyk[k]+2*j*nxk[k]+2*h+1]+
				//	 r[level[k]+(2*i+1)*nxk[k]*nyk[k]+(2*j+1)*nxk[k]+2*h]+
				//	 r[level[k]+(2*i+1)*nxk[k]*nyk[k]+(2*j+1)*nxk[k]+2*h+1])/4;
			}
		}
	}
}

static finergridinterpolation2d(REAL *u,
							  INT *level, 
							  INT k, 
							  INT *nxk, 
							  INT *nyk)
{
	INT i,j;
//	int nthreads = 1;

	for (i = 0; i < nyk[k]-1; i++) {
//	#pragma omp parallel for private(j) num_threads(nthreads)
        for (j = 0; j < nxk[k]-1; j++) {
			u[level[k-1]+2*(j+i*nxk[k-1])] += u[level[k]+j+i*nxk[k]];
            u[level[k-1]+2*(j+i*nxk[k-1])+1] += (u[level[k]+j+i*nxk[k]]+
                                                 u[level[k]+j+1+i*nxk[k]])/2;
            u[level[k-1]+(2*i+1)*nxk[k-1]+2*j] += (u[level[k]+i*nxk[k]+j]+u[level[k]+(i+1)*nxk[k]+j])/2;
			u[level[k-1]+(2*i+1)*nxk[k-1]+2*j+1] += (u[level[k]+i*nxk[k]+j]+u[level[k]+(i+1)*nxk[k]+1+j])/2;
		}
	}
}

static finergridinterpolation3d(REAL *u, 
							  INT *level, 
							  INT k, 
							  INT *nxk, 
							  INT *nyk, 
							  INT *nzk)
{
	INT i, j, h;
	int i0,j0,j1,j2,j3,k0,k1,k2,k3,k4,k5,k6;
	int i01,j01,j11,j21,j31,k01,k11,k21,k31;
    const int levelk = level[k];
	const int levelk1 = level[k-1];
	const int nxkk = nxk[k];
	const int nxkk1 = nxk[k-1];
	const int nykk = nyk[k];
	const int nykk1 = nyk[k-1];
	const int nzkk = nzk[k];
	const int nzkk1 = nzk[k-1];
	// OpenMP settings
	//INT nthreads = omp_get_thread_num();
//	INT nthreads = 1;

	int nxyk = nxkk*nykk;
	int nxyk1 = nxkk1*nykk1;
//#pragma omp parallel for num_threads(nthreads)
    for (i = 0; i < nzkk-1; i++) {
		i0 = levelk+i*nxkk*nykk;
		i01 = levelk1+2*i*nxkk1*nykk1; 
        for (j = 0; j < nykk-1; j++) {
			j0 = i0+j*nxkk;
			j01 = i01 + 2*j*nxkk1;
			j1 = j0+nxkk;
			j11 = j01 + nxkk1;
			j2 = j0 + nxyk;
			j21 = j01+nxyk1;
			j3 = j2+nxkk;
			j31 = j21+nxkk1;
			//level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h
			//level[k-1]+(2*i+1)*nxk[k-1]*nyk[k-1]+2*j*nxk[k-1]+2*h
    //      u[level[k-1]+(2*i+1)*nxk[k-1]*nyk[k-1]+(2*j+1)*nxk[k-1]+2*h] += 
    //                 (u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //                  u[level[k]+(i+1)*nxk[k]*nyk[k]+(j+1)*nxk[k]+h])/2;

            for (h = 0; h < nxkk-1; h++) {
				k01 = j01+2*h;
				k0 = j0+h;
				u[k01] += u[k0];
				u[k01+1] += (u[k0]+u[k0+1])/2;
/*                 u[level[k-1]+2*i*nxk[k-1]*nyk[k-1]+2*j*nxk[k-1]+2*h] += 
			 	     u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h];
                 u[level[k-1]+2*i*nxk[k-1]*nyk[k-1]+2*j*nxk[k-1]+2*h+1] += 
                     (u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
                      u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h+1])/2;
 */             
				k11 = j11+2*h;
				k1 = j1+h;
				u[k11] +=(u[k0]+u[k1])/2; 
				u[k11+1] +=(u[k0]+u[k1+1])/2;
				//u[level[k-1]+2*i*nxk[k-1]*nyk[k-1]+(2*j+1)*nxk[k-1]+2*h] += 
    //                 (u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //                  u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h])/2;
    //             u[level[k-1]+2*i*nxk[k-1]*nyk[k-1]+(2*j+1)*nxk[k-1]+2*h+1] += 
    //                 (u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //                  u[level[k]+i*nxk[k]*nyk[k]+(j+1)*nxk[k]+h+1])/2;
				k21 = j21+2*h;
				k2 = j2+h;
				u[k21] += (u[k0]+u[k2])/2;
				u[k21+1] += (u[k0]+u[k2+1])/2;
/*                 u[level[k-1]+(2*i+1)*nxk[k-1]*nyk[k-1]+2*j*nxk[k-1]+2*h] += 
                     (u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+    
                      u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h])/2;    
                 u[level[k-1]+(2*i+1)*nxk[k-1]*nyk[k-1]+2*j*nxk[k-1]+2*h+1] += 
                     (u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+       
                      u[level[k]+(i+1)*nxk[k]*nyk[k]+j*nxk[k]+h+1])/2;    
 */             
				k31 = j31+2*h;
				k3 = j3+h;
				u[k31] += (u[k0]+u[k3])/2;
				u[k31+1] += (u[k0]+u[k3+1])/2;
				//u[level[k-1]+(2*i+1)*nxk[k-1]*nyk[k-1]+(2*j+1)*nxk[k-1]+2*h] += 
    //                 (u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
    //                  u[level[k]+(i+1)*nxk[k]*nyk[k]+(j+1)*nxk[k]+h])/2;
    //             u[level[k-1]+(2*i+1)*nxk[k-1]*nyk[k-1]+(2*j+1)*nxk[k-1]+2*h+1] += 
    //                 (u[level[k]+i*nxk[k]*nyk[k]+j*nxk[k]+h]+
				//	  u[level[k]+(i+1)*nxk[k]*nyk[k]+(j+1)*nxk[k]+h+1])/2;
             }
		}
	}
}

