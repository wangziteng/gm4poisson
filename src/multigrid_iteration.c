#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>

#include "multigrid_functs.h"

void multigrid_vcycle_1d(REAL *u,
                         REAL *b,
                         INT *level,
                         INT maxlevel)
{
    REAL atol = 1.0E-15;   
	REAL rtol = 1.0E-8;
	REAL error;
    REAL *r0;
    REAL norm_r0, norm_r;
    INT done, max_itr_num,count,i;

    r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	for(i=0;i<level[maxlevel];i++){
		r0[i] = 0;
	}
    compute_r_1d(b, u, r0, 0, level);
    norm_r0 = computenorm(r0, level, 0);
	printf("norm_r0=%lf\n",norm_r0);
    if (norm_r0 < atol) {
        free(r0);
//        return 0;
    }

	done = 0; max_itr_num = 15;count=0;

    while ((!done) && (count < max_itr_num)) {
        count++;
        multigriditeration1d(u, b, level, maxlevel);
        compute_r_1d(b, u, r0, 0, level);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
		printf("%f,%d\n",error,count);
        if (error < rtol || norm_r < atol) {
            done = 1;                
        }
    }    

    if (count > max_itr_num) {
        printf("Multigrid failed to converge.\n");
        exit(EXIT_FAILURE);
    }
    
    free(r0);
    
}

void multigrid_vcycle_2d(REAL *u,
                         REAL *b,
                         INT *level,
                         INT maxlevel,
						 INT nx,
						 INT ny)
{
    REAL atol = 1.0E-15;   
	REAL rtol = 1.0E-8;
	REAL error;
    REAL *r0;
    REAL norm_r0, norm_r,time0;
    INT done, max_itr_num,count,i,nxk0[1], nyk0[1];
	clock_t t1,t2;

		t1 = clock();
    r0 = (REAL *)malloc(3*nx*ny*sizeof(REAL));
	for(i=0;i<3*nx*ny;i++){
		r0[i] = 0;
	}
	
	nxk0[0] = nx+1;
	nyk0[0] = ny+1;
	compute_r_2d(u, b, r0, 0, level, nxk0, nyk0);
    norm_r0 = computenorm(r0, level, 0);
	printf("norm_r0 = %f\n",norm_r0);
    if (norm_r0 < atol) {
        free(r0);
//        return 0;
    }

	done = 0; max_itr_num = 50;count=0;


    while ((!done) && (count < max_itr_num)) {
        count++;
        multigriditeration2d(u, b, level, maxlevel, nx, ny);
         compute_r_2d(u, b, r0, 0, level, nxk0, nyk0);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
		printf("%f,%d\n",error,count);
        if (error < rtol) {
            done = 1;                
        }
    }    

    if (count > max_itr_num) {
        printf("Multigrid failed to converge.\n");
        exit(EXIT_FAILURE);
    }
    t2 = clock();
	time0 = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("%lf seconds\n", time0);
    free(r0);
    
}

void multigrid_vcycle_3d(REAL *u,
                         REAL *b,
                         INT *level,
                         INT maxlevel,
						 INT nx,
						 INT ny,
						 INT nz)
{
    REAL atol = 1.0E-15;   
	REAL rtol = 1.0E-8;
	REAL error;
    REAL *r0;
    REAL norm_r0, norm_r, time0;
    INT done, max_itr_num,count,i,nxk0[1], nyk0[1], nzk0[1];
	clock_t t1,t2;

	t1 = clock();

    r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	for(i=0;i<level[maxlevel];i++){
		r0[i] = 0;
	}
	nxk0[0] = nx+1;
	nyk0[0] = ny+1;
	nzk0[0] = nz+1;

	compute_r_3d(u, b, r0, 0, level, nxk0, nyk0, nzk0);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) {
        free(r0);
        return;
    }

	done = 0; max_itr_num = 50;count=0;

    while ((!done) && (count < max_itr_num)) {
        count++;
        multigriditeration3d(u, b, level, maxlevel, nx, ny, nz);
        compute_r_3d(u, b, r0, 0, level, nxk0, nyk0, nzk0);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
		printf("%f,%d\n",error,count);
        if (error < rtol) {
            done = 1;                
        }
    }    

    if (count > max_itr_num) {
        printf("Multigrid failed to converge.\n");
        exit(EXIT_FAILURE);
    }
    t2 = clock();
	time0 = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("%lf seconds\n", time0);
    free(r0);
    
}

void multigrid_full_3d(REAL *u,
                       REAL *b,
                       INT *level,
                       INT maxlevel,
			   		   INT nx,
					   INT ny,
					   INT nz)
{
    REAL atol = 1.0E-15;   
	REAL rtol = 1.0E-8;
	REAL error;
    REAL *r0;
    REAL norm_r0, norm_r;
    INT done, max_itr_num,count,i,nxk0[1], nyk0[1], nzk0[1];

    r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	for(i=0;i<level[maxlevel];i++){
		r0[i] = 0;
	}
	nxk0[0] = nx;
	nyk0[0] = ny;
	nzk0[0] = nz;
    compute_r_3d(u, b, r0, 0, level, nxk0, nyk0, nzk0);
    norm_r0 = computenorm(r0, level, 0);
	printf("norm_r0 = %f\n",norm_r0);
    if (norm_r0 < atol) {
        free(r0);
//        return 0;
    }

	done = 0; max_itr_num = 50;count=0;

    while ((!done) && (count < max_itr_num)) {
        count++;
        fullmultigrid_3d(u, b, level, maxlevel, nx, ny, nz, rtol);
        compute_r_3d(u, b, r0, 0, level, nxk0, nyk0, nzk0);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
		printf("%f,%d\n",error,count);
        if (error < rtol) {
            done = 1;                
        }
    }    

    if (count > max_itr_num) {
        printf("Multigrid failed to converge.\n");
        exit(EXIT_FAILURE);
    }
    
    free(r0);
    
}

void multigrid_pcg_3d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
	       	  		  INT nx,
				      INT ny,
					  INT nz)

{    
	REAL atol = 1.0E-15;   
	REAL rtol = 1.0E-8;
	REAL error;
    REAL *r0;
    REAL norm_r0, norm_r;
    INT done, max_itr_num,count,i,nxk0[1], nyk0[1], nzk0[1];

    r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	for(i=0;i<level[maxlevel];i++){
		r0[i] = 0;
	}
	nxk0[0] = nx;
	nyk0[0] = ny;
	nzk0[0] = nz;

    compute_r_3d(u, b, r0, 0, level, nxk0, nyk0, nzk0);
    norm_r0 = computenorm(r0, level, 0);
	printf("norm_r0 = %f\n",norm_r0);
    if (norm_r0 < atol) {
        free(r0);
//        return 0;
    }

	done = 0; max_itr_num = 50;count=0;

    while ((!done) && (count < max_itr_num)) {
        count++;
        pcg_3d(u, b, level, maxlevel, nx, ny, nz, rtol, max_itr_num);
        compute_r_3d(u, b, r0, 0, level, nxk0, nyk0, nzk0);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
		printf("%f,%d\n",error,count);
        if (error < rtol) {
            done = 1;                
        }
    }    
}

void multigrid_pcg_2d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
	       	  		  INT nx,
					  INT ny)
{    
	REAL atol = 1.0E-15;   
	REAL rtol = 1.0E-8;
	REAL error;
    REAL *r0;
    REAL norm_r0, norm_r;
    INT done, max_itr_num,count,i,nxk0[1], nyk0[1], nzk0[1];

    r0 = (REAL *)malloc(3*nx*ny*sizeof(REAL));
	for(i=0;i<3*nx*ny;i++){
		r0[i] = 0;
	}
	nxk0[0] = nx+1;
	nyk0[0] = ny+1;

    compute_r_2d(u, b, r0, 0, level, nxk0, nyk0);
    norm_r0 = computenorm(r0, level, 0);
	printf("norm_r0 = %f\n",norm_r0);
    if (norm_r0 < atol) {
        free(r0);
		return;
    }

	done = 0; max_itr_num = 50;count=0;

    while ((!done) && (count < max_itr_num)) {
        count++;
        pcg_2d(u, b, level, maxlevel, nx, ny, rtol, max_itr_num);
        compute_r_2d(u, b, r0, 0, level, nxk0, nyk0);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
		printf("%f,%d\n",error,count);
        if (error < rtol) {
            done = 1;                
        }
    }    
}

void multigrid_full_2d(REAL *u,
                       REAL *b,
                       INT *level,
                       INT maxlevel,
			   		   INT nx,
					   INT ny)
{
    REAL atol = 1.0E-15;   
	REAL rtol = 1.0E-8;
	REAL error;
    REAL *r0;
    REAL norm_r0, norm_r;
    INT done, max_itr_num,count,i,nxk0[1], nyk0[1], nzk0[1];

    r0 = (REAL *)malloc(3*nx*ny*sizeof(REAL));
	for(i=0;i<3*nx*ny;i++){
		r0[i] = 0;
	}
	nxk0[0] = nx+1;
	nyk0[0] = ny+1;
    compute_r_2d(u, b, r0, 0, level, nxk0, nyk0);
    norm_r0 = computenorm(r0, level, 0);
	printf("norm_r0 = %f\n",norm_r0);
    if (norm_r0 < atol) {
        free(r0);
        return;
    }

	done = 0; max_itr_num = 50;count=0;

    while ((!done) && (count < max_itr_num)) {
        count++;
        fullmultigrid_2d(u, b, level, maxlevel, nx, ny, rtol);
        compute_r_2d(u, b, r0, 0, level, nxk0, nyk0);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
		printf("%f,%d\n",error,count);
        if (error < rtol) {
            done = 1;                
        }
    }
}

void multigrid_full_1d(REAL *u,
                       REAL *b,
                       INT *level,
                       INT maxlevel,
			   		   INT nx)
{
    REAL atol = 1.0E-15;   
	REAL rtol = 1.0E-8;
	REAL error;
    REAL *r0;
    REAL norm_r0, norm_r;
    INT done, max_itr_num,count,i,nxk0[1] ;

    r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	for(i=0;i<level[maxlevel];i++){
		r0[i] = 0;
	}
	nxk0[0] = nx;
    compute_r_1d(u, b, r0, 0, level);
    norm_r0 = computenorm(r0, level, 0);
	printf("norm_r0 = %f\n",norm_r0);
    if (norm_r0 < atol) {
        free(r0);
        return;
    }

	done = 0; max_itr_num = 50;count=0;

    while ((!done) && (count < max_itr_num)) {
        count++;
        fullmultigrid_1d(u, b, level, maxlevel, nx, rtol);
        compute_r_1d(u, b, r0, 0, level);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
		printf("%f,%d\n",error,count);
        if (error < rtol) {
            done = 1;                
        }
    }
}