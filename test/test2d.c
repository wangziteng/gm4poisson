#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "multigrid_functs.h"

int main(int argc, const char *argv[])
{
    INT nx,maxlevel,i,ny,j;
    INT *level;
    REAL *u,*b,h;

    // temporalily nx is 2^n
	maxlevel = 9;
    nx = pow(2,maxlevel);
	ny = pow(2,maxlevel);
	printf("%d\n",nx);
    
	level = (INT *)malloc((maxlevel+1)*sizeof(INT));
	for (i = 0; i < maxlevel+1; i++) {
        level[i] = 0;
    }	
    level[0] = 0;
	level[1] = (nx+1)*(ny+1);
	u = (REAL *)malloc(3*nx*ny*sizeof(REAL));
    b = (REAL *)malloc(3*nx*ny*sizeof(REAL));
	h = 1.0/nx;
    for (i = 0; i < 3*nx*ny; i++) {
        u[i] = 0.0;
        b[i] = 0.0;
    }
	for (i = 1; i < ny; i++){
		for (j = 1; j < nx; j++){
			b[i*(nx+1)+j] = h*h;
		}
	}

    multigrid_vcycle_2d(u, b, level, maxlevel, nx, ny);
	//multigrid_pcg_2d(u, b, level, maxlevel, nx, ny);
	//multigrid_full_2d(u, b, level, maxlevel, nx, ny);

    free(b);
    free(u);

    return 0;
}
