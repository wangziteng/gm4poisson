#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "multigrid_functs.h"

int main(int argc, const char *argv[])
{
    INT nx,maxlevel,i,ny,nz,j,k;
    INT *level;
    REAL *u,*b,h;

    // temporalily nx is 2^n
	maxlevel = 8;
    nx = pow(2,maxlevel);
	ny = pow(2,maxlevel);
	nz = pow(2,maxlevel);

	level = (INT *)malloc((maxlevel+2)*sizeof(INT));
    level[0] = 0;
	level[1] = (nx+1)*(ny+1)*(nz+1);
	h = 1.0/(nx);
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(nx/pow(2,i)+1)*(ny/pow(2,i)+1)*(nz/pow(2,i)+1);
    }
	level[maxlevel+1] = level[maxlevel]+1;
	u = (REAL *)malloc(level[maxlevel+1]*sizeof(REAL));
    b = (REAL *)malloc(level[maxlevel+1]*sizeof(REAL));
	// set the initial vectors
    for (i = 0; i < level[maxlevel+1]; i++) {
        u[i] = 0.0;
		b[i] = 0.0;
    }
	for (i = 1; i < nz; i++){
		for (j = 1; j < ny; j++){
			for (k = 1; k < nx; k++){
				b[k+j*(nx+1)+i*(nx+1)*(ny+1)] = h*h;
			}
		}
	}

    multigrid_vcycle_3d(u, b, level, maxlevel, nx, ny, nz);
	//multigrid_full_3d(u, b, level, maxlevel, nx, ny, nz);
	//multigrid_pcg_3d(u, b, level, maxlevel, nx, ny, nz);
    free(b);
    free(u);

    return 0;
}
