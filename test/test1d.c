#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "malloc.h"

#include "multigrid_functs.h"

int main(int argc, const char *argv[])
{
    INT nx,maxlevel,i;
    INT *level;
    REAL *u,*b,*x,h;

    // temporalily nx is 2^n
	maxlevel = 18;
    nx = pow(2,maxlevel);

	level = (INT *)malloc((maxlevel+1)*sizeof(INT));
    for (i = 0; i < maxlevel+1; i++) {
        level[i] = 0;
    }	
    level[0] = 0;
	level[1] = nx;
	h = 1.0/nx;
    for (i = 1; i < maxlevel; i++) {
        level[i+1] += level[i]*3/2-level[i-1]/2;
    }
	u = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	x = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	x[0] = 0.0;
	for (i = 1; i < level[1]; i++){
		x[i] = x[i-1] + h;
	}
    for (i = 0; i < level[maxlevel]; i++) {
        u[i] = 0.0;
        b[i] = F(x[i])*h*h*pow(2,15);
    }
	b[0] = 0.0;

    multigrid_vcycle_1d(u, b, level, maxlevel);
	//multigrid_full_1d(u, b, level, maxlevel, nx);
    free(b);
    free(u);

    return 0;
}
