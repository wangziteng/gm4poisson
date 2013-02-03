/**
 * \brief For Fortran compatibilty 
 */
#define SHORT            short  /**< short integer type */
#define INT              int    /**< regular integer type */
#define LONG             long   /**< long integer type */ 
#define REAL             double /**< float type */

/** 
 * \struct ivector
 * \brief Vector with n entries of INT type.
 */
typedef struct ivector{
	
    //! number of rows
	INT row;
    //! actual vector entries
	INT *val;
	
} ivector; /**< Vector of INT type */

/** 
 * \struct dvector
 * \brief Vector with n entries of REAL type.
 */
typedef struct dvector{
	
    //! number of rows
	INT row;
    //! actual vector entries
	REAL *val;
	
} dvector; /**< Vector of REAL type */

void multigriditeration1d(REAL *u,
                          REAL *b,
                          INT *level,
                          INT maxlevel);

void multigriditeration2d(REAL *u,
                          REAL *b,
                          INT *level,
                          INT maxlevel,
                          INT nx,
                          INT ny);

void multigriditeration3d(REAL *u,
						   REAL *b,
                           INT *level,
                           INT maxlevel,
                           INT nx,
                           INT ny,
                           INT nz);

void multigrid_full_3d(REAL *u,
                       REAL *b,
                       INT *level,
                       INT maxlevel,
			   		   INT nx,
					   INT ny,
					   INT nz);

void gsiteration3d(REAL *u,
				   REAL *b, 
				   INT *level, 
				   INT k, 
				   INT maxlevel, 
				   INT *nxk, 
				   INT *nyk, 
				   INT *nzk);

void gsiteration_2color_2d(REAL *u,
						   REAL *b,
						   INT *level,
                           INT k,
                           INT maxlevel,
                           INT *nxk, 
                           INT *nyk);

void gsiteration_2color_3d(REAL *u,
                           REAL *b,
                           INT *level,
                           INT k,
                           INT maxlevel,
                           INT *nxk, 
                           INT *nyk, 
                           INT *nzk);

void coarsergridrestriction2d(REAL *u, 
							  REAL *b, 
							  REAL *r, 
							  INT *level, 
							  INT k, 
							  INT *nxk, 
							  INT *nyk);

void coarsergrid7pointrestriction2d(REAL *b, 
									REAL *r, 
									INT *level, 
									INT k, 
									INT *nxk, 
									INT *nyk);

void coarsergridrestriction3d(REAL *u, 
							  REAL *b,
							  REAL *r, 
							  INT *level, 
							  INT k, 
							  INT maxlevel,
							  INT *nxk, 
							  INT *nyk, 
							  INT *nzk);

void coarsergrid7pointrestriction3d(REAL *b, 
									REAL *r, 
									INT *level, 
									INT k, 
									INT *nxk, 
									INT *nyk,
									INT *nzk);

void finergridinterpolation2d(REAL *u, 
							  INT *level, 
							  INT k, 
							  INT *nxk, 
							  INT *nyk);

void finergridinterpolation3d(REAL *u, 
							  INT *level, 
							  INT k, 
							  INT *nxk, 
							  INT *nyk, 
							  INT *nzk);

void fullmultigrid_1d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
                      INT nx,
					  REAL rtol);

void fullmultigrid_2d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
                      INT nx,
                      INT ny,
					  REAL rtol);

void fullmultigrid_3d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
                      INT nx,
                      INT ny,
                      INT nz,
					  REAL rtol);

void vcycle_at_levelk_1d(REAL *u, 
					     REAL *b, 
						 REAL *r, 
						 INT *level, 
						 INT k,
						 INT maxlevel);

void vcycle_at_levelk_2d(REAL *u, 
					     REAL *b, 
						 REAL *r, 
						 INT *level, 
						 INT k, 
						 INT maxlevel, 
						 INT *nxk, 
						 INT *nyk);

void vcycle_at_levelk(REAL *u, 
					  REAL *b, 
					  REAL *r, 
					  INT *level, 
					  INT k, 
					  INT maxlevel, 
					  INT *nxk, 
					  INT *nyk, 
					  INT *nzk);

void pcg_2d(REAL *u,
            REAL *b,
            INT *level,
			INT maxlevel,
            INT nx,
            INT ny,
            REAL rtol,
			INT maxiteration);

void pcg_3d(REAL *u,
            REAL *b,
            INT *level,
			INT maxlevel,
            INT nx,
            INT ny,
            INT nz, 
            REAL rtol,
			INT maxiteration);

void multigrid_vcycle_1d(REAL *u,
                         REAL *b,
                         INT *level,
                         INT maxlevel);

void multigrid_vcycle_2d(REAL *u,
                         REAL *b,
                         INT *level,
                         INT maxlevel,
						 INT nx,
						 INT ny);

void multigrid_vcycle_3d(REAL *u,
                         REAL *b,
                         INT *level,
                         INT maxlevel,
						 INT nx,
						 INT ny,
						 INT nz);

void multigrid_full_1d(REAL *u,
                       REAL *b,
                       INT *level,
                       INT maxlevel,
			   		   INT nx);

void multigrid_full_2d(REAL *u,
                       REAL *b,
                       INT *level,
                       INT maxlevel,
			   		   INT nx,
					   INT ny);

void multigrid_full_3d(REAL *u,
                       REAL *b,
                       INT *level,
                       INT maxlevel,
			   		   INT nx,
					   INT ny,
					   INT nz);

void multigrid_pcg_2d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
	       	  		  INT nx,
					  INT ny);

void multigrid_pcg_3d(REAL *u,
                      REAL *b,
                      INT *level,
                      INT maxlevel,
	       	  		  INT nx,
				      INT ny,
					  INT nz);

void compute_r_1d(REAL *b,
                  REAL *u,
                  REAL *r,
                  INT   k,
                  INT *level);

void compute_r_2d(REAL *u,
                  REAL *b,
                  REAL *r,
                  INT k,
                  INT *level,
                  INT *nxk,
                  INT *nyk);

void compute_r_3d(REAL *u,
                  REAL *b,
                  REAL *r,
                  INT k,
                  INT *level,
                  INT *nxk,
                  INT *nyk,
                  INT *nzk);

REAL computenorm(REAL *u,
                 INT *level,
                 INT k);

void xequaly(REAL *x,
             REAL *y,
             INT *level,
             INT k);

void xequalypcz(REAL *x,
                REAL *y,
                REAL c,
                REAL *z,
                INT *level,
                INT k);

REAL compute_quad_norm_1d(REAL *u,
                          INT *level, 
                          INT k,
                          INT *nxk);

REAL compute_quad_norm_2d(REAL *u,
                          INT *level,
                          INT k,
                          INT *nxk,
                          INT *nyk);

REAL compute_quad_norm_3d(REAL *u,
                          INT *level,
                          INT k,
                          INT *nxk,
                          INT *nyk,
                          INT *nzk);

void xequalay_1d(REAL *x,
                 REAL *y,
                 INT *level,
                 INT k);

void xequalay_2d(REAL *x,
                 REAL *y,
                 INT *level,
                 INT k,
                 INT *nxk,
                 INT *nyk);

void xequalay_3d(REAL *x,
                 REAL *y,
                 INT *level,
                 INT k,
                 INT *nxk,
                 INT *nyk,
                 INT *nzk);

REAL innerproductxy(REAL *x, 
					REAL *y, 
					INT *level, 
					INT k);

REAL F(REAL x);