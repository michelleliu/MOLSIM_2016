#include<stdlib.h>
#include<stdio.h>
#include "math.h"
#include "string.h"

#ifndef DANGLEMAX
#define DANGLEMAX 2.0  /* a random rotation (in degrees) of the velocity vector is taken from [-DANGLEMAX,DANGELMAX] */
#endif

#define MAXLINE 128
#define MAXFRAME 10000
#define NATOMS 15
#define RBOX 3.0
#define MD 0
#define TPS 1
#define COMMITTOR 2
#define NBINS 200
#define MAXLAMBDA 2.0

#define PBC_OFF    /* with PBC: use square box with pbc; without: use cicle with hard wall */
#define LBOX 5.0

#define XYZOUTPUT "vmd_movie.xyz"
#define ENEROUTPUT "energies.out"
#define LAMBDAOUTPUT "lambda.out"
#define DISTROUTPUT "distribution.out"
#define RESTARTOUT "restart.out"
#define RWCA pow(2.0,1.0/6.0)

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

#ifndef M_SQRT2
#define M_SQRT2         1.41421356237309504880
#endif

int npaths;           /* number of paths to generate */
int npathlength;      /* number of frames in paths to generate */
int nlengthrestart;   /* number of frames in restartfile */
int nprintframe;      /* frame interval to write frames to file */
int nprintpath;       /* path interval for writing paths to file */
int maxtrail;         /* max trails to find accepted path in TPS */
int ishiftshoot;      /* if 1 then shift else if 2 then shoot else do not perturb */
int istartframe;      /* if larger than 0 shoot from this frame of a restart trajectory */
int imethod;          /* imethod is MD or TPS or COMMITTOR */
int lrstrtout;        /* if 1 then write restart.out */
double dt;            /* MD timestep */
double leftbound;
double rightbound;
double height;
double width;
double k_restraint;
double req_restraint;
double lambda_distr[NBINS];
double epsilon;
double *lambda;
double *ekin, *vpot;
double (*px)[NATOMS],(*py)[NATOMS],(*px_try)[NATOMS],(*py_try)[NATOMS];
double (*vx)[NATOMS],(*vy)[NATOMS],(*vx_try)[NATOMS],(*vy_try)[NATOMS];
char restartfilename[MAXLINE];
FILE *fpout_p, *fpout_e, *fpout_l, *fpout_d;

/* initialize */
void initialize();
void init_pos(double xpos[NATOMS],double ypos[NATOMS]);
void init_vel(double xvel[NATOMS],double yvel[NATOMS]);

/* utils.c */
double myrandom(void);
double gauss1(void);
void usage(char *prog);
void die(char *line);
double get_ekin(double xvel[NATOMS], double yvel[NATOMS]);
void store_trailpath();
void scale_vel(double vpot,double etot,double xvel[NATOMS],double yvel[NATOMS]);
void apply_reflective_boundaries(double xpos[NATOMS],double ypos[NATOMS],
				 double xvel[NATOMS],double yvel[NATOMS]);
double compute_lambda(double px1,double py1, double px2, double py2);
void compute_lambda_distr(double *lambda,int npathlength);

/* force.c */
void force(double xpos[NATOMS], double ypos[NATOMS],
	   double xforce[NATOMS], double yforce[NATOMS], double *vpot);

/* integrate.c */
int integrate(
		double *vpot,double *ekin,double *lambda,
    int istart, int nsteps, int isign, int icheckstate);

/* read.c */
void read_input(char *filename,double *etot,int *lrestart);
void read_restart(double xpos[NATOMS],double ypos[NATOMS],double xvel[NATOMS],double yvel[NATOMS]);

/* write.c */
void write_restart();
void write_frame(FILE *fpout, double xpos[NATOMS], double ypos[NATOMS],
    int ipath, int istep);
void write_lambda(FILE *fpout, double lambda, int ipath, int istep);
void write_energy(FILE *fpout, double vpot, double ekin,int ipath,int istep);
void write_histogram();
void write_output(double (*px)[NATOMS],double (*py)[NATOMS], double *vpot,double *ekin,
  double *lambda,int ipath);
