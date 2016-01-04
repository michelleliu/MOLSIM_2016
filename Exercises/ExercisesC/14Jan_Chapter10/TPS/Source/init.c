#include"tps.h"

/* ================================================================== */
void initialize(){
  int nlen, ibin;

  /* allocate some memory to store trajectories */
  nlen = npathlength+1;
  ekin = (double *)malloc(nlen*sizeof(double));
  vpot = (double *)malloc(nlen*sizeof(double));
  lambda = (double *)malloc(nlen*sizeof(double));
  px = (double (*)[NATOMS])malloc(nlen*sizeof(double [NATOMS]));
  py = (double (*)[NATOMS])malloc(nlen*sizeof(double [NATOMS]));
  vx = (double (*)[NATOMS])malloc(nlen*sizeof(double [NATOMS]));
  vy = (double (*)[NATOMS])malloc(nlen*sizeof(double [NATOMS]));
  px_try = (double (*)[NATOMS])malloc(nlen*sizeof(double [NATOMS]));
  py_try = (double (*)[NATOMS])malloc(nlen*sizeof(double [NATOMS]));
  vx_try = (double (*)[NATOMS])malloc(nlen*sizeof(double [NATOMS]));
  vy_try = (double (*)[NATOMS])malloc(nlen*sizeof(double [NATOMS]));
  
  /* open output files */
  fpout_p = fopen(XYZOUTPUT,"w");
  fpout_e = fopen(ENEROUTPUT,"w");
  fpout_l = fopen(LAMBDAOUTPUT,"w");
  fpout_d = fopen(DISTROUTPUT,"w");
  fprintf(fpout_e,"# time Etot Ekin Vpot\n");
  fprintf(fpout_l,"# time lambda\n");
  fprintf(fpout_d,"# lambda P(lambda)\n");
  
  /* initialize histogram to zero */
  for(ibin=0;ibin<NBINS;ibin++) lambda_distr[ibin] = 0.0;
  
  return;
}

/* ================================================================== */
/* place particles on circle */
void init_pos(double xpos[NATOMS],double ypos[NATOMS]){
  
  int iat;
#ifdef PBC
  int ngrid;
  double xgrid, ygrid, dgrid;
#else
  double dangle, angle;
#endif

#ifdef PBC
  /* place particles on a lattice */
  ngrid = (int)sqrt((double)NATOMS) + 1;
  dgrid = LBOX/(double)ngrid;
  xgrid = ygrid = 0.0;
  for(iat=0;iat<NATOMS;iat++){
    xpos[iat] = xgrid;
    ypos[iat] = ygrid;
    xgrid = xgrid + dgrid;
    if(xgrid >= LBOX){
      xgrid = 0.0;
      ygrid = ygrid + dgrid;
    }
  }
#else
  dangle = 2.0 * M_PI / (double)NATOMS;
  angle = 0.0;
  for(iat=0;iat<NATOMS;iat++){
    xpos[iat] = RBOX * cos(angle);
    ypos[iat] = RBOX * sin(angle);
    angle += dangle;
  }
#endif    
  return;
}


/* ================================================================== */
/* draw random velocities from a Boltzman distribution */
void init_vel(double xvel[NATOMS],double yvel[NATOMS]){
  
  int iat;

  /* draw velocities */
  for(iat=0;iat<NATOMS;iat++){
    xvel[iat] = gauss1();
    yvel[iat] = gauss1();
  }
  return;
}
