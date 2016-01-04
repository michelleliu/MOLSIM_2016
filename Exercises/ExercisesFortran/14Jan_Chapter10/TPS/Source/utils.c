
#include"tps.h"
double dseed = 428769461.0; /* random number seed */

/* ================================================================== */
/* draw random number between 0 and 1 */
double myrandom(void){

  double random;
  const double d2p31m=2147483647.0, d2p31=2147483648.0, d7p5=16807.0;

  dseed=fmod(d7p5*dseed,d2p31m);
  random=dseed/d2p31;
  return random;
}

/* ================================================================== */
/* draw random number from Gaussian distribution */
double gauss1(void)
{
  double phi,al;
  phi = 2.*M_PI*myrandom();
  al = sqrt(-log(myrandom()))*M_SQRT2;
  return (al*cos(phi));
}

/* ================================================================== */
void usage(char *prog){
  printf("  Program for transition path sampling of model system\n");
  printf("  Usage: %s <inputfile>\n  Program stopped\n",prog);
  exit(1);
} 

/* error message and close */
/* ================================================================== */
void die(char *line){
  printf("  ERROR: %s\n  Program stopped\n",line);
  fclose(fpout_p);
  fclose(fpout_e);
  fclose(fpout_l);
  exit(1);
}
  
/* ================================================================== */
double get_ekin(double xvel[NATOMS], double yvel[NATOMS]){

    int iat;
    double energy;

    energy = 0.0;
    for(iat=0;iat<NATOMS;iat++){
      energy += xvel[iat] * xvel[iat] + yvel[iat] * yvel[iat];
    }
    energy *= 0.5;

    return energy;
}


/* ================================================================== */
/* copy trail path to "old" path                                      */
void store_trailpath(){
  double (*ptmp)[NATOMS];
  
  /* swap array pointers */
  ptmp = px; px = px_try; px_try = ptmp;
  ptmp = py; py = py_try; py_try = ptmp;
  ptmp = vx; vx = vx_try; vx_try = ptmp;
  ptmp = vy; vy = vy_try; vy_try = ptmp;
}

/* ================================================================== */
/* scale velocities to obtain target total energy given vpot */
void scale_vel(double vpot, double etot, 
	       double xvel[NATOMS],double yvel[NATOMS]){
  
  int iat;
  double ekin, scalefact;

  /* measure ekin */
  ekin = get_ekin(xvel,yvel);

  /* scale velocities to aim for target energy */
  if(etot < vpot) die("Failed to scale Ekin to obtain target energy (Potential energy to high)");
  scalefact = sqrt( (etot-vpot) / ekin );
  for(iat=0;iat<NATOMS;iat++){
    xvel[iat] *= scalefact;
    yvel[iat] *= scalefact;
  }
  ekin = get_ekin(xvel,yvel);

  return;
}

/* ================================================================== */
void apply_reflective_boundaries(double xpos[NATOMS],double ypos[NATOMS],
				 double xvel[NATOMS],double yvel[NATOMS]){
  int iat;
  double dr, dx, dy, proj, rbox2;

  rbox2 = RBOX*RBOX;
  /* check if atom is outside the spherical box */
  for(iat=0;iat<NATOMS;iat++){
    dr = xpos[iat]*xpos[iat] + ypos[iat]*ypos[iat];
    if(dr >= rbox2){
      /* reflect velocity */
      proj = (xvel[iat]*ypos[iat] - yvel[iat]*xpos[iat]) / dr;
      dx = ypos[iat] * proj;
      dy = -xpos[iat] * proj;
      xvel[iat] -= 2.0 * dx;
      yvel[iat] -= 2.0 * dy;
      xvel[iat] *= -1.0;
      yvel[iat] *= -1.0;
    }
  }
  return;
}

      
/* ================================================================== */
double compute_lambda(double px1,double py1, double px2, double py2){

  double dx, dy;

  /* compute distance between particles 1 and 2 */
  dx = px1 - px2;
  dy = py1 - py2;
  return sqrt(dx*dx+dy*dy);
}

/* ================================================================== */
/* make a histogram of lambda                                         */
void compute_lambda_distr(double *lambda,int npathlength){
  int iframe, ibin;
  double binsize;
    
  binsize = MAXLAMBDA / (double)NBINS;
  for(iframe=0;iframe<npathlength;iframe++){
    ibin = (int)(lambda[iframe] / binsize);
    if(ibin<NBINS) lambda_distr[ibin] += 1.0;
  }
  
  /* normalize */
  for(ibin=0; ibin<NBINS; ibin++) lambda_distr[ibin] /= ((double) npathlength * binsize);
  
  return;
}

