#include "tps.h"

/* ================================================================== */
void force(double xpos[NATOMS], double ypos[NATOMS],
	   double xforce[NATOMS], double yforce[NATOMS], double *vpot){
  
  int iat, jat;
  double dx, dy, dr, rwca, pot, deriv, rsix;
  double ff, gg, lambda;


  /* init */
  rwca = RWCA;
  pot = 0.0;
  for(iat=0;iat<NATOMS;iat++){
    xforce[iat] = yforce[iat] = 0.0;
  }

  /* sum WCA pair potentials (skipping pair 0 and 1) */
  for(iat=0;iat<NATOMS-1;iat++){
    for(jat=iat+1;jat<NATOMS;jat++){
      if((iat==0)&&(jat==1)) continue;
      dx = xpos[iat] - xpos[jat];
      dy = ypos[iat] - ypos[jat];
#ifdef PBC
      dx  -= LBOX*((int)(2.0*dx/LBOX)-(int)(dx/LBOX));
      dy  -= LBOX*((int)(2.0*dy/LBOX)-(int)(dy/LBOX));
//      dx -= LBOX*rint(dx/LBOX);
//      dy -= LBOX*rint(dy/LBOX);
#endif
      dr = sqrt(dx*dx+dy*dy);
      if(dr <= rwca ){
	rsix = pow(1.0/dr,6.0);
	pot += epsilon*(1.0 + 4.0*(rsix*rsix - rsix));
	deriv = 2.0 * epsilon * 4.0 * (6.0*rsix - 12.0*rsix*rsix) / (dr*dr);
	xforce[iat] -= dx * deriv;
	yforce[iat] -= dy * deriv;
	xforce[jat] += dx * deriv;
	yforce[jat] += dy * deriv;
      }
    }
  }

  /* double well potential between particles 1 and 2 */
  iat = 0;
  jat = 1;
  dx = xpos[iat] - xpos[jat];
  dy = ypos[iat] - ypos[jat];
#ifdef PBC
      dx  -= LBOX*( (int)(2.0*dx/LBOX)-(int)(dx/LBOX));
      dy  -= LBOX*( (int)(2.0*dy/LBOX)-(int)(dy/LBOX));
//      dx -= LBOX*rint(dx/LBOX);
//      dy -= LBOX*rint(dy/LBOX);
#endif
  lambda = sqrt(dx*dx+dy*dy);
  ff = (lambda - width - rwca)/width;
  gg = 1.0 - ff*ff;
  pot += height * gg *gg;
  deriv = -8.0 * height * gg * ff/width / lambda;
  xforce[iat] -= dx * deriv;
  yforce[iat] -= dy * deriv;
  xforce[jat] += dx * deriv;
  yforce[jat] += dy * deriv;

  /* harmonic restraint potential between particles 1 and 2 */
  ff = (lambda - req_restraint);
  gg = 0.5 * k_restraint * ff*ff;
  if( (k_restraint != 0.0) && ((gg > 1000.0)||(abs(ff)>0.5)) ){
    printf("  WARNING: RESTRAINT POTENTIAL OR RESTRAINT DISTANCE VERY LARGE (simulation may crash)\n");
    printf("  - lambda= %lf  Req= %lf  k= %lf  Vrestraint= %lf\n",lambda,req_restraint,k_restraint,gg);
  }
  pot += gg;
  deriv = k_restraint * ff / lambda;
  xforce[iat] -= dx * deriv;
  yforce[iat] -= dy * deriv;
  xforce[jat] += dx * deriv;
  yforce[jat] += dy * deriv;

  *vpot = 2.0*pot;
  return;
}

