#include "tps.h"

/* ================================================================== */
/* integrate equations of motion to create (transition) path          */
/*   - using velocity verlet integrator                               */
/*   - "isign" parameter sets forward or backward in time integration */
int integrate(
		double *vpot,double *ekin,double *lambda,
    int istart, int nsteps, int isign, int icheckstate){

  int istep, iat, iframe, icheck;
  double fx[NATOMS], fy[NATOMS];
  double delta;
  double (*xpos)[NATOMS] = &px_try[istart];
  double (*ypos)[NATOMS] = &py_try[istart];
  double (*xvel)[NATOMS] = &vx_try[istart];
  double (*yvel)[NATOMS] = &vy_try[istart];
  
  iframe = 0;
  icheck = -1;
  delta = (double)isign*dt;
  ekin[0] = get_ekin(xvel[0],yvel[0]);
  force(xpos[0],ypos[0],fx,fy,&vpot[0]);
  lambda[0] = compute_lambda(xpos[0][0],ypos[0][0],xpos[0][1],ypos[0][1]);

  for(istep=1;istep<nsteps;istep++){
    iframe += isign;
    /* velocity verlet step */
    for(iat=0;iat<NATOMS;iat++){
      xvel[iframe][iat] = xvel[iframe-isign][iat] + 0.5 * delta * fx[iat];
      yvel[iframe][iat] = yvel[iframe-isign][iat] + 0.5 * delta * fy[iat];
    }
    for(iat=0;iat<NATOMS;iat++){
      xpos[iframe][iat] = xpos[iframe-isign][iat] + 0.5 * delta * xvel[iframe][iat];
      ypos[iframe][iat] = ypos[iframe-isign][iat] + 0.5 * delta * yvel[iframe][iat];
    }
#ifndef PBC
    apply_reflective_boundaries(xpos[iframe],ypos[iframe],xvel[iframe],yvel[iframe]);
#endif
    for(iat=0;iat<NATOMS;iat++){
      xpos[iframe][iat] = xpos[iframe][iat] + 0.5 * delta * xvel[iframe][iat];
      ypos[iframe][iat] = ypos[iframe][iat] + 0.5 * delta * yvel[iframe][iat];
    }
    force(xpos[iframe],ypos[iframe],fx,fy,&vpot[iframe]);
    for(iat=0;iat<NATOMS;iat++){
      xvel[iframe][iat] = xvel[iframe][iat] + 0.5 * delta * fx[iat];
      yvel[iframe][iat] = yvel[iframe][iat] + 0.5 * delta * fy[iat];
    }
    /* compute kinetic energy */
    ekin[iframe] = get_ekin(xvel[iframe],yvel[iframe]);

    lambda[iframe] = compute_lambda(xpos[iframe][0],ypos[iframe][0],
				    xpos[iframe][1],ypos[iframe][1]);
    if(lambda[iframe] < leftbound){
      icheck = 0;
      if(icheckstate==COMMITTOR) return icheck;
    }else if(lambda[iframe] > rightbound){
      icheck = 1;
      if(icheckstate==COMMITTOR) return icheck;
    }else{
      icheck = -1;
    }
  }
  return icheck;
}

