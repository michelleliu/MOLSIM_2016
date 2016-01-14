#include "tps.h"

/* ================================================================== */
void finish(){

  /* write last path to restart file  */
  write_restart();

  /* finalize and exit  */
  printf(" Simulation finished\n");
  fclose(fpout_p);
  fclose(fpout_e);
  fclose(fpout_l);
  fclose(fpout_d);
}

/* ================================================================== */
/* perform transition path sampling                                   */
void do_tps(){

  int ipath, ichk_f,ichk_b, nforward, nbackward, iframe, iat, istart;
  int ntrails, nshootmoves, nshiftmoves, nacceptshoot, nacceptshift;
  double angle;
  char line[MAXLINE];

  ipath = ntrails = nshootmoves = nshiftmoves = nacceptshoot = nacceptshift = 0;
  nforward = nbackward = istart = 0;

  /* main TPS loop */
  while(ipath < npaths){

    if(ishiftshoot == 1){
      /* try random shooting move */
      ishiftshoot = 2;
      nshootmoves++;
      /* choose random frame */
      iframe = (int) (npathlength * myrandom());
      /* make random rotation of all particle velocities */
      for(iat=0;iat<NATOMS;iat++){
        angle = 2.0 * DANGLEMAX * (myrandom() - 0.5);
        angle *= M_PI / 180.0;
        px_try[iframe][iat] = px[iframe][iat];
        py_try[iframe][iat] = py[iframe][iat];
        vx_try[iframe][iat] = vx[iframe][iat] * cos(angle) - vy[iframe][iat] * sin(angle) ;
        vy_try[iframe][iat] = vy[iframe][iat] * cos(angle) + vx[iframe][iat] * sin(angle);
      }
      istart = iframe;
      nforward = npathlength - iframe;
      nbackward = iframe;
    }else if(ishiftshoot == 2){
      /* try random shift move */
      ishiftshoot = 1; /* continue with shooting move hereafter */
      nshiftmoves++;
      istart = npathlength / 2 - 1;
      iframe = (int) (npathlength * myrandom());
      for(iat=0;iat<NATOMS;iat++){
        px_try[istart][iat] = px[iframe][iat];
        py_try[istart][iat] = py[iframe][iat];
        vx_try[istart][iat] = vx[iframe][iat];
        vy_try[istart][iat] = vy[iframe][iat];
      }
      nforward = npathlength / 2;
      nbackward = npathlength - nforward;
    }

    /* run new trial path forward and backward */
    ichk_b = integrate(&vpot[istart],&ekin[istart],&lambda[istart],istart,nbackward,-1,TPS);
    ichk_f = integrate(&vpot[istart],&ekin[istart],&lambda[istart],istart,nforward,1,TPS);

    /* if trial path is accepted store path and write some output */
    if( (ichk_b == 0) && (ichk_f == 1) ){
      if(ishiftshoot==1) nacceptshift++;
      if(ishiftshoot==2) nacceptshoot++;
      ipath++;
      store_trailpath();

      /* write output */
      if(ipath%nprintpath == 0){
        printf("  Path %d accepted in %d trails. Acceptance Ratios: %6.3lf (shifting) %6.3lf (shooting)\n",
          ipath,ntrails,(double)nacceptshift/(double)nshiftmoves,(double)nacceptshoot/(double)nshootmoves);
        write_output(px,py,vpot,ekin,lambda,ipath);
      }
      ntrails = 0;
    }else{
      ntrails++;
      if(ntrails > maxtrail){
        sprintf(line," Failed to find acceptable path in %d trails; choose better initial path or increase MAX_TPS_TRAILS ",maxtrail);
        die(line);
      }
    }
  }
  return;
}

/* ================================================================== */
/* return committor value of this frame by shooting random pathways   */
double do_committor(int istart, double etot){
  int iat, ipath, icheck;
  double xforce[NATOMS], yforce[NATOMS];
  double p;

  p = 0.0;
  for(iat=0;iat<NATOMS;iat++){
    px_try[0][iat] = px[istart][iat];
    py_try[0][iat] = py[istart][iat];
  }
  force(px_try[0],py_try[0],xforce,yforce,vpot);

  for(ipath=0;ipath<npaths;ipath++){
    init_vel(vx_try[0],vy_try[0]);
    scale_vel(vpot[0],etot,vx_try[0],vy_try[0]);
    icheck = integrate(&vpot[0],&ekin[0],&lambda[0],0,npathlength,1,MD);
    p += (double)icheck;
    write_output(px_try,py_try,vpot,ekin,lambda,ipath);
  }
  return p/(double)npaths;
}
/* ================================================================== */

int main(int nvar, char **cvar){
  int lrestart;
  double etot, xforce[NATOMS], yforce[NATOMS];
  printf("DANGLEMAX: %f\n",DANGLEMAX);
  /* read input file */
  if(nvar != 2) usage(cvar[0]);
  read_input(cvar[1],&etot,&lrestart);

  /* initialize */
  initialize();

  /* get initial configuration */
  if(lrestart){
    /* either from restart file */
    read_restart(px_try[0],py_try[0],vx_try[0],vy_try[0]);
  }else{
    /* or generate path from particles on circle and random velocities */
    init_pos(px_try[0],py_try[0]);
    init_vel(vx_try[0],vy_try[0]);
  }

  /* scale velocites to get target total energy and generate (initial) trajectory */
  if( imethod != COMMITTOR ){
    force(px_try[0],py_try[0],xforce,yforce,&vpot[0]);
    if(k_restraint == 0.0){   /* don't try setting total energy when there is a restraint potential */
      scale_vel(vpot[0],etot,vx_try[0],vy_try[0]);
    }
    integrate(&vpot[0],&ekin[0],&lambda[0],0,npathlength,1,MD);
    store_trailpath();
    write_output(px,py,vpot,ekin,lambda,0);
  }else{
    store_trailpath();
  }

  /* if we're doing normal MD, we can stop here */
  if(imethod == MD){
    compute_lambda_distr(lambda,npathlength);
    write_histogram();
    finish();
    return 0;
  }else if(imethod == TPS){
    do_tps();
  }else if(imethod == COMMITTOR){
    lambda[0] = compute_lambda(px[0][0],py[0][0],px[0][1],py[0][1]);
    printf("  - lambda[frame %d] = %12.6lf\n",istartframe,lambda[0]);
    printf(" p = %g\n",do_committor(0,etot));
  }


  /* finish and exit */
  finish();
  return 0;
}
