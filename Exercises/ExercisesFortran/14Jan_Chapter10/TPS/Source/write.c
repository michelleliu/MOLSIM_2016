#include "tps.h"

/* ================================================================== */
void write_restart(){
    
  int iat, istep;
  FILE *fpout;
  
  if(lrstrtout==0) return;
  printf("  Writing last path to restart file\n");
  
  fpout = fopen(RESTARTOUT,"w");
  if(fpout==NULL) die("Failed to open restart output file");
  
  fprintf(fpout,"%d %d\n",NATOMS,npathlength);
  for(istep=1;istep<npathlength;istep++){
    for(iat=0;iat<NATOMS;iat++){
      fprintf(fpout,"%20.10lf%20.10lf%20.10lf%20.10lf\n",
        px[istep][iat],py[istep][iat],vx[istep][iat],vy[istep][iat]);
    }
  }
  fclose(fpout);
  return;
}

/* ================================================================== */
void write_frame(FILE *fpout, double xpos[NATOMS], double ypos[NATOMS],
		 int ipath, int istep){

  int iat;

  fprintf(fpout," %d\n Path %d: Time = %12.6lf step = %d\n",NATOMS,ipath,istep*dt,istep);
#ifdef PBC
  for(iat=0;iat<2;iat++){
    fprintf(fpout,"O  %12.6lf%12.6lf 0.0\n",xpos[iat]-LBOX*rint(xpos[iat]/LBOX),
      ypos[iat]-LBOX*rint(ypos[iat]/LBOX));
  }  for(iat=2;iat<NATOMS;iat++){
    fprintf(fpout,"Ar  %12.6lf%12.6lf 0.0\n",xpos[iat]-LBOX*rint(xpos[iat]/LBOX),
      ypos[iat]-LBOX*rint(ypos[iat]/LBOX));
  }
#else
  for(iat=0;iat<2;iat++){
    fprintf(fpout,"O  %12.6lf%12.6lf 0.0\n",xpos[iat],ypos[iat]);
  }  for(iat=2;iat<NATOMS;iat++){
    fprintf(fpout,"Ar %12.6lf%12.6lf 0.0\n",xpos[iat],ypos[iat]);
  }
#endif  
   return;
}

/* ================================================================== */
void write_lambda(FILE *fpout, double lambda, int ipath, int istep){

  fprintf(fpout,"%d %12.6lf %lf\n",ipath,istep*dt,lambda);
  return;
}

/* ================================================================== */
void write_energy(FILE *fpout, double vpot, double ekin,int ipath,int istep){
  
  fprintf(fpout,"%d %12.6lf %lf %lf %lf\n",ipath,istep*dt,ekin+vpot,ekin,vpot);
  return;
}

/* ================================================================== */
/* write histogram                                                    */
void write_histogram(){
  int ibin;
  double binsize, dr;

  binsize = MAXLAMBDA / (double) NBINS;
  for(ibin=0;ibin<NBINS;ibin++){
    dr = (0.5 + ibin) * binsize;
    fprintf(fpout_d," %12.6lf %20.12lf\n",dr,lambda_distr[ibin]);
  }
  return;
}

/* ================================================================== */
void write_output(double (*px)[NATOMS],double (*py)[NATOMS], double *vpot,
  double *ekin, double *lambda,int ipath){
  
  int istep;
  
  /* write output */
  for(istep=1;istep<npathlength;istep++){
    if(istep%nprintframe == 0){
      write_energy(fpout_e,vpot[istep],ekin[istep],ipath,istep); 
      write_frame(fpout_p,px[istep],py[istep],ipath,istep);
      write_lambda(fpout_l,lambda[istep],ipath,istep);
    }
  }

  /* write empty lines between blocks for plotting with gnuplot */
  fprintf(fpout_e,"\n\n");
  fprintf(fpout_l,"\n\n");
  
  /* flush */
  fflush(fpout_e);
  fflush(fpout_p);
  fflush(fpout_l);
  return;
}

