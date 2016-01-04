#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

#define MAXIMUM_NUMBER_OF_BINS 500

// samples the radial distribution function
void Sample(int Ichoise)
{
  int i,j;
  static double Distribution[MAXIMUM_NUMBER_OF_BINS];
  static double Count,Delta=0.0;
  double r2;
  VECTOR dr;
  FILE *FilePtr;

  switch(Ichoise)
  {
    case INITIALIZE:
      for(i=0;i<MAXIMUM_NUMBER_OF_BINS;i++)
        Distribution[i]=0.0;

      Count=0.0;
      // Delta  = Binsize
      Delta=0.5*BOXSIZE/(MAXIMUM_NUMBER_OF_BINS-1.0);
      break;
    case SAMPLE:
      // Sample The Radial Distribution Function
      // Loop Over All Particle Pairs
      // See Frenkel/Smit P. 77

      // Start Modification

      // End   Modification
      break;
    case WRITE_RESULTS:
      // write results to disk
      FilePtr=fopen("rdf.dat","w");
      for(i=0;i<MAXIMUM_NUMBER_OF_BINS-1;i++)
      {
        r2=M_PI*SQR(Delta)*NumberOfParticles*(NumberOfParticles-1)*(SQR(i+1)-(SQR(i)));
        fprintf(FilePtr,"%f %f\n",(i+0.5)*Delta,Distribution[i]*SQR(BOXSIZE)/(r2*Count));
      }
      fclose(FilePtr);
  }
}
