#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

#define MAXIMUM_NUMBER_OF_BINS 1500

// samples the radial distribution function
void Sample(int Option)
{
  int i,j;
  static double Distribution[MAXIMUM_NUMBER_OF_BINS];
  static double Count,Delta=0.0;
  double r2;
  VECTOR dr;
  FILE *FilePtr;

  switch(Option)
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
      // See Frenkel/Smit p. 86

      // Start Modification
      Count += 1;
      for(i=0;i<NumberOfParticles-1;i++) {
        for(j=i+1;j<NumberOfParticles;j++) {
          dr.x=Positions[i].x-Positions[j].x;
          dr.y=Positions[i].y-Positions[j].y;

          // apply boundary conditions
          if(PBC)
          {
            dr.x-=BOXSIZE*rint(dr.x/BOXSIZE);
            dr.y-=BOXSIZE*rint(dr.y/BOXSIZE);
          }

          r2=SQR(dr.x)+SQR(dr.y);
          //    bin particles in Distribution[]
          if( sqrt(r2) < BOXSIZE/2.0 ) {
            Distribution[(int)( sqrt(r2)/Delta )] += 2.0;
          }
        }
      }

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
