#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

// real Part; also direct calculation
// loop over all particle pairs
void Realspace(double *Ureal)
{
  int i,j;
  VECTOR dr;
  double U,r,r2;

  U=0.0;

  // start modification
  // 1. For all particle pairs calculate the distance in x, y, and z.
  for(i=0;i<NumberOfParticles-1;i++) {
    for(j=i+1;j<NumberOfParticles;j++) {
        // FIXME do something

    }
  }
  // 2. Apply periodic boundary conditions where necessary.
  // 3. Calculate the real-space contribution to the energy.
  // (Use the "ErrorFunctionComplement(x)" function.)


  // end  modification
  *Ureal=U;
}
