/**
 * \file   Design2D.hpp
 * \date   Oct 2010
 * \author Sanne Abeln
 * \brief  Defines Design2D 
*/


#ifndef _DESIGN_H_
#define _DESIGN_H_

#include <cstring>
#include <cstdio>
#include <iomanip>

#include "AA.hpp"
using namespace std;

class Lattice2D;



/// class  defines functions and a framework to design a new sequence
class Design2D{
public:
  
  Design2D(string fn_aa);
  void initStructure(string fn_in);
  void initAADistribution(string fn_in);
  void designProcedure(double beta, long numSteps);
  void writeDesignStats(string fn_pdb);
  void setAADistribution(string fn_in, AA* aai);

  // DEPRECATED
  void initOldFormat(string fn_in);
  // DEPRECATED
  double getMaxSequenceVariation();
  // DEPRECATED
  double getMinSequenceVariation();

  /// design temperature
  double q_variationStrength;

private:
  double getSystemVar();
  double getSystemVar2();
  void changeAA(int chain);
  void initCountAA();
  void checkVariance();

  /// pointer to the lattice
  Lattice2D * mainLattice;

  /// keeps the number instances of each AA occuring in the sequence 
  int countAA[MAX_NUM_RES];

  /// keeps an observed amino acid distribution from pdb
  double expAA[MAX_NUM_RES];

  /// sequence length
  int seqLength;

  /// keeps the incremented sequence variation, updated by changeAA
  double systemVar;
};

#endif
