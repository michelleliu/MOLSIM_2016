#ifndef _Chain_H_
#define _Chain_H_

#include "Residue2D.hpp"

#define MAX_RES 30


class Lattice2D;

class Chain2D{
public:
  Chain2D(Lattice2D *l,int chnum);
  ~Chain2D();
  int localMove();
  int globalMove();
  int rotatePoint();
  //void setBeta(double beta);
  //double getBeta();
  void setChainNum(int n);
 

  int getCext();
  bool hasCext();
  void setRandomSpins();
  void newRandomSpin(Residue2D * res);
  int shuffleSpin();
  void setStateResidues();
  int N;
  Residue2D * residues[MAX_RES];
  bool frozen;
  bool locked;
  
  //int Cext;
private:
  Lattice2D * l;
  int crankshaft(int n1, int n2);
 
  int translate();
  Pos2D rotationPosition(Pos2D posOld, int rotationDir,Pos2D posPivot);
  bool strandPossible(Residue2D * res);
  int chainNum;
  //  double betaMoves;
 
};


#endif
