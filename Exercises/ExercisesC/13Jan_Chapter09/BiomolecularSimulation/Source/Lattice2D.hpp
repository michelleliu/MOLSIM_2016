#ifndef _Lattice_H_
#define _Lattice_H_

#define MAX_CHAINS 2


#include "Residue2D.hpp"
#include "Stats.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

using namespace std;

class Chain2D;
class Stats;
class Native2D;
//class EnergyMap;
class AA;
//struct Config;

class Lattice2D{
public:
  Lattice2D();

  Residue2D * getResidue(const Pos2D & p);
  void emptyLatticePos(const Pos2D & p);
  void setResidue(const Pos2D &  p,Residue2D * r);

  void readPDBMultiChain(string s);
  void writePDBMultiChain(string s);
  void printPeriodicPDB();
  //int insertChain(Config * config, int * sequence, double beta);
  int deleteChain(int chainNum);
  void deleteAllChains();
  void setAA();
  void setAA(string fnaa);
  bool checkStats();
  double getBetaMoves();
  double getRealBeta();
  void setBetaMoves(double b);
  void setNative(string fn_native);
  void setNative();
  void applyPeriodicBoundaries();
  Chain2D * chains[MAX_CHAINS];
  Residue2D * r[LX][LY];
  Stats stats;
  int nChains;
  Native2D * native;
  //EnergyMap * energyMap;
  AA * aaInt;
private:
  double betaMoves;

  //int freeChains;
 
 
};



inline
Residue2D * Lattice2D::getResidue(const Pos2D & p){
  return r[p.x][p.y];
}


inline
void Lattice2D::emptyLatticePos(const Pos2D & p){
  r[p.x][p.y]=NULL;
}

inline
void Lattice2D::setResidue(const Pos2D & p,Residue2D * res){
  r[p.x][p.y]=res;
}

inline
void Lattice2D::setBetaMoves(double m){
  betaMoves=m;
}


inline
double Lattice2D::getBetaMoves(){
  return betaMoves;
}

inline
double Lattice2D::getRealBeta(){
  return 100*betaMoves;
}


#endif
