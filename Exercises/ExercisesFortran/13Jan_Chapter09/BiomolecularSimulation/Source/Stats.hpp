#ifndef _Stats_H_
#define _Stats_H_

#include "Residue2D.hpp"

//#include "EnergyMap.hpp"

#define NATIVE

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>



//#define NATIVE
//#define AGGR
//#define EXCL_THR_C

using namespace std;

class Lattice2D;
class Chain2D;


class Stats{
public:
  //inline StatsMethods.hpp
  Stats(); 
  static void setLoopUpCrossNbs();
  Stats operator+(const Stats & s) const;
  Stats & operator=(const Stats & s);
  Stats & operator+= (const Stats& s);
  Stats & operator-= (const Stats& s);
  const bool  operator!= (const Stats& s) const;
  int getCext(){return Cext;};
  int getEtot(){return Etot;};
  int getNtot(){return Ntot;};
  Stats  delta(const Stats & Snew, const Stats & Sold)const;
  int getDeltaE(const Stats & old)const;  

  void localStats(Residue2D * res,const Pos2D & pos,int resAA, Lattice2D * l);
  void localStats(Residue2D * res,const Pos2D & pos,Lattice2D * l);
  void localStatsExclude(Residue2D * res,const Pos2D & pos, Lattice2D * l,int start,int end);
  
  
  //Stats.cpp
  void getLatticeStats(Lattice2D * l);
  void get_Eint_Cint(Chain2D *  c,Lattice2D *l);
  void printCout();
  string print2string();
  

  //static void setLattice(Lattice * l);
  
private:
  int Etot;
  int Ctot;  
  int Cext;
  int Eint;
  int Cint;
  int Eext;
  //  friend class EnergyMap;
  
  
  //int Ctot;  
#ifdef NATIVE
  int Nint;
  int Next;
  int Ntot;
#endif
 
  //static Lattice * l;
};

/*** NON-MEMBER FUNCTIONS ***/





#endif
