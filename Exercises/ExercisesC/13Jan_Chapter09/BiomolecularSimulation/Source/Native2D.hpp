#ifndef _Native_H_
#define _Native_H_


#include "Lattice2D.hpp"
#include "Chain2D.hpp"
#include "Stats.hpp"



#ifdef NATIVE
class Native2D{
public:
  Native2D(Lattice2D *l);
  ~Native2D();
  void setNative(Lattice2D *l);
  bool isNativeContact(int chainA, int chainB, int resA,  int resB);
  bool isBelowNativeEnergy(int Etot);
  bool hasNativeStructure(int Ntot);
  int getNativeEnergy();
  int getNumNativeContacts();
  void getNativeStructure ();
private:
  bool ****contacts;
  int totCnat;
  int nativeEnergy;
  int numNativeChains;
  int nativeChainLengths[MAX_CHAINS];
  Residue2D nativeRes[MAX_CHAINS][MAX_RES];
};



inline
bool Native2D::isNativeContact(int chainA, int chainB,int resA, int resB){
  //check if right chainLengths
  if(resA>= nativeChainLengths[chainA]) return false;
  if(resB>= nativeChainLengths[chainB]) return false;
  //
  return contacts[chainA][chainB][resA][resB];
}

inline
bool Native2D::isBelowNativeEnergy(int Etot){
  return (Etot < nativeEnergy);
}

inline
bool Native2D::hasNativeStructure(int Ntot){
  return (Ntot == totCnat);
}

inline
int Native2D::getNativeEnergy(){
  return (nativeEnergy);
}

inline
int Native2D::getNumNativeContacts(){
  return (totCnat);
}


#endif
#endif
