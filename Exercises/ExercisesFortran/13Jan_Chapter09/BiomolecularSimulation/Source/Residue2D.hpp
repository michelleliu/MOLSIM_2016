#ifndef _Residue_H_
#define _Residue_H_

#include "Pos2D.hpp"


#define XDIR 0
#define YDIR 1



extern Pos2D local2D[4];


struct Residue2D{
  int aa;
  Pos2D pos;
  int n; // number in chain
  int chainNum;
};


#endif
