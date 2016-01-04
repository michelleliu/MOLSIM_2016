#ifndef _AA_H_
#define _AA_H_

#define MAX_NUM_RES 20


#include <cstdio>
#include <cstring>
#include <iostream>
using namespace std;

class AA{
public:
  AA(string filename);
  int getInteraction(int a1,int a2);
  int stringtoAA(string s);
  string getAAfrom(int aa);
  int getNumberOfAA();
  int numberOfAAs;
private:
  char int2aa[MAX_NUM_RES][4]; 
  int aaInt[MAX_NUM_RES][MAX_NUM_RES];
  
};


inline
int AA::getInteraction(int a1,int a2){
  return aaInt[a1][a2];
}

inline
int AA::getNumberOfAA(){
  return numberOfAAs;
}

#endif
