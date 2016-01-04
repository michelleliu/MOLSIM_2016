#include "Pos2D.hpp"


#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>


using namespace std;


void Pos2D::printCout(){
  cout<<"x:"<<x<<" y:"<<y<<endl;
}


string Pos2D::toString(){
  std::stringstream myStream;
  myStream<<"x:"<<x<<" y:"<<y;
  return myStream.str();
}
