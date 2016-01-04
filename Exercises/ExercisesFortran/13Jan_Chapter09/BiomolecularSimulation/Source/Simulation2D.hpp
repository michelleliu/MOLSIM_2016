#ifndef __LATTICESIM_HPP_
#define __LATTICESIM_HPP_

#include <cstring>
#include <iostream>
#include <cstdio>

using namespace std;
#define SEQUENCE_LENGTH 30
#define NATIVE
#define MAX_TEMPS 100
#define MAX_SAMPLES 10000

class Lattice2D;

/// Contains the simulation framework
/// edit the Simulation2D.cpp file for sampling
class Simulation2D 
{
public: 
  Simulation2D(string fn_aa);
  void init(string fn_pdb);
  int  simulate ( int	nMcSteps, double temperature, int step );
  void printCv(string fn_out, int total_steps, int totalSamples);
  void printCn(string fn_out, int total_steps);
  void printPDB(string fn_pdb);


  double  getEnergy();
  double getNativeContacts();
  double getTotalNativeContacts();
  void setTemperature	(double dTemperature);
  void stretch();
  void fold();
  void checkStats();
protected:
  void setInitialConfig();
  
  /// lattice that holds the chain 
  Lattice2D * lattice;
private:

  /// keeps a list with temperatures
  double temperatureTable[MAX_TEMPS];

  /// keeps ensemble average of internal energy <E>
  double energyTable[MAX_TEMPS][MAX_SAMPLES];

  /// keeps ensemble avarage of nativeContacts <Cn>
  double nativeContactsTable[MAX_TEMPS];
};



#endif
