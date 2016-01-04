#ifndef __LATTICESIM_HPP_
#define __LATTICESIM_HPP_

#include <string.h>
#include <iostream>
#include <stdio.h>


using namespace std;




#define SEQUENCE_LENGTH 30
#define NATIVE



class Lattice2D;

class LatticeSim : public ISimInterface
{
  public: 
  //... ISimInterface methods ...
  LatticeSim();
  ~LatticeSim();

  void init();
  void free();
  void  setNativeStructure ( CAABase& rBases );
  void  getNativeStructure ( CAABase& rBases );
  double  getMaxEnergy();
  void  setStructure ( CAABase& rBases );
  void	getStructure ( CAABase& rBases );
  void	getStructureAndContacts ( CAABase& rBases, CAAContact &  rContacts);
  double  getEnergy();
  double getDegreeOfFolding();
  double getVariation();
  double getTemperature();
  void setTemperature	( double dTemperature );
  void optimizeCode();
  void	simulate ( int	nMcSteps );
  void stretch();
  void fold();
  void checkStats();
  protected:
  void setInitialConfig();
  void getScreenCoords(CAABase& rBases);
  void translate2screenCoords(CAABase& rBases);
  void rotate2screenCoords(CAABase& rBases);
  void mirror2screenCoords(CAABase& rBases);
  //   ... your stuff ...
 
  Lattice2D * lattice;

};




#endif
