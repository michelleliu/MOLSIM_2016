#include <cstring>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <fstream>


#include "Simulation2D.hpp"
#include "Lattice2D.hpp"
#include "Chain2D.hpp"
#include "Residue2D.hpp"
#include "Native2D.hpp"
#include "AA.hpp"
#include "Design2D.hpp"


/// create new instance
Simulation2D::Simulation2D(string fn_aa){
  lattice= new Lattice2D;
  lattice->setAA(fn_aa);
}


/// initialise simulation and create lattice
void Simulation2D::init(string fn_pdb){
  
  
  lattice->readPDBMultiChain(fn_pdb);
  setTemperature(0.2);
  lattice->stats.getLatticeStats(lattice);
  lattice->setNative();
  lattice->stats.getLatticeStats(lattice);
}




/// starts a simulation round
/// \param nMcSteps the number of Monte Carlo steps
/// \param temperature the temperature at which the simulation is performed
/// \param step iteration number for statistics
/// returns the number of energy samples taken
int Simulation2D::simulate ( int nMcSteps, double temperature, int step){

  // get a pointer to the chain to simulate (first chain on the lattice)
  Chain2D * chain= lattice->chains[0];
 
 // set the temperature
  setTemperature(temperature);

  // set the fraction of global moves (point rotations)
  double pGlobal= 0.1;

  // keep the accumalitive statistics for the 
  // and number of native contacts (sumCn) 
  double sumCn=0;

  // calculate the the fraction of sampling
  int modSampling = nMcSteps / MAX_SAMPLES;
  int sample=0;

  // perform monte carlo moves
  for(int i=0;i<nMcSteps;i++){
    chain->localMove();
    if(drand48() < pGlobal){
      chain->globalMove();
    }
    //sample
    // update the number of native contacts
    sumCn +=  getNativeContacts();

    // store energy as often as 'modSampling' allows
    if(((i % modSampling) == 0) && sample < MAX_SAMPLES){
          energyTable[step][sample] = getEnergy();
      sample++;
    }
  }
  // store statistics
  temperatureTable[step]= temperature;
  // get fraction of native contacts
  nativeContactsTable[step]= sumCn/(nMcSteps * getTotalNativeContacts());

  // return the number of samples taken
  return sample;
}

/// prints the heat capacity statistics for all temperatures
/// \param fn_out the filename for output
void Simulation2D::printCv(string fn_out, int total_steps, int totalSamples){
  // open file
  ofstream outs;
  cout<< "printing to "<< fn_out<<endl;
  outs.open(fn_out.c_str()); 

  // loop over the temperatures
  for(int step =0; step< total_steps;step++){

    // START CODING HERE                                                       
    // calculate the heat capacity (Cv*kB)                                     
    // by using energyTable and temperatureTable  
    double Cv=1.0;

    //first calculate the ensemble average over the energy
    //double sumEnergy =0;
    
    // loop over all samples stored in energyTable, upto totalSamples
    //for (int sample=0;sample < totalSamples;sample++)...
      

    // then calculate the average energy fluctuation 
    //double sumSqDiff=0;
    //for (int sample=0;sample < totalSamples;sample++)...
   
    // finally calculate the Cv using temperatureTable
    double temperature = temperatureTable[step];
    // Cv =  ...

    // END CODING HERE

    // print the temperature and heat capacity
    outs <<temperature<<" ";
    outs <<Cv<<endl;
  }
  // close file
  outs.close();
}

/// prints the native contact statistics for all temperatures
/// \param fn_out the filename for output
void Simulation2D::printCn(string fn_out, int total_steps){

  // open file
  ofstream outs;
  cout<< "printing to "<< fn_out<<endl;
  outs.open(fn_out.c_str()); 

 // loop over the temperatures
  for(int step =0; step< total_steps;step++){
    // print the temperature and the ensemble average over
    // the number of native contacts
    outs << temperatureTable[step] <<" ";
    outs << nativeContactsTable[step] <<endl;
  }
  // close file
  outs.close();

}

void Simulation2D::printPDB(string fn_pdb){
  lattice->writePDBMultiChain(fn_pdb);
}




double  Simulation2D::getEnergy(){
  return (double) lattice->stats.getEtot() *0.01;
}


double Simulation2D::getNativeContacts(){
  return (double) lattice->stats.getNtot();
}

double Simulation2D::getTotalNativeContacts(){
  return (double) lattice->native->getNumNativeContacts();
}


/*double Simulation2D::getTemperature(){
  return(0.01/ lattice->getBetaMoves());
  }*/


void Simulation2D::setTemperature	(double dTemperature){
  lattice->setBetaMoves(0.01/dTemperature);
}
 

/*void Simulation2D::stretch(){

  Chain2D * ch = lattice->chains[0];
  int nBases= ch->N;

  for(int n=0;n<nBases;n++){
    //cout<<n<<endl;
    Residue2D * res = ch->residues[n];
    lattice->emptyLatticePos(res->pos);
    res->pos.x =  (LX/2) + n - (nBases/2);
    res->pos.y = LY/2;
    lattice->setResidue(res->pos,res);
   }

  lattice->applyPeriodicBoundaries();
  lattice->stats.getLatticeStats(lattice);
  //lattice->writePDBMultiChain("stretch.pdb");
  }*/


/*void Simulation2D::fold(){
  CAABase  pB;
  getNativeStructure(pB);
  setNativeStructure(pB);
  }*/

void Simulation2D::checkStats(){
  lattice->checkStats();
  lattice->stats.printCout();
}


