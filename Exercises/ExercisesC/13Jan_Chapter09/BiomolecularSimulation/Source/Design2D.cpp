/**
 * \file   Design2D.cpp
 * \date   Oct 2010
 * \author Sanne Abeln
 * \brief  Implements Design2D
*/



#include "Design2D.hpp"
#include "Lattice2D.hpp"
#include "Chain2D.hpp"
#include "Residue2D.hpp"
#include "AA.hpp"
#include "StatsMethods2D.hpp"

#include <fstream>

/// creates a new instance, while assigning a lattice
Design2D::Design2D(string fn_aa){
  mainLattice = new Lattice2D;
  mainLattice->setAA(fn_aa);
  q_variationStrength = 0.05;
}



/// initializes the instance, with a new sequence
/// \param fn_in should be the filename of a PDB style file
void Design2D::initStructure(string fn_in){
  mainLattice->readPDBMultiChain(fn_in);
  initCountAA();
}

/// initializes the instance, with a new sequence AA distribution
/// \param fn_in should be the filename of a PDB style file
void Design2D::initAADistribution(string fn_in){
  setAADistribution(fn_in,mainLattice->aaInt);
}


/// read an amino acid distributions
/// which is used by the design process to keep the variance similar
/// Note, the sequence file and AA interaction matrix
/// need to be read in first
/// \param fn_in file name containing aa distribution
/// \param aai reference to AA interaction matrix
void Design2D::setAADistribution(string fn_in, AA* aai){
  int numaa = aai->numberOfAAs;
  int sumFreq=0;
  int freqPDB[MAX_NUM_RES];

  //// read file
  ifstream aaFile (fn_in.c_str());
  if (!aaFile.is_open()){
    cout << "could not open file: "<<fn_in.c_str()<<  endl;
    exit(1);
  }
  int i=0;
  char line[500];
  while(aaFile.getline(line,200)){
    //// split line into fields
    char aa_c[5];
    int freq;
    sscanf(line,"%3s%d",aa_c,&freq);
    string aa = string(aa_c);

    cout<< aa<<endl;
    //// check if aa order is same as in neteraction matrix
    if(aa != "HOH"){
      if( aai->stringtoAA(aa)!=i){
        cout<< "amino acid order in file "<<fn_in<<" does not match interaction matrix";
        exit(1);
      }
      freqPDB[i]=freq;
      sumFreq+=freq;
    }
    i++;
  }

  for(int i=0;i<numaa;i++){
    expAA[i]=freqPDB[i]*seqLength/sumFreq;
  }
}

/// the procedure to design a new sequence
/// \param temperature sets the temperature for the energy
/// \param numSteps the number of MC steps
void Design2D::designProcedure(double temperature, long numSteps, int designcase){


  mainLattice->setBetaMoves(0.01/temperature);

  mainLattice->stats.getLatticeStats(mainLattice);
  systemVar = getSystemVar();
  printf("\ninitial systemVar: %f\n",systemVar);
  for(long i=0;i<numSteps;i++){
    // change AA
    changeAA((int) floor(mainLattice->nChains * drand48()),designcase);
    //checkVariance();
  }
  mainLattice->checkStats();
  checkVariance();
}



/// check the calculated variance, with the recorded variance
void Design2D::checkVariance(){
  double calcVar = getSystemVar();
  if(systemVar*1.01 < calcVar || systemVar*0.99 > calcVar){
    cout <<"WARNING! in design process accummulative variance:"<<systemVar;
    cout <<" does not match calc variance:" <<getSystemVar()<<endl;
    //exit(1);
 }
}


/// substitutes an amino acid, and updates the statistics
/// \param chain integer of the chain to design
void Design2D::changeAA(int chain, int designcase){
  // get pointer to chain (and sequence)
  Chain2D * c = mainLattice->chains[chain];

  // choose random residues
  int n = (int) floor(c->N * drand48());

  // get pointer to residue
  Residue2D * res = c->residues[n];

  // get old amino acid (integer)
  int oldAA = res->aa;

  // choose random new amino acid
  int newAA =  (int) floor(mainLattice->aaInt->getNumberOfAA() * drand48());

  // if nothing changed return
  if(oldAA==newAA)return;


  /* START CODING HERE */
  // calculate the incremtal variance
  // the array countAA holds the counts of each amino acid
  // the array expAA holds an observed distribution of AAs

  // localVar should be the change in variance
  // accVar should be your acceptance ratio


  double localVar=1.0;
  double accVar=1.0;

  if(designcase==1) {
  // case 1 using variance
  // try to avoid using the factorial function
  localVar*=(double)(countAA[oldAA]/(countAA[newAA]+1.0));
  accVar*=pow(localVar,1.0/q_variationStrength);
  }
  else if(designcase==2) {
  // case 2 using an amino acid distribution in expAA
  int numOfRes = mainLattice->aaInt->getNumberOfAA();
  localVar=0.0;
  localVar+=-pow(expAA[oldAA]-countAA[oldAA]/numOfRes,2.0);
  localVar+=pow(expAA[oldAA]-(countAA[oldAA]-1.0)/numOfRes,2.0);
  localVar+=pow(expAA[newAA]-(countAA[newAA]+1.0)/numOfRes,2.0);
  localVar+=-pow(expAA[newAA]-countAA[newAA]/numOfRes,2.0);
  accVar=exp(-1.0/(q_variationStrength*localVar));
  }

  /*  END CODING HERE */




  // accept move according to variance
  if(!(drand48()<accVar)){
    return;
  }

  // calculate new and old statistics
  // (don't worry about the details)
  Stats oldLocalStats;
  Stats newLocalStats;
  oldLocalStats.localStats(res,res->pos,oldAA,mainLattice);
  newLocalStats.localStats(res,res->pos, newAA,mainLattice);
  Stats newLatticeStats =  mainLattice->stats.delta(newLocalStats,oldLocalStats);
  int dE= newLocalStats.getDeltaE(oldLocalStats);

  // calcute boltzman factor
  double boltz = exp(((-(float) (dE)) )*  mainLattice->getBetaMoves());
  bool accept=  dE<=0 || (drand48() < boltz);

  // accept move according to energy
  if(accept){
    // update the amino acids in the sequence
    res->aa=newAA;
    // update energy
    mainLattice->stats = newLatticeStats;
    // update AAcount array
    countAA[oldAA]--;
    countAA[newAA]++;
    // update the incremental variance
    if(oldAA!=newAA){
      /* START CODING HERE */

      if(designcase==1) {
      // case 1
      systemVar *= localVar;
      }
      else if(designcase==2) {
      // case 2
      systemVar += localVar;
      }

      /*  END CODING HERE */

    }
  }

}


/// set initial counts of amino acids
void Design2D::initCountAA(){
  for(int i=0;i<MAX_NUM_RES;i++){
    countAA[i]=0;
  }
  seqLength=0;
  for(int nc=0;nc<mainLattice->nChains;nc++){
    Chain2D * ch = mainLattice->chains[nc];
    for (int n=0;n<ch->N;n++){
      Residue2D * res =ch->residues[n];
      countAA[res->aa]++;
      seqLength++;
    }
  }
  systemVar = getSystemVar();
}

/// calculates factorial
long factorial(long number) {
  long temp;
  if(number <= 1) return 1;
  temp = number * factorial(number - 1);
  return temp;
}

/// calculates the system variation
/// for case 1
double Design2D::getSystemVar(){
  double ni=1.0;
  int numOfRes = mainLattice->aaInt->getNumberOfAA();
  for(int aa=0;aa<numOfRes ;aa++){
    if(countAA[aa]){
      ni *= factorial(countAA[aa]);
    }
  }
  return (1.0)/((double) ni);
}

/// calculates the system variation
/// for case 2
double Design2D::getSystemVar2(){
  double dist=0;
  int numOfRes = mainLattice->aaInt->getNumberOfAA();
  for(int i=0;i<numOfRes;i++){
    if(countAA[i]){
      dist += (countAA[i] - expAA[i])*(countAA[i] - expAA[i]);
    }
  }
  return dist;
}


double Design2D::getMaxSequenceVariation(){
  // int numAA=mainLattice->aaInt->getNumberOfAA();
  double numAA= (double) mainLattice->aaInt->getNumberOfAA();
  return (1.0/ pow(factorial( seqLength/numAA),numAA));
}

double Design2D::getMinSequenceVariation(){
  //int numAA=mainLattice->aaInt->getNumberOfAA();
  return 1.0/(factorial( seqLength));
}

/// write statistics on the Design process to the standard output
void Design2D::writeDesignStats(string fn_pdb){
  mainLattice->writePDBMultiChain(fn_pdb);
  //q_variationStrength
  cout<<"variation strength (q):"<<q_variationStrength<<endl;
  //energyT
  cout<<"energyT:"<<mainLattice->getBetaMoves()<<endl;
  //sequence
  for(int nc=0;nc<mainLattice->nChains;nc++){
    Chain2D * ch = mainLattice->chains[nc];
    cout<<"chain no:"<< nc<<endl;
    for (int n=0;n<ch->N;n++){
      Residue2D * res =ch->residues[n];
      cout<<mainLattice->aaInt->getAAfrom(res->aa);
    }
    cout<<endl;
  }
  //Stats internal energy
  cout<< mainLattice->stats.print2string();
  //sequence variance
  cout<<"systemsVar:"<<systemVar<<endl;
  //AA count
  cout<<"AA count:"<< endl;
  int numOfAA = mainLattice->aaInt->getNumberOfAA();
  for(int aa=0;aa<numOfAA;aa++){
    cout <<mainLattice->aaInt->getAAfrom(aa)<<" "<<countAA[aa]<<", ";
  }
  cout<<endl;
}
