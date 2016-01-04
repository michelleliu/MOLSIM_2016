#include <string.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>


#include "LatticeSim.hpp"
#include "Lattice2D.hpp"
#include "Chain2D.hpp"
#include "Residue2D.hpp"
#include "Native2D.hpp"
#include "AA.hpp"
#include "Design2D.hpp"

string fn_aa = "input.dat";
string fn_pdb="default.pdb";

// function to create the simulation
ISimInterface* createLatticeSim()
{
  LatticeSim * program = new LatticeSim();
  program->init();
	return (ISimInterface*) program;
}

LatticeSim::LatticeSim(){
  lattice= new Lattice2D;
}

LatticeSim::~LatticeSim(){
  cout << "deletion not yet implemented - LatticeSim " <<endl;
}

void LatticeSim::init(){
  lattice->setAA(fn_aa);
  setInitialConfig();
  setTemperature(0.2);
  lattice->stats.getLatticeStats(lattice);
  lattice->setNative();
  lattice->stats.getLatticeStats(lattice);
}

void LatticeSim::free(){
  cout<<"free - not yet implemented - LatticeSim"<<endl;
}

void LatticeSim::setInitialConfig(){
  lattice->readPDBMultiChain(fn_pdb);
  lattice->setNative();
} 

void  LatticeSim::setNativeStructure (  CAABase&  rBases ){
  setStructure (rBases);
  lattice->setNative();
}


void LatticeSim::getNativeStructure  ( CAABase& rBases){
  lattice->native->getNativeStructure(rBases);
  getScreenCoords(rBases);
}

double LatticeSim::getMaxEnergy(){
  return (double) lattice->native->getNativeEnergy() *0.01;
}


void LatticeSim::setStructure ( CAABase & pBases ){
  
  lattice->deleteAllChains();
 
  lattice->chains[0]= new Chain2D(lattice,0);

  int nBases = pBases.size();
  for(int n=0;n<nBases;n++){
    //cout<<n<<endl;
    
    Residue2D * res = new Residue2D; 
    res->pos.x = pBases[n].x + LX/2;
    res->pos.y = pBases[n].y + LY/2;
    res->aa = pBases[n].t;
    res->chainNum =0;
    res->n =n;
   //insert in chain
    res->pos.periodicBoundary();
    lattice->chains[0]->residues[n]=res;
    //insert in lattice
    lattice->setResidue(res->pos,res);
  }
  lattice->chains[0]->N=nBases;

  lattice->nChains = 1;

  lattice->applyPeriodicBoundaries();
  
  lattice->stats.getLatticeStats(lattice);
 
}


void LatticeSim::getStructure ( CAABase&  pBases ){
  //TBase * pBases = new TBase[SEQUENCE_LENGTH];
  Chain2D * ch = lattice->chains[0];
  int nBases= ch->N;
  pBases.resize(nBases);
  //centre on screen
  for(int n=0;n<ch->N;n++){
    Residue2D * res = lattice->chains[0]->residues[n];  
    pBases[n].x= res->pos.x ;
    pBases[n].y=res->pos.y;
    pBases[n].t=res->aa;
    //cout<<pBases[n].t<<endl;
  }
  getScreenCoords(pBases);
  //printTBase(nBases,pBases);
  //pB = pBases;
  return;
}


void LatticeSim::getStructureAndContacts  ( CAABase & pB,
					    CAAContact & pContacts){
  
  getStructure(pB);
  int nB=pB.size();
  int numContacts=0;
  pContacts.resize(2*SEQUENCE_LENGTH);
  for (int i=0;i<nB;i++){
    for(int j=i+2;j<nB;j++){
      int diffx=pB[i].x - pB[j].x;
      int diffy=pB[i].y - pB[j].y;
      if((diffx==0 && (diffy==1 || diffy==-1))||
	 (diffy==0 && (diffx==1 || diffx==-1))){
	pContacts[numContacts].idA=i;
	pContacts[numContacts].idB=j;
	pContacts[numContacts].energy = 0.01 * 
	  lattice->aaInt->getInteraction(pB[i].t, pB[j].t);
	numContacts++;
	
      }	 
    }
  }
  pContacts.resize(numContacts);
}
 

double  LatticeSim::getEnergy(){
  return (double) lattice->stats.getEtot() *0.01;
}


double LatticeSim::getDegreeOfFolding(){
  double nativeC= (double) lattice->native->getNumNativeContacts();
  if(nativeC==0)return 0;
  return ((double) lattice->stats.getNtot())/(nativeC);
}

double LatticeSim::getVariation(){
  Design2D design(lattice);
  // double var= design.getSequenceVariation();
  
  double var = log(design.getSequenceVariation());
  double min= log(design.getMinSequenceVariation());
  double max = log(design.getMaxSequenceVariation());
  return (var -min)/(max-min);
}


double LatticeSim::getTemperature(){
  return(0.01/ lattice->getBetaMoves());
}


void LatticeSim::setTemperature	(double dTemperature){
  lattice->setBetaMoves(0.01/dTemperature);
}
 

void LatticeSim::optimizeCode(){
  double oldT=lattice->getBetaMoves();
  Design2D design(lattice);
  design.designProcedure(0.01, 1000);
  design.designProcedure(0.05, 100000);
  design.designProcedure(0.5, 100000);
  lattice->setBetaMoves(oldT);
  CAABase rBases;
  rBases.resize(30);
  getStructure(rBases);
  setNativeStructure(rBases);
}



void LatticeSim::simulate ( int	nMcSteps ){
  Chain2D * c= lattice->chains[0];
  double pGlobal= 0.1;
  for(int i=0;i<nMcSteps;i++){
    c->localMove();
    if(drand48() < pGlobal){
      c->globalMove();
    }
  }
}



void LatticeSim::stretch(){

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
}


void LatticeSim::fold(){
  CAABase  pB;
  getNativeStructure(pB);
  setNativeStructure(pB);
}

void LatticeSim::checkStats(){
  lattice->checkStats();
  lattice->stats.printCout();
}



/*********  PRIVATE **************/

void LatticeSim::getScreenCoords(CAABase& rBases){
  translate2screenCoords(rBases);
  rotate2screenCoords(rBases);
  mirror2screenCoords(rBases);
}


void LatticeSim::translate2screenCoords(CAABase& rBases){

  int nBases= rBases.size();
  //centre on screen
  int midx = rBases[nBases/2].x;
  int midy = rBases[nBases/2].y;
  
  int loBound = -SEQUENCE_LENGTH;
  int hiBound  = SEQUENCE_LENGTH; 
  for(int n=0;n<nBases;n++){  
    int diffx = rBases[n].x - midx;
    int diffy = rBases[n].y - midy;
    if( diffx < loBound){
      rBases[n].x= LX + diffx;
    }else if(diffx > hiBound){
      rBases[n].x= - LX + diffx;
    }else{
      rBases[n].x=diffx;
    }
    if( diffy < loBound){
      rBases[n].y= LY + diffy;
    }else if(diffy > hiBound){
      rBases[n].y= - LY + diffy;
    }else{
      rBases[n].y=diffy;
    } 
  }
}



void LatticeSim::rotate2screenCoords(CAABase& rBases){
  int nBases=rBases.size();
  int bm = nBases/2-1;
  int beforemidx = rBases[bm].x;
  int beforemidy = rBases[bm].y;
  if(beforemidx==1){
    return;
  }else if(beforemidx==-1){
    for(int n=0;n<nBases;n++){
      //int oldy=rBases[n].y;
      rBases[n].x= - rBases[n].x;
      rBases[n].y= - rBases[n].y;
    }
  }else if(beforemidy==1){
     for(int n=0;n<nBases;n++){
      int oldy=rBases[n].y;
      rBases[n].y= - rBases[n].x;
      rBases[n].x=  oldy;
    }
  }else if(beforemidy==-1){
    for(int n=0;n<nBases;n++){
      int oldy=rBases[n].y;
      rBases[n].y= rBases[n].x;
      rBases[n].x= - oldy;
    }
  }else{
    cout<< "haha dit mag helemaal niet\n";
  }
}


void LatticeSim::mirror2screenCoords(CAABase& rBases){
  int nBases=rBases.size();
  int mid = nBases/2;
  bool flip = false;
  bool turnfound =false;
  for(int n=mid; (n>=0) && (! turnfound);n--){
    if(rBases[n].y==-1){
      flip = true;
      turnfound=true;
    }else if(rBases[n].y==1){
      turnfound=true;
    }
  }
  for(int n=mid; (n<nBases) && (! turnfound);n++){
    if(rBases[n].y==1){
      flip = true;
      turnfound=true;
    }else if(rBases[n].y==-1){
      turnfound=true;
    }
  } 
  if(flip){
    for(int n=0;n<nBases;n++){
      rBases[n].y= - rBases[n].y;
    }
  }
}
