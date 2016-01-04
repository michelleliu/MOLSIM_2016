#include "LatticeSim.hpp"



//////// PROTECTED STUFF
//   ... your stuff ...


void printTBase(CAABase& pBases){
  int nBases=pBases.size();
  cout<<endl;
  for(int i=0;i<nBases;i++){
    cout<<"res: "<<i;
    cout<<" x: "<<pBases[i].x;
    cout<<" y: "<<pBases[i].y;
    cout<<" type: "<<pBases[i].t;
    cout<<endl;
  }
}

void printTContact(CAAContact& pContact){
  int nBases=pContact.size();
  for(int i=0;i<nBases;i++){
    cout<<"contact: "<< i;
    cout<<" res a: "<<pContact[i].idA;
    cout<<" res b: "<<pContact[i].idB;
    cout<<" E: "<<pContact[i].energy;
    cout<<endl;
  }
}



void printInfo(LatticeSim & program){
  cout<<"E: "<<program.getEnergy()<<endl;
  cout<<"maxE: :"<<program.getMaxEnergy()<<endl;
  cout<<"T: "<<program.getTemperature()<<endl;
  cout<<"N: "<<program.getDegreeOfFolding()<<endl; 
  cout<<"Var: "<<(double) program.getVariation()<<endl;
}

/*int main(){
  LatticeSim * program = new LatticeSim;
  program->init();
  

  CAABase aBases;
  aBases.resize(12);



  aBases[0].x = -3; aBases[0].y =   1; aBases[0].t =  0;
  aBases[1].x = -2; aBases[1].y =   1; aBases[1].t =  0;
  aBases[2].x = -1; aBases[2].y =   1; aBases[2].t =  0;
  aBases[3].x =  0; aBases[3].y =   1; aBases[3].t =  0;
  aBases[4].x =  1; aBases[4].y =   1; aBases[4].t =  0;
  aBases[5].x =  2; aBases[5].y =   1; aBases[5].t =  0;
  aBases[6].x =  2; aBases[6].y =   0; aBases[6].t =  0;
  aBases[7].x =  1; aBases[7].y =   0; aBases[7].t =  0;
  aBases[8].x =  0; aBases[8].y =   0; aBases[8].t =  0;
  aBases[9].x = -1; aBases[9].y =   0; aBases[9].t =  0;
  aBases[10].x =-2; aBases[10].y =  0; aBases[10].t = 0;
  aBases[11].x =-3; aBases[11].y =  0; aBases[11].t = 0;
  
  program->setNativeStructure(aBases);
  program->optimizeCode();
  program->getStructure(aBases);
  printTBase(aBases);

  }*/

int main(){
  LatticeSim * program = new LatticeSim;
  program->init();
  //int nB;
  CAABase  pB;
  CAAContact pC;
  program->getStructure(pB);
  program->setStructure(pB);
  printTBase(pB);
  printInfo(*program);
  program->checkStats();
  program->setTemperature(0.2);
  //program->setNativeStructure(pB);
  for(int i=0;i<10000;i++){
    //cout<<"round "<<i<<endl;
    program->simulate(100);
    program->getStructureAndContacts(pB,pC);
    
    program->setStructure(pB);
    // printTBase(pB);
    //printInfo(*program);
    //program->checkStats();
  }

  program->getStructure(pB);
  printTBase(pB);
  printInfo(*program);
  program->checkStats();

  program->optimizeCode();
  program->setNativeStructure(pB);
  program->stretch();
  program->checkStats();

  program->getStructure(pB);
  printTBase(pB);
  printInfo(*program);
  program->checkStats();
 
  for(int i=0;i<100;i++){
    //cout<<"round "<<i<<endl;
    program->simulate(1);
     program->getStructure(pB);
    
    // program->setStructure(pB);
     printTBase(pB);
    //printInfo(*program);
    //program->checkStats();
  }

  for(int i=0;i<10000;i++){
    program->simulate(10);
    program->stretch();
    program->simulate(10);
    program->fold();
  }
  program->getStructureAndContacts(pB,pC);
  printTBase(pB);
  printTContact(pC);
  printInfo(*program);
  program->checkStats();

 

  exit(1);
  }
