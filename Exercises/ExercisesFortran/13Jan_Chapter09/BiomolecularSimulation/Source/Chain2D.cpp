#include "Chain2D.hpp"
#include "Lattice2D.hpp"
#include "Stats.hpp"
#include "StatsMethods2D.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

//#define DEBUG

/******************Standard Move **************************/
/* Generate new configuration                        */
/* Check if newConfig illegal,if illegal return -1   */
/* Calculate energy difference between states        */
/* Sample boths states for statistics                */
/* accept or reject move                             */
/* reject: return 0                                  */
/* accept: update lattice and chain,return 1         */
/**********************************************************/



Chain2D::Chain2D(Lattice2D *l1,int chnum){
  l=l1;
  chainNum=chnum;
  frozen =false;
  locked=false;
}


Chain2D::~Chain2D(){
  for(int i=0;i<N;i++){
    delete residues[i];
  }
}

void Chain2D::setChainNum(int n){
  chainNum=n;
}


int Chain2D::localMove(){
  Pos2D posNew;
  int px,py; 
  //choose Residue to move
  int n = (int) floor(N * drand48());
  Residue2D * res=residues[n];

  if(n<0 ||n >= N){ cout<<"residue not in range"<<n<<endl; exit(-1);} 
  
  if(n==0 || n == N-1){ // try end move
    int nn =1;
    if(n==(N-1)){nn = N-2;} // nearest neighbour   
    int rr= (int)floor(4 * drand48()); 
    // OPT slow, Behnaz: put links instead
    posNew = local2D[rr] + residues[nn]->pos;
    posNew.periodicBoundary();
  }else{ // try cornerflip
    px = (res->pos.x == residues[n-1]->pos.x) & (res->pos.x == residues[n+1]->pos.x);
    py = (res->pos.y == residues[n-1]->pos.y) & (res->pos.y == residues[n+1]->pos.y);
    if(px+py==1){
      // on a straight line
      return -1;
    }else{ 
      posNew.x =  residues[n-1]->pos.x + residues[n+1]->pos.x - res->pos.x;
      posNew.y =  residues[n-1]->pos.y + residues[n+1]->pos.y - res->pos.y;
    }
  }
  //check if no steric clash main chain + check if cranshaft possible
  if(l->getResidue(posNew)!=NULL){
    if(l->getResidue(posNew)->n == n-2 && l->getResidue(posNew)->chainNum == chainNum){
      return crankshaft(n-1,n);
    }else if(l->getResidue(posNew)->n == n+2 && l->getResidue(posNew)->chainNum== chainNum){
      return crankshaft(n,n+1);
    }else{
      return(-1);
    }
  }


  Stats oldLocalStats; 
  Stats newLocalStats; 
  oldLocalStats.localStats(res,res->pos,l);
  newLocalStats.localStats(res,posNew,l);

  Stats newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);
  
  int dE= newLocalStats.getDeltaE(oldLocalStats);
  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  //double bb = boltz/(boltz+1.0);
  /*if(l->energyMap != NULL){
    l->energyMap->mapStats(&l->stats,1.0 - bb);
    l->energyMap->mapStats(&newLatticeStats,bb);
    }*/
  bool accept=false;
  if(dE<=0 || drand48() < boltz){
    accept=true;
  }
  if(accept){
    // update position and lattice 
    l->emptyLatticePos(res->pos);
    res->pos = posNew;
    l->setResidue(posNew,res);
    // update lattice stats
    l->stats = newLatticeStats;


    if(dE==0){
      return(1);
    }else{
      return(2);
    }
  }else{
    //delete newLatticeStats;
    //delete newRes;
    return(0);
  }
}



int Chain2D::crankshaft(int n1, int n2){
  //cout << "cranshaft "<<n1<<" "<<n2<<endl;
  int n0,n3;
  Residue2D * resn1 = residues[n1];
  Residue2D * resn2 = residues[n2];
  
  // TO DO rewrite this:
  if(n1<n2){
    n0=n1-1;
    n3=n1+2;
  }else{
    n0=n1+1;
    n3=n1-2;
    //swap all
    int tmp =n0;
    n0=n3;
    n3=tmp;
    tmp=n1;
    n1=n2;
    n2=tmp;
  }
  Residue2D * resn0 = residues[n0];
  Residue2D * resn3 = residues[n3];
 
  if(n0<0 || n3<0 || n0>=N || n3>=N){
    return -3;
  }
  // OPT only need dirOrthogonal
  int dirStalk,dirLegs;
  
  if((resn1->pos.x - resn2->pos.x)==0){
    dirStalk=XDIR;
    dirLegs=YDIR;
  }else{
    dirStalk=XDIR;
    dirLegs=YDIR;
  }


  Pos2D old0 =Pos2D(resn0->pos.x,resn0->pos.y);
  Pos2D old3 =Pos2D(resn3->pos.x,resn3->pos.y);

  Pos2D new1;
  Pos2D new2;

  new1[dirStalk]= old0[dirStalk] +  old0[dirStalk] - resn1->pos[dirStalk] ;
  new2[dirStalk]= old3[dirStalk] +  old3[dirStalk] - resn2->pos[dirStalk];
 
  new1[dirLegs]= old0[dirLegs];
  new2[dirLegs]= old3[dirLegs];

  new1.periodicBoundary();
  new2.periodicBoundary();

  // test for collision
  if(l->getResidue(new1)!=NULL || 
     l->getResidue(new2)!=NULL){
    return -1;
  }

  Stats oldS1,oldS2,newS1,newS2;
  oldS1.localStatsExclude(resn1,resn1->pos,l,n2,n2);
  oldS2.localStatsExclude(resn2,resn2->pos,l,n1,n1);
  newS1.localStatsExclude(resn1,new1,l,n2,n2);
  newS2.localStatsExclude(resn2,new2,l,n1,n1);


  Stats oldLocalStats = oldS1 +oldS2;
  Stats newLocalStats = newS1 +newS2;

  Stats newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);
  int dE= newLocalStats.getDeltaE(oldLocalStats);

  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  //double bb = boltz/(boltz+1.0);
  /*if(l->energyMap != NULL){
    l->energyMap->mapStats(&l->stats,1.0 - bb);
    l->energyMap->mapStats(&newLatticeStats, bb);
    }*/
  bool accept=false;
  if(dE<=0 || drand48() < boltz){
    accept=true;
  }
  if(accept){
    
    // update position and lattice 
    l->emptyLatticePos(resn1->pos);
    l->emptyLatticePos(resn2->pos);
    resn1->pos = new1;
    resn2->pos = new2;
    l->setResidue(new1,resn1);
    l->setResidue(new2,resn2);
    // update lattice stats
    l->stats = newLatticeStats;
    
    if(dE==0){
      return(1);
    }else{
      return(2);
    }
  }else{
    
    return(0);
  }


}



int Chain2D::globalMove(){
  return rotatePoint();
  /*if(drand48() < 0.5){
    return translate();
  }else{
    return rotatePoint();
    }*/
}


int Chain2D::translate(){
  int dir= (int)floor(4 * drand48());

  Pos2D posNew[MAX_RES];

  for(int n=0; n<N;n++){
    Residue2D * res = residues[n];
    posNew[n] = res->pos+ local2D[dir];
    posNew[n].periodicBoundary();
    
    Residue2D * testRes= l->getResidue(posNew[n]);
    //check for collision
    if(testRes!=NULL && testRes->chainNum != res->chainNum ){
#ifdef DEBUG
      cout << "clash with ligand "<<endl;;
#endif
      return -1;
    }  
  }

  Stats oldLocalStats,newLocalStats;
  for(int n=0; n<N;n++){
    Stats oldS, newS; 
    Residue2D * res=residues[n];
    oldS.localStatsExclude(res,res->pos,l,0,N-1);
    newS.localStatsExclude(res,posNew[n],l,0,N-1);
    oldLocalStats += oldS;
    newLocalStats += newS;
  }

  Stats newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);
  int dE= newLocalStats.getDeltaE(oldLocalStats);
  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  //double bb = boltz/(boltz+1.0);
  /*if(l->energyMap != NULL){
    l->energyMap->mapStats(&l->stats,1.0 - bb);
    l->energyMap->mapStats(&newLatticeStats, bb);
    }*/
  bool accept=false;
  if(dE<=0 || drand48() < boltz){
    accept=true;
  }
  
  if(accept){
    //Cext+= newLocalStats.Cext - oldLocalStats.Cext;
    for(int n=0;n<N;n++){
      l->emptyLatticePos(residues[n]->pos);
    }
    for(int n=0;n<N;n++){
      residues[n]->pos=posNew[n];
      l->setResidue(posNew[n],residues[n]);
    }
    l->stats = newLatticeStats;
    if(dE==0){
      return(1);
    }else{
      return(2);
    }
  }else{
    return(0);
  }
}


int Chain2D::rotatePoint(){
  //cout<<"rotate point"<<endl;
  int rotationalDir= (int)floor(2 * drand48());
  int pivot= (int)floor(N * drand48());
  int startRes=0;
  int endRes=N-1;
  // OPT don't need to include pivot
  if(!frozen) {// if frozen only allow global moves
    if(drand48()<0.5){
      startRes= pivot;
    }else{
      endRes=pivot;
    }
  }
  //cout<<"pivot"<<pivot<<endl;
  //cout<<"startRes"<<startRes<<endl;
  //cout<<"endRes"<<endRes<<endl;
  Pos2D posNew[MAX_RES];
  // test for collisions
  for(int n=startRes; n<=endRes;n++){
    Residue2D * res = residues[n];
    posNew[n]= rotationPosition(res->pos,rotationalDir,residues[pivot]->pos);
#ifdef DEBUG
    cout<<"move "<<n<<" from " <<res->pos.toString()<<" to (" << posNew[n].toString() <<")"<<endl;
#endif
    Residue2D * testRes= l->getResidue(posNew[n]);
    if(testRes!=NULL && (testRes->chainNum != res->chainNum || testRes->n<startRes || testRes->n> endRes  )){
      return -1;
    }  
  }
  
  Stats oldLocalStats, newLocalStats;
  int exclStart =startRes;
  int exclEnd = endRes;
  /*if(startRes==pivot){
    exclStart = startRes +1;
  }else{
    exclEnd = endRes -1;
    }*/
#ifdef DEBUG
  cout<< "exclSt "<<exclStart<<" exclEnd "<<exclEnd<<endl;
#endif

  for(int n=exclStart; n<=exclEnd;n++){
    Residue2D * res=residues[n];
    Stats oldS, newS; 
    oldS.localStatsExclude(res,res->pos,l,exclStart,exclEnd);
    newS.localStatsExclude(res,posNew[n],l,exclStart,exclEnd);
    oldLocalStats += oldS;
    newLocalStats += newS;
  }
  
  Stats newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);

  int dE= newLocalStats.getDeltaE(oldLocalStats);

  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  // double bb = boltz/(boltz+1.0);
  /*if(l->energyMap != NULL){
    l->energyMap->mapStats(&l->stats,1.0 - bb);
    l->energyMap->mapStats(&newLatticeStats, bb);
    }*/
  bool accept=false;
  if(dE<=0 || drand48() < boltz){
    accept=true;
  }
  
  if(accept){
     //Cext+= newLocalStats.Cext - oldLocalStats.Cext;
    //cout <<"ACCEPT"<<endl;
    for(int n=startRes;n<=endRes;n++){
      l->emptyLatticePos(residues[n]->pos);
    }
    for(int n=startRes;n<=endRes;n++){
      residues[n]->pos=posNew[n];
      l->setResidue(posNew[n],residues[n]);
    }  
    l->stats = newLatticeStats;
   // set spins
  
   
    
    if(dE==0){
      return(1);
    }else{
      return(2);
    }
  }else{
    return(0);
  }
}


Pos2D Chain2D::rotationPosition(Pos2D posOldi, int rotationDir,Pos2D posPivot){
  Pos2D output;
  // int prim = rotationDir % 3;
  int sign =2*rotationDir - 1;
  Pos2D posOld =posOldi - posPivot ;
  
  output.x = sign*posOld.y + posPivot.x;
  output.y = -1*sign*posOld.x + posPivot.y ;
  
  output.periodicBoundary();
  return output  ;
}




