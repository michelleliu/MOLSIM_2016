#include "Native2D.hpp"

#include <cmath>
#include <cstdlib>


#ifdef NATIVE


Native2D::Native2D(Lattice2D *l){
  numNativeChains=l->nChains;
  //cout<<"creating new native"<<endl;
  contacts = new  bool  *** [numNativeChains];
  for(int i=0;i<l->nChains;i++){
    contacts[i] = new bool ** [numNativeChains];
    int length_i =l->chains[i]->N;
    nativeChainLengths[i]=length_i;
    for(int j=0;j<l->nChains;j++){
      int length_j =l->chains[j]->N;
      contacts[i][j] = new  bool* [length_i];
      for(int k=0;k<length_i;k++){
	contacts[i][j][k]=new bool[length_j];
	for(int l=0;l<length_j;l++){
	  contacts[i][j][k][l]=false;
	}
      }
    }
  }
}


Native2D::~Native2D(){
  //cout<<"deleting native"<<endl;
  for(int i=0;i<numNativeChains;i++){
    for(int j=0;j<numNativeChains;j++){
      for(int k=0;k< nativeChainLengths[i];k++){
	delete[] contacts[i][j][k];
      }
      delete[] contacts[i][j];
    }
    delete[] contacts[i];
  }
  delete[] contacts;
  //cout<<"finished"<<endl;
}



void Native2D::setNative(Lattice2D * l){
  // cout<<"setting native"<<endl;
  //TO DO should check if not out of bounds
  int totalC=0;
  for(int nc=0;nc<l->nChains;nc++){
    for(int n=0;n<l->chains[nc]->N;n++){
      Residue2D * res = l->chains[nc]->residues[n];
      for(int k=0;k<4;k++){
	Pos2D posNB =local2D[k]+ res->pos;
	posNB.periodicBoundary();
	Residue2D * resNB = l->getResidue(posNB);
	if(resNB!=NULL){
	  if(resNB->chainNum !=  res->chainNum || abs(resNB->n - res->n) !=1){
	    contacts[res->chainNum][resNB->chainNum][res->n][resNB->n]=true;
	    contacts[resNB->chainNum][res->chainNum][resNB->n][res->n]=true;
	    totalC++;
	  }
	}
      }
    }
  }

  totalC =totalC/ 2;
  totCnat=totalC;
  
  for(int nc=0;nc<l->nChains;nc++){
    for(int i=0;i<l->chains[nc]->N;i++){
      Residue2D * res = l->chains[nc]->residues[i];
      nativeRes[nc][i].aa=res->aa;
      nativeRes[nc][i].pos.x =res->pos.x;
      nativeRes[nc][i].pos.y =res->pos.y;
      nativeRes[nc][i].chainNum=res->chainNum;
      nativeRes[nc][i].n=res->n;
    }
  }
  //cout<<"TOTAL NATIVE CONTACTS:" <<totCnat<<endl;
  //l->stats.getLatticeStats(l); //should be set after reading in
  nativeEnergy = l->stats.getEtot();
  //cout<<"NATIVE INTERNAL ENERGY:" <<nativeEnergy<<endl;
}
 
void Native2D::getNativeStructure  ( ){
  // only get first chain in box
  /*int nBases=nativeChainLengths[0];
  pBases.resize(nBases);
  for(int n=0;n<nativeChainLengths[0];n++){
    pBases[n].x=nativeRes[0][n].pos.x;
    pBases[n].y=nativeRes[0][n].pos.y;
    pBases[n].t=nativeRes[0][n].aa;
    }*/
}
#endif
