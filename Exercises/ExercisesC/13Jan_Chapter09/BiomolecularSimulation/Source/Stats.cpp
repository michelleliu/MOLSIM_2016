#include "Stats.hpp"
#include "StatsMethods2D.hpp"
#include "Lattice2D.hpp"
#include "Native2D.hpp"
#include "AA.hpp"


#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>


int getInteraction(int,int);


void Stats::getLatticeStats(Lattice2D * l){
  //cout<< "getting Lattice stats"<<endl;
  Eext=0;
  Eint=0;
  Cint=0;
  Cext=0;
  Ctot=0;
  Etot=0;
#ifdef NATIVE	 
  Nint=0;
  Next=0;
  Ntot=0;
#endif 
  for(int cn=0;cn< l->nChains ; cn++){
    Chain2D * c=l->chains[cn];
    for(int n=0;n< c->N;n++){
      // cout<<n<<endl;
      Residue2D * res= c->residues[n];
      Stats tmp;
      tmp.localStats(res,res->pos,l);
      (*this) += tmp;
    }
  }
  Eext /=2 ;
  Eint /=2;
  Cint /=2;
  Cext /=2;

  Etot=Eext+Eint;
  Ctot=Cint+Cext;

#ifdef NATIVE	 
  Nint /=2;
  Next /=2;
  Ntot=Nint+Next;
#endif 
}






void Stats::get_Eint_Cint(Chain2D*  c, Lattice2D *l){
  Eint =0;
  Cint =0;
  for(int n=0;n< c->N;n++){
      Residue2D * res= c->residues[n];
      Stats tmp;
      tmp.localStats(res,res->pos,l);
      (*this) += tmp;
  }

  Eint /=2;
  Cint /=2;
  Etot = Eint;
  Ctot = Cint;
  // set all others to zero
  Eext=0;
  Cext=0;
#ifdef NATIVE
  Nint/=2;  
  Next=0;
  Ntot=Nint;
#endif  
}



void Stats::printCout(){
  cout<<"Eint "<< Eint;
  cout<<" Eext "<< Eext;
  cout<<" Etot "<< Etot<<endl;

  cout<<"Cint "<<Cint;
  cout<<" Cext "<<Cext;
  cout<<" Ctot "<<Ctot<<endl;
  

#ifdef NATIVE	 
  cout<<"Nint "<<Nint;
  cout<<" Next "<<Next;
  cout<<" Ntot "<<Ntot<<endl;
#endif

}

string Stats::print2string(){
  std::stringstream myStream;
  myStream<<"Eint "<< Eint;
  myStream<<" Eext "<< Eext;
  myStream<<" Etot "<< Etot<<endl;

  myStream<<"Cint "<<Cint;
  myStream<<" Cext "<<Cext;
  myStream<<" Ctot "<<Ctot<<endl;
  
#ifdef NATIVE	 
  myStream<<"Nint "<<Nint;
  myStream<<" Next "<<Next;
  myStream<<" Ntot "<<Ntot<<endl;
#endif
  
  return myStream.str();
}




