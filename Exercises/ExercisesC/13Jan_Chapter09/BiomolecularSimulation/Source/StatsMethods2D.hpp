#ifndef _Stats_Methods_H_
#define _Stats_Methods_H_

#include <cmath>
#include <cstdlib>

#include "Stats.hpp"
#include "Lattice2D.hpp"
#include "Pos2D.hpp"
#include "AA.hpp"
#include "Chain2D.hpp"
#include "Native2D.hpp"



/////////////////////////////////////////////////////////////////////////
inline
Stats::Stats(){
  Eext=0;
  Eint=0;
  Etot=0;
  Cint=0;
  Cext=0;
  Ctot=0;
#ifdef NATIVE
  Nint=0;
  Next=0;
  Ntot=0;
#endif
};



inline
Stats Stats::operator+ (const Stats& s) const{
  Stats result;
  result.Cint = (Cint + s.Cint);
  result.Cext = (Cext + s.Cext);
  result.Ctot = (Ctot + s.Ctot);
  result.Eint = (Eint + s.Eint);
  result.Eext = (Eext + s.Eext);
  result.Etot = (Etot + s.Etot);
#ifdef NATIVE
  result.Nint = (Nint + s.Nint);
  result.Next = (Next + s.Next);
  result.Ntot = (Ntot + s.Ntot);
#endif
  return result;
}


inline
Stats & Stats::operator+= (const Stats& s){
  Stats result;
  Cint += ( s.Cint);
  Cext += ( s.Cext);
  Ctot += ( s.Ctot);
  Eint += ( s.Eint);
  Eext += ( s.Eext);
  Etot += ( s.Etot);
#ifdef NATIVE
  Nint += ( s.Nint);
  Next += ( s.Next);
  Ntot += ( s.Ntot);
#endif
  return (*this);
}

inline
Stats & Stats::operator-= (const Stats& s){
  Stats result;
  Cint -= ( s.Cint);
  Cext -= ( s.Cext);
  Ctot -= ( s.Ctot);
  Eint -= ( s.Eint);
  Eext -= ( s.Eext);
  Etot -= ( s.Etot);
#ifdef NATIVE
  Nint -= ( s.Nint);
  Next -= ( s.Next);
  Ntot -= ( s.Ntot);
#endif
  return (*this);
}

inline
Stats & Stats::operator=(const Stats & s){
  Cint =  s.Cint;
  Cext =  s.Cext;
  Ctot =  s.Ctot;
  Eint =  s.Eint;
  Eext =  s.Eext;
  Etot =  s.Etot;
#ifdef NATIVE
  Nint =  s.Nint;
  Next =  s.Next;
  Ntot =  s.Ntot ;
#endif
  return (*this);
}

inline
const bool  Stats::operator != (const Stats& s) const{
  bool ans = (Cint !=  s.Cint ||
	      Cext !=  s.Cext ||
	      Ctot !=  s.Ctot ||
	      Eint !=  s.Eint ||
	      Eext !=  s.Eext ||
	      Etot !=  s.Etot );
#ifdef NATIVE
  ans = (ans ||
	 Nint !=  s.Nint ||
	 Next !=  s.Next ||
	 Ntot !=  s.Ntot );
#endif
  return ans;
}

inline
Stats  Stats::delta(const Stats & Snew, const  Stats & Sold)const{  
  Stats out; 
  out.Eext = Eext + Snew.Eext - Sold.Eext;
  out.Eint = Eint + Snew.Eint - Sold.Eint;
  out.Etot = Etot + Snew.Etot - Sold.Etot;

  out.Cext = Cext + Snew.Cext - Sold.Cext;
  out.Cint = Cint + Snew.Cint - Sold.Cint;
  out.Ctot = Ctot + Snew.Ctot - Sold.Ctot;
#ifdef NATIVE	 
  out.Next = Next + Snew.Next - Sold.Next;
  out.Nint = Nint + Snew.Nint - Sold.Nint;
  out.Ntot = Ntot + Snew.Ntot - Sold.Ntot;
#endif
  return out;
}


inline
void Stats::localStats(Residue2D * res,const Pos2D & pos, Lattice2D * l){
  localStats(res,pos,res->aa,l);
}


inline
void Stats::localStats(Residue2D * res,const Pos2D & pos,int resAA, Lattice2D * l){ 
  // for Int(AA) and h-bond energy
  for(int k=0;k<4;k++){
    Pos2D posNB =local2D[k] + pos;
    posNB.periodicBoundary();
    Residue2D * resNB = l->getResidue(posNB);
    if(resNB!=NULL){
      bool sameChain = resNB->chainNum==res->chainNum;
      if(abs(resNB->n - res->n)!=1 || ! sameChain){
	// Neighbouring sites on lattice, excluding neighbours in chain
	if(sameChain){
	  Cint++; 
	  Eint +=  l->aaInt->getInteraction(resAA, resNB->aa);
	}else{
	  Cext++;
	  Eext +=  l->aaInt->getInteraction(resAA, resNB->aa);
	}
#ifdef NATIVE
	if(l->native->isNativeContact(res->chainNum,resNB->chainNum,res->n,resNB->n)){
	  if(sameChain){Nint++;}else{Next++;}
	}
#endif
      }
    }
  }
  Etot=Eext+Eint;
  Ctot=Cint+Cext;
#ifdef NATIVE
  Ntot=Nint+Next;
#endif
}

inline
void Stats::localStatsExclude(Residue2D * res,const Pos2D & pos,
			      Lattice2D * l,int start,int end){
  
  // for Int(AA) and h-bond energy
  for(int k=0;k<4;k++){
    Pos2D posNB =local2D[k]+ pos;
    posNB.periodicBoundary();
    Residue2D * resNB = l->getResidue(posNB);
    if(resNB!=NULL){
      bool sameChain = resNB->chainNum==res->chainNum;
      if(!sameChain ||( resNB->n <start || resNB->n > end)){ // exluded residues
	if(abs(resNB->n - res->n)!=1 || ! sameChain){
	  // Neighbouring sites on lattice, excluding neighbours in chain
	  if(sameChain){
	    Cint++;
	    Eint +=  l->aaInt->getInteraction(res->aa, resNB->aa);
	  }else{
	    Cext++;
	    Eext +=  l->aaInt->getInteraction(res->aa, resNB->aa);
	  }
#ifdef NATIVE
	  if(l->native->isNativeContact(res->chainNum,resNB->chainNum,res->n,resNB->n)){
	    if(sameChain){Nint++;}else{Next++;}
	  }
#endif
	}
      }
    }
  }
  Etot=Eext+Eint;
  Ctot=Cint+Cext;
#ifdef NATIVE
  Ntot=Nint+Next;
#endif
}



inline
int Stats::getDeltaE(const Stats & old)const{
  return Etot -old.Etot;
}





#endif


