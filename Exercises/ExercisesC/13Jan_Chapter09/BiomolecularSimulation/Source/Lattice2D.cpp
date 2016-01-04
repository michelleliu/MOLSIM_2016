#include "Lattice2D.hpp"
#include "Stats.hpp"
#include "StatsMethods2D.hpp"
//#include "EnergyMap.hpp"
#include "AA.hpp"
#include "Native2D.hpp"
//#include "MC.hpp"

#include <fstream>


Pos2D local2D[4]={Pos2D(1,0),
		  Pos2D(0,1),
		  Pos2D(-1,0),
		  Pos2D(0,-1)
	      
};

extern string fn_aa;

Lattice2D::Lattice2D(){
  //Stats::setLattice(this);
  //  energyMap = new EnergyMap(this);
  for (int x=0;x<LX;x++){
    for (int y=0;y<LY;y++){
      r[x][y]=NULL;
    }
  }
}



void Lattice2D::setAA(string fnaa){
  aaInt = new AA(fnaa);
}




void Lattice2D::readPDBMultiChain(string s){
  //TO DO check if correct format
  cout<<"reading file "<<s<<endl;
  ifstream infile(s.c_str());
  if (!infile.is_open()){
    cout << "can't open file "<< s <<endl;
    exit(0);
  }
  int resnum;
  int cn =0;
  char line[500];
  char type[5],restype[4], atom[4];;
  float xcoord,ycoord,zcoord;
  float occupancy,bfactor;
  // bool lig=false;
  chains[cn]=new Chain2D(this,cn);
  
  while(!infile.eof()){
    if(cn>= MAX_CHAINS){
      cout<<"Too many chains in "<< s <<endl;
      exit(0);
    }
    infile.getline(line,2000);
    type[0]=' ';
    sscanf(line,"%4s%*9c%2s%*2c%3s%*2c%4d%*4c%9f%9f%9f%4f%4f\n",type,atom,restype,&resnum,&xcoord,&ycoord,&zcoord,&occupancy,&bfactor);
    if ((strcmp(type,"ATOM")==0) &&(strcmp(line,"")!=0)){

      if(resnum>= MAX_RES){
	cout<<"Too many residues in chain "<<cn<< " in file "<< s <<endl;
	exit(0);
      }
      cout<<type<<"*"<<atom<<"*"<<restype<<"*"<<resnum<<"*"<<xcoord<<"*"<<ycoord<<"*"<<zcoord<<"*"<<occupancy<<"*"<<bfactor<<endl;

      if(strcmp(atom,"CA")==0){
	
	Residue2D * res = new Residue2D; 
	res->pos.x =(int) (xcoord/3.0);
	res->pos.y = (int) (ycoord/3.0);
	//res->pos.z =(int) (zcoord/3.0);
	res->aa = aaInt->stringtoAA(restype);
	res->n =resnum;
	res->chainNum = cn;
	chains[cn]->residues[resnum] = res;
	r[res->pos.x][res->pos.y]= res;
	
      }else{
	cout<<"atom type not recognised: "<<atom<<endl;
	exit(1);
      }

    }else if(strcmp(type,"TER")==0){
      if(resnum>= LX ||resnum>= LY ){
	cout<<"warning: chain "<<cn<< " in file "<< s <<" may be too long for box "<<endl;
      }
      chains[cn]->N = resnum+1;
      
     
      // chains[cn]->setRandomSpins();
      

      cn++;
      //cout <<"cn "<< cn<< endl;
      //cout<<line<<endl;
      chains[cn]=new Chain2D(this,cn);
    }else if(strcmp(type,"LIG")==0){
      chains[cn]->frozen = true;
      chains[cn]->locked = true;
    }else if(strcmp(type,"RIG")==0){
      chains[cn]->frozen = true;
    }
  }
  nChains = cn;
  applyPeriodicBoundaries();
#ifdef NATIVE
  native = new Native2D(this);
#endif
  stats.getLatticeStats(this);
  cout << "number of chains: "<< nChains<< endl;
}


void Lattice2D::applyPeriodicBoundaries(){
  for(int nc=0;nc<nChains;nc++){
    Chain2D * ch = chains[nc];
    for (int n=0;n<ch->N;n++){
      Residue2D * res =ch->residues[n];
      emptyLatticePos(res->pos);
      res->pos.periodicBoundary();  
      setResidue(res->pos,res);
    }
  }
}


void Lattice2D::writePDBMultiChain(string s){
  ofstream outs;
  cout<< "printing to "<< s<<endl;
  outs.open(s.c_str()); 
  outs << "MODEL "<<endl;
  //outs <<  setw(7) << moviestep<<endl;
  for(int nc=0;nc<nChains;nc++){
    Chain2D * ch = chains[nc];
    if(ch->frozen && ch->locked){
      outs << "LIG"<<endl ;
    }else if(ch->frozen){
      outs << "RIG"<<endl;
    }
    for (int n=0;n<ch->N;n++){
      Residue2D * res =ch->residues[n];
      outs << "ATOM  " ;
      outs << setw(5)<< n;
      outs << "  CA  ";
      outs << setw(3)<< aaInt->getAAfrom(res->aa);
      outs << " "<< (char)('A' +(res->chainNum)%('A' - 'z'));    // " A";
      outs << setw(4)<<n;
      outs << "   ";
      outs << setw(8)<< 3.0* (float) res->pos.x ;
      outs << setw(8)<< 3.0*(float) res->pos.y ;
      outs << setw(8)<< 0.0;  //3.0*(float) res->pos.z ;
     
      outs << "   1.00 22.00";
     
      outs << endl;
    }
    outs <<"TER   " << endl; 
  }
  outs.close();
}

void Lattice2D::deleteAllChains(){
  for(int nc=0;nc<nChains;nc++){
    for(int n=0;n<chains[nc]->N;n++){
      Residue2D * res=chains[nc]->residues[n];
      emptyLatticePos(res->pos);
    }
    delete chains[nc];
  }
  nChains=0;
}


int Lattice2D::deleteChain(int chain){
  //cout<<"deleting chain: "<< chain<<endl;
  // delete residues from lattice
  if (chains[chain]==NULL) return -1;
  Stats s;
  s.get_Eint_Cint(chains[chain],this);
  stats -= s;

  for(int n=0;n<chains[chain]->N;n++){
    Residue2D * res=chains[chain]->residues[n];
    r[res->pos.x][res->pos.y]=NULL;
  }
  delete chains[chain];
  nChains--;
  //freeChains--;

  // replace last chain, to new position, changeing chainNum in residues
  if(chain != nChains){
    for(int n=0;n<chains[nChains]->N;n++){
      chains[nChains]->residues[n]->chainNum=chain;
    }
    chains[chain] =chains[nChains];
    chains[chain]->setChainNum(chain);
  }
  //cout<<"finished deleting chain: "<< chain<<endl;
  return 0;
}


bool Lattice2D::checkStats(){
  Stats newStats;
  newStats.getLatticeStats(this);
  if(stats != newStats){
    cout<< "ERROR IN STATS"<<endl;
    cout<< "calculated stats:"<<endl;
    newStats.printCout();
    cout<< "additive stats:"<<endl;
    stats.printCout();
    writePDBMultiChain("errorInStats.pdb");
    exit(1);
    return false;
  }else{
    return true;
  }

}

void Lattice2D::setNative(){
  if(native !=NULL)delete native;
  native = new Native2D(this);
  native->setNative(this);
  stats.getLatticeStats(this);
}



int mainOLD(){

  Lattice2D* l=new Lattice2D();
  l->readPDBMultiChain("Seq2D.pdb");
  double pGlobal=0.1;
  for(int i=0;i<20000000;i++){
    if(drand48() >= pGlobal){
      if(!l->chains[0]->frozen)	l->chains[0]->localMove();
    }else{
     
    }
   
  }
  l->stats.printCout();
  Stats actualStats;actualStats.getLatticeStats(l);
  actualStats.printCout();
  l->writePDBMultiChain("out.pdb");
  //l->energyMap->printGnuPlotDataWR("stats",1,0);
  return(0);
}


