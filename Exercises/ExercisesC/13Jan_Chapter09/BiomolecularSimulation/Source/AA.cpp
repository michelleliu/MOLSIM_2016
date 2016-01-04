#include "AA.hpp"

#include <cmath>
#include <fstream>
#include <cstdlib>


//#define DEBUG

//int averageInt[20];

//char  AA::int2aa[20][4] = {"Cys","Phe","Leu","Trp","Val","Ile","Met","His","Tyr","Ala","Gly","Pro","Asn","Thr","Ser","Arg","Gln","Asp","Lys","Glu"};

AA::AA(string filename){
  
 
  ifstream aaFile (filename.c_str());
  if (!aaFile.is_open()){
    cout << "could not open file: "<<filename.c_str()<<  endl;
    exit(1);
  }
  int i=0;
  while(!aaFile.eof()){
    char line[500];
    char aaName[4];
    float values[MAX_NUM_RES];
    aaFile.getline(line,2000);
    sscanf(line,"%3s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",aaName,values,values+1,values+2,values+3,values+4,values+5,values+6,values+7,values+8,values+9,values+10,values+11,values+12,values+13,values+14,values+15,values+16,values+17,values+18,values+19);
    // needs to read from diagonal and above 
    // cout<<aaName<<"*"<<endl;
    strncpy(int2aa[i],aaName,4);
    //  cout<<int2aa[i]<<"*"<<endl;
    for(int j=i;j<MAX_NUM_RES;j++){
      aaInt[i][j]=(int) rint(100*values[j]);
      aaInt[j][i]=(int) rint(100*values[j]);      
    }
    i++;
  }
  numberOfAAs=i-1;
#ifdef DEBUG
  for(int i=0;i<MAX_NUM_RES;i++){
    cout << int2aa[i] <<"\t";
    for(int j=0;j<MAX_NUM_RES;j++){
      cout<< aaInt[i][j]<<"\t";
    }
    cout<<endl;
  }
#endif
  
  cout<<numberOfAAs<<" different residues"<<endl; 
 // OPT would need reordering with other matrix input
  aaFile.close();
}


string AA::getAAfrom(int aa){

  return string(int2aa[aa]);
}

int AA::stringtoAA(string s){
  for(int i=0;i<numberOfAAs;i++){
    if(!s.compare(int2aa[i])){
      return i;
    }
  }
  cout << "aa not found:" <<s<<endl;
  exit(1);
  return(0);
}
