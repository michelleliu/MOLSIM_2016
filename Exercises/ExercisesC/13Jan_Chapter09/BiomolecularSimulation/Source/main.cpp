#include "Simulation2D.hpp"
#include "Design2D.hpp"

using namespace std;
#include <cstdlib>

int main(){
  bool design = true;

  // default input structure
  string fn_pdb_default = "../input/default.pdb";

  // file name of amino acid interaction matrix
  string fn_aa = "../input/matrix4.dat";

  // file name of amino acid distribution
  string fn_aa_distr = "../input/distribution4.dat";

  // ask user to design or simulate
  char answer;
  printf("Design a sequence (y/n)?");
  cin >> answer;
  cin.ignore();
  if(answer=='y'){
    design = true;
  }else if(answer=='n'){
    design = false;
  }else{
    printf("%c No valid option\n",answer);
    exit(1);
  }

  // ask user for a pdb structure
  string fn_pdb;
  cout<< "Give file name of pdb structure ";
  cout<< "(default "<<fn_pdb_default<<")"<<endl;
  getline(cin,fn_pdb);
  if(fn_pdb == ""){
    fn_pdb = fn_pdb_default;
  }
  cout << "using " <<fn_pdb <<endl;


  // start design
  if(design){
    printf("Starting design program\n");

    // START CODING
    // change your design temperature here
    double temperature, Dummy1, Dummy2, Dummy3;
    int designcase;
    FILE *FilePtr;
    FilePtr=fopen("input","r");
    fscanf(FilePtr,"%*s %*s %*s %*s\n");
    fscanf(FilePtr,"%lf %d %lf %lf",&temperature, &designcase, &Dummy2, &Dummy3);
    fclose(FilePtr);

    // END CODING

    int mc_steps =10000000;

    // initialise designProgram
    // fn_aa gives the pair potential
    Design2D * designProgram = new Design2D(fn_aa);



    // initialise PDB amino acid distribution for design
    designProgram->initStructure(fn_pdb);

    // design procedure using distribution of amino acids
    designProgram->initAADistribution(fn_aa_distr);

    // run design MC, see Design2D.cpp
    designProgram->designProcedure(temperature, mc_steps, designcase);

    // write design to design.pdb and statistics to stdout
    designProgram->writeDesignStats("../output/design.pdb");


    // start simulation
  }else{
    printf("Starting simulation \n");




    // Simulation
    Simulation2D * program = new Simulation2D(fn_aa);
    int mc_steps=1000000;
    int step=0;

    program->init(fn_pdb);
    // simulate the chain for several temperatures t
    int totalSamples;
    for(double t=0.1;t<0.4;t+=0.05){
      totalSamples = program->simulate(mc_steps,t,step);
      program->checkStats();
      step++;
    }
    int totalSteps = step;
    // print data (see Simulation2D.cpp)
    program->printCv("../output/Cv.txt", totalSteps,totalSamples);
    program->printCn("../output/Cn.txt",totalSteps);
    program->printPDB("../output/finalConfiguration.pdb");
    exit(1);
  }
}

