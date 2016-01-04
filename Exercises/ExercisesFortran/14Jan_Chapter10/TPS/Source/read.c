#include "tps.h"

/* read input file */
/* ================================================================== */
void read_input(char *filename,double *etot,int *lrestart){
  
  char line[MAXLINE],*pt,method[MAXLINE];
  FILE *fpin;
  
  /* initialize */
  *etot = -1.0;
  npathlength = -1;
  dt = -1.0;
  nprintframe = 10;   /* by default skip every 9 frames when writing a path */
  nprintpath = 1;     /* by default skip 0 paths when writing to output     */
  leftbound = -999.0;
  rightbound = 999.0;
  height = -1.0;
  width = -1.0;
  epsilon = -1.0;
  k_restraint = 0.0;
  req_restraint = 0.0;
  npaths = 1;
  *lrestart = 0;
  lrstrtout = 0;
  maxtrail = 100;    /* default max number of TPS trails */
  istartframe = 0;   /* by default start from a randomly chosen frame */ 
  strcpy(restartfilename,"-1");
  
  /* open input file */
  fpin = fopen(filename,"r");
  if(fpin == NULL){
    printf("  ERROR: could not open inputfile %s\n",filename);
    exit(1);
  }
  
  /* read and parse input parameters */
  while(fgets(line,MAXLINE,fpin) != NULL){
    pt = strtok(line," ");
    if(strcmp(pt,"METHOD")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%s",method);
    }else if(strcmp(pt,"ETOT")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%lf",etot);
    }else if(strcmp(pt,"PATHLENGTH")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%d",&npathlength);
    }else if(strcmp(pt,"TIMESTEP")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%lf",&dt);
    }else if(strcmp(pt,"NPRINTFRAME")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%d",&nprintframe);
    }else if(strcmp(pt,"NPRINTPATH")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%d",&nprintpath);
    }else if(strcmp(pt,"RIGHTBOUND")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%lf",&rightbound);
    }else if(strcmp(pt,"LEFTBOUND")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%lf",&leftbound);
    }else if(strcmp(pt,"EPSILON")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%lf",&epsilon);
    }else if(strcmp(pt,"BARRIERHEIGHT")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%lf",&height);
    }else if(strcmp(pt,"BARRIERWIDTH")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%lf",&width);
    }else if(strcmp(pt,"K_RESTRAINT")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%lf",&k_restraint);
    }else if(strcmp(pt,"REQ_RESTRAINT")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%lf",&req_restraint);
    }else if(strcmp(pt,"NPATHS")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%d",&npaths);
    }else if(strcmp(pt,"MAX_TRAIL_TPS")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%d",&maxtrail);
    }else if(strcmp(pt,"RESTART")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%d",lrestart);
    }else if(strcmp(pt,"RESTARTOUT")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%d",&lrstrtout);
    }else if(strcmp(pt,"RESTARTFILE")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%s",restartfilename);
    }else if(strcmp(pt,"RESTARTFROMFRAME")==0){
      pt = strtok(NULL," ");
      sscanf(pt,"%d",&istartframe);
      if(istartframe > 0) ishiftshoot = -1;
    }else{
      printf("  WARNING, unknown keyword: %s\n",pt);
    }
  }

  /* close file and perform some checks */
  fclose(fpin);
  if(strcmp(method,"MD")==0){
    imethod = MD;
    printf("  Performing normal MD simulation\n");
    if(npaths != 1){
      printf("  Funny, you set NPATHS=%d; resetting NPATHS to 1\n",npaths);
      npaths = 1;
      ishiftshoot = 0;  /* don't perform shooting or shifting move */
    }
  }else if(strcmp(method,"TPS")==0){
    imethod = TPS;
    printf("  Performing Transition Path Sampling simulation\n");
    ishiftshoot = 1;  /* start with a shooting move */
  }else if(strcmp(method,"COMMITTOR")==0){
    imethod = COMMITTOR;
    printf("  Performing Committor Analysis simulation\n");
    printf("  - frame number to shoot from: %d\n",istartframe);
    printf("  - number of shooting paths  : %d\n",npaths);
  }else{
    die("in input METHOD should be MD or TPS or COMMITTOR");
  }
  
  if(*etot == -1.0) die("Total energy (ETOT) was not set in input");
  if(npathlength == -1) die("Path length (PATHLENGTH) was not set in input");
  if(dt == -1.0) die("Time step (TIMESTEP) was not set in input");
  if(height == -1.0) die("h (BARRIERHEIGHT) was not set in input");
  if(width == -1.0) die("w (BARRIERWIDTH) was not set in input");
  if(epsilon == -1.0) die("epsilon (EPSILON) was not set in input");
  if(*lrestart){
    if(istartframe==0){
      die("  When using RESTART choose frame with RESTARTFROMFRAME larger than zero");
    }else if(strcmp(restartfilename,"-1") == 0){
      die("  No restart inputfile name specified. Please use RESTARTFILE keyword");
    }else{
      printf("  - Generating initial configuration from restart file %s, frame %d\n",
        restartfilename,istartframe);
    }
  }else{
    printf("  Starting from particles on circle with random velocities\n");
  }
}

/* ================================================================== */
/* read initial transition path from restart file                     */
void read_restart(double xpos[NATOMS],double ypos[NATOMS],double xvel[NATOMS],double yvel[NATOMS]){
    
  int iat, istep, natoms;
  char line[MAXLINE];
  FILE *fpin;
  
  /* open restart file */
  fpin = fopen(restartfilename,"r");
  if(fpin==NULL){
    sprintf(line,"Failed to open restart input file: %s",restartfilename);
    die(line);
  }
  
  /* read header */
  fscanf(fpin,"%d%d\n",&natoms,&nlengthrestart);
  if(natoms!=NATOMS){
    sprintf(line,"Number of atoms in restart input file not equal to %d",NATOMS);
    die(line);
  }
  if(nlengthrestart < istartframe){
    sprintf(line,"RESTARTFROMFRAME is %d but there are only %d frames in restart file",istartframe,nlengthrestart);
    die(line);
  }
  
  /* read frames until istartframe */
  for(istep=0;istep<istartframe;istep++){
    for(iat=0;iat<NATOMS;iat++){
      fscanf(fpin,"%lf%lf%lf%lf",&xpos[iat],&ypos[iat],&xvel[iat],&yvel[iat]);
    }
  }
    
  fclose(fpin);
  return;
}


