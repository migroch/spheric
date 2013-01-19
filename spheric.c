/* spherIC_0.1
**
**  This is a code for generating two component spherical systems in equilibrium. 
**  It can generate haloes from the alpha-beta-gamma family with a stellar component 
**  that follows a King, Plummer or Hernquist density profile. 
**
**  This is an extended version of the halogen4muse code written by Marcel Zemp (mzemp@pku.edu.cn, http://kiaa.pku.edu.cn/~mzemp/). 
**
**  halogen4muse was extended by Miguel Rocha (rocham@uci.edu, http://mrocha.org) to embeed stellar
**  components into the generated haloes.
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "definitions.h"
#include "functions.h"
#include "routines.h"
#include "io.h"

int main(int argc, char **argv) {

  /*
  ** Variables
  */
  
  INT i;
  INT outputgridr,outputgriddf,output_gadget_binary,output_tipsy_binary,output_ifrit,output_profiles;
  DOUBLE randomseed;
  DOUBLE t0, t1, t2, t3, t4, t5, t6, t7,t8;
  PARTICLE *bh;
  SI *si;
  TIPSY_STRUCTURE *ts;
  CHAR FILENAME[STRINGSIZE], INPUTNAME[STRINGSIZE];
  FILE *file;
  
  randomseed = time(NULL);
  
  t0 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
  
  /*
  ** Initialise fixed structures
  */
  
  bh = malloc(sizeof(PARTICLE));
  assert(bh != NULL);  
  si = malloc(sizeof(SI));
  assert(si != NULL);
  si->sp = malloc(sizeof(SP));
  assert(si->sp != NULL);
  si->griddf = malloc(sizeof(GRIDDF));
  assert(si->griddf != NULL);
  si->stuff = malloc(sizeof(STUFF));
  assert(si->stuff != NULL);
  si->gridr = malloc(sizeof(GRIDR));
  assert(si->gridr != NULL);
  ts = malloc(sizeof(TIPSY_STRUCTURE));
  assert(ts != NULL);

  /*
  ** Initialise standard flag values
  */

  outputgridr = 0;
  outputgriddf = 0;
  output_gadget_binary = 0;
  output_tipsy_binary = 0;
  output_ifrit = 0;
  output_profiles = 0;
	
  si->dorvirexact = 0;
  si->halo_flag = 0;
  si->king_flag = 0;
  si->hernquist_flag = 0;
  si->plummer_flag = 0;
  si->stars_flag = 0;
  si->nostarpot_flag = 0;
	
  sprintf(INPUTNAME,"none");

  initialise_black_hole(bh);
  initialise_parameters(si);

  /*
  ** Read in and calculate model parameters
  */

  i = 1;
  while (i < argc) {
    /*
    ** General flags
    */
    if (strcmp(argv[i],"-halo") == 0) {
      si->halo_flag = 1;
      i++;
    }
    else if (strcmp(argv[i],"-king") == 0) {
      si->king_flag = 1;
      si->stars_flag = 1;
      i++;
    }
    else if (strcmp(argv[i],"-hernquist") == 0) {
      si->hernquist_flag = 1;
      si->stars_flag = 1;
      i++;
    }
    else if (strcmp(argv[i],"-plummer") == 0) {
      si->plummer_flag = 1;
      si->stars_flag = 1;
      i++;
    }	
    else if (strcmp(argv[i],"-nostarpot") == 0) {
      si->nostarpot_flag = 1;
      i++;
    }
    /*
    ** Halo parameters
    */
    else if (strcmp(argv[i],"-a") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->alpha = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-b") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->beta = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-c") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->gamma = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-Mhalo") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->M = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-rs") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->rs = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-rcutoff") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->rcutoff = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-Nhalo") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->N0 = (INT) (atof(argv[i]));
      i++;
    }
    else if (strcmp(argv[i],"-soft0") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->soft0 = atof(argv[i]);
      i++;
    }
    /*
    ** Stellar parameters
    */

    else if (strcmp(argv[i],"-Nstar") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->Nstar = atof(argv[i]);
      i++;
    }	
    else if (strcmp(argv[i],"-Mstar") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->Mstar = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-rc") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->rc = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-rt") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->rt = atof(argv[i]);
      i++;
    }	
    else if (strcmp(argv[i],"-rhern") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->rhern = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-rp") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->sp->rp = atof(argv[i]);
      i++;
    }

    /*
    ** Black hole parameters
    */

    else if (strcmp(argv[i],"-MBH") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      bh->mass = atof(argv[i]);
      i++;
    }
	
    /*
    ** Model name
    */

    else if (strcmp(argv[i],"-name") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      sprintf(INPUTNAME,"%s",argv[i]);
      i++;
    }
	
    /*
    ** Output parameters
    */

    else if (strcmp(argv[i],"-ogr") == 0) {
      outputgridr = 1;
      i++;
    }
    else if (strcmp(argv[i],"-ogdf") == 0) {
      outputgriddf = 1;
      i++;
    }
    else if (strcmp(argv[i],"-ogb") == 0) {
      output_gadget_binary = 1;
      i++;
    }
    else if (strcmp(argv[i],"-otb") == 0) {
      output_tipsy_binary = 1;
      i++;
    }
    else if (strcmp(argv[i],"-oift") == 0) {
      output_ifrit = 1;
      i++;
    }
    else if (strcmp(argv[i],"-opfs") == 0) {
      output_profiles = 1;
      i++;
    }	
    else if (strcmp(argv[i],"-dx") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->deltapos[0] = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-dy") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->deltapos[1] = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-dz") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->deltapos[2] = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-dvx") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->deltavel[0] = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-dvy") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->deltavel[1] = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-dvz") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      si->deltavel[2] = atof(argv[i]);
      i++;
    }
    /*
    ** Special parameters
    */

    else if (strcmp(argv[i],"-randomseed") == 0) {
      i++;
      if (i >= argc) {
	usage();
      }
      randomseed = atof(argv[i]);
      i++;
    }
    else if (strcmp(argv[i],"-dorvirexact") == 0) {
      si->dorvirexact = 1;
      i++;
    }
    /*
    ** Version
    */		
    else if (strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"-version") == 0 ){
      fprintf(stderr,"\nversion 1.0");
      exit(1);
    }
    /*
    ** Help or failure
    */
    else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
      usage();
    }
    else {
      fprintf(stderr,"\nERROR: The given command line argument %s was not recognized",argv[i]);
      usage();
    }
  }

  fprintf(stderr,"\n\nChecking parameters, calculating halo properties and initialising grid in r... \n\n");
  srand(randomseed);

  /*
  ** Check main input parameters
  */

  if (strcmp(INPUTNAME,"none") == 0) {
    fprintf(stderr,"You have not set a name for the output model.\n");
    usage();
  }
  if ((NGRIDR-1) % (NGRIDDF-1) != 0) {
    fprintf(stderr,"Bad choice of NGRIDR and NGRIDDF!\n");
    fprintf(stderr,"These numbers have to fulfill the condition (NGRIDR-1) mod (NGRIDDF-1) == 0.\n");
    usage();
  } 
		
  check_main_parameters(si);

  /*
  ** Initialise gridr
  */

  calculate_parameters(si);
  initialise_gridr(bh,si);

  if (outputgridr == 1) {
    sprintf(FILENAME,"%s-gridr.txt",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    write_gridr(file,si);
    fclose(file);
  }

  /*
  ** Initialise griddf
  */

  t1 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
  fprintf(stderr,"Done in "OFD1" seconds.\nInitialising grid for the distribution functions... \n",t1-t0);

  initialise_griddf(si);
  if (outputgriddf == 1) {
    sprintf(FILENAME,"%s-griddf.txt",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    write_griddf(file,si);
    fclose(file);
  }
	
  /*
  ** Initialise structure
  */
  if (si->halo_flag == 1) initialise_haloStructure(si);
  if (si->stars_flag == 1 && si->nostarpot_flag == 0) initialise_starStructure(si);

  /*
  ** Set particle positions
  */

  t2 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
  fprintf(stderr,"Done in "OFD1" seconds.\n\nSetting particle positions... \n",t2-t1);

  if (si->halo_flag == 1) set_haloPositions(si);
  if (si->stars_flag == 1 && si->nostarpot_flag == 0) set_starPositions(si);

  /*
  ** Set particle velocities
  */

  t3 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
  fprintf(stderr,"Done in "OFD1" seconds.\n\nSetting particle velocities... \n",t3-t2);

  if (si->halo_flag == 1) set_haloVelocities(si);
  if (si->stars_flag == 1 && si->nostarpot_flag == 0) set_starVelocities(si);

  /*
  ** Set remaining attributes
  */

  t4 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
  fprintf(stderr,"Done in "OFD1" seconds.\n\nSetting remaining particle attributes... \n",t4-t3);

  set_attributes(si);

  /*
  ** Calculate a few things and do center of mass correction
  */

  t5 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
  fprintf(stderr,"Done in "OFD1" seconds\n\nCalculating a few things, doing mass scaling and correct center of mass position and velocity... \n",t5-t4);
    
  if (si->halo_flag == 1) double_haloParticles(si);	
  if (si->stars_flag == 1 && si->nostarpot_flag == 0) double_starParticles(si);    
  set_index(si);
  calculate_stuff(bh,si);

  /*
  ** Selecting Stars if -nostarpot option is given
  */

  t6 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
  fprintf(stderr,"Done in "OFD1" seconds\n\nSetting displacements... \n",t6-t5);
	
  if (si->stars_flag == 1 && si->nostarpot_flag == 1) {
    fprintf(stderr,"Setting mass to light ratios... \n");	
    set_lighttomass(si);
  }
  displace(si);	

  /*
  ** Write Output
  */

  t7 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
  fprintf(stderr,"Done in "OFD1" seconds\n\nWriting Output... \n",t7-t6);

  if (output_gadget_binary == 1) {
    sprintf(FILENAME,"%s-gadget.bin",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    write_gadget(file,bh,si);
    fclose(file);
  }

  if (output_tipsy_binary == 1) {
    sprintf(FILENAME,"%s-tipsy.bin",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    transfer_particles(bh,si,ts);
    write_tipsy(file,ts);
    fclose(file);
  }

  if (output_ifrit == 1){  
    sprintf(FILENAME,"%s-ifrit.bin",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    write_ifrit(file,si);
    fclose(file);
  }
	
  if (output_profiles == 1){
    sprintf(FILENAME,"%s-profiles.txt",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    write_profiles(file,si,ngrid_profile,rin_profile,rout_profile);
    fclose(file);
  }
	
  if (si->halo_flag == 1){
    sprintf(FILENAME,"%s-haloIC.txt",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    write_ics(file,si,1);
    fclose(file);	
  }

  if (si->stars_flag == 1 && si->nostarpot_flag == 0){
    sprintf(FILENAME,"%s-starsIC.txt",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    write_ics(file,si,2);
    fclose(file);	
  }

	   
  /* 
  ** Print some output in file
  */

  sprintf(FILENAME,"%s.out",INPUTNAME);
  file = fopen(FILENAME,"w");
  assert(file != NULL);
  fprintf(file,"\n\n");
  fprintf(file,"Command line\n\n");
  for (i = 0; i < argc; i++) fprintf(file,"%s ",argv[i]);
  fprintf(file,"\n\n");
  if (si->halo_flag == 1){ 
    fprintf(file,"Halo model properties\n\n");	
    fprintf(file,"alpha = "OFD1"\n",si->sp->alpha);
    fprintf(file,"beta  = "OFD1"\n",si->sp->beta);
    fprintf(file,"gamma = "OFD1"\n",si->sp->gamma);
    fprintf(file,"rho0  = "OFD3" MU LU^-3\n",si->sp->rho0);
    fprintf(file,"rs    = "OFD3" LU\n",si->sp->rs);
    fprintf(file,"rhalf = "OFD3" LU\n",si->sp->rhalf);
    if (si->sp->rcutoff != SBI) {
      fprintf(file,"rcutoff = "OFD3" LU\n",si->sp->rcutoff);
      fprintf(file,"rdecay  = "OFD3" LU\n",si->sp->rdecay);
    }
    fprintf(file,"\n");
  }
  if (si->stars_flag == 1){ 
    fprintf(file,"Stellar model properties\n\n");
    fprintf(file,"K        = "OFD3" MU LU^-3\n",si->sp->K);
    if (si->king_flag == 1){
      fprintf(file,"rt       = "OFD3" LU\n",si->sp->rt);
      fprintf(file,"rc       = "OFD3" LU\n",si->sp->rc);
    } else if (si->plummer_flag == 1) {
      fprintf(file,"rp       = "OFD3" LU\n",si->sp->rp);
    } else if (si->hernquist_flag == 1) {
      fprintf(file,"rhernquist = "OFD3" LU\n",si->sp->rhern);
    }		
    fprintf(file,"rhalfStar= "OFD3" LU\n",si->sp->rhalfStar);
    fprintf(file,"\n");
  }
  if (bh->mass > 0) {
    fprintf(file,"MBH = "OFD3" MU\n",bh->mass);
    fprintf(file,"\n");
  }
  fprintf(file,"System properties\n\n");
  fprintf(file,"rvir    = "OFD3" LU\n",si->sp->rvir);
  fprintf(file,"rinner  = "OFD3" LU\n",si->rinner);
  fprintf(file,"router  = "OFD3" LU\n",si->router);
		
  fprintf(file,"\n");
  if (si->halo_flag == 1){ 
    fprintf(file,"Mhalo(rs)      = "OFD3" MU\n",MencHalo(si->sp->rs,si));
    fprintf(file,"Mhalo(rhalf)   = "OFD3" MU\n",MencHalo(si->sp->rhalf,si));
    fprintf(file,"Mhalo(rvir)    = "OFD3" MU\n",MencHalo(si->sp->rvir,si));
    fprintf(file,"Mhalo(rinner)  = "OFD3" MU\n",si->Mmin);
    fprintf(file,"Mhalo(router)  = "OFD3" MU\n",si->Mmax);	
    if (si->sp->rcutoff != SBI) {
      fprintf(file,"Mhalo(rcutoff) = "OFD3" MU\n",MencHalo(si->sp->rcutoff,si));
    }
    fprintf(file,"MparticleHalo  = "OFD3" MU\n",si->sp->M/si->N);
  }
  if (si->stars_flag == 1) {
    if (si->king_flag == 1) { 
      fprintf(file,"Mstar(rt)      = "OFD3" MU\n",si->MmaxStar);
      fprintf(file,"Mstar(rc)      = "OFD3" MU\n",MencStar(si->sp->rc,si));
    } else if (si->plummer_flag == 1){
      fprintf(file,"Mstar(rp)      = "OFD3" MU\n",MencStar(si->sp->rp,si));
    } else if (si->hernquist_flag == 1){
      fprintf(file,"Mstar(rhernquist) = "OFD3" MU\n",MencStar(si->sp->rhern,si));
    }
    fprintf(file,"Mstar(rhalfStar) = "OFD3" MU\n",MencStar(si->sp->rhalfStar,si));
    fprintf(file,"Mtotal(rhalfStar) = "OFD3" MU\n",Menc(si->sp->rhalfStar,si));
    fprintf(file,"MparticleStar  = "OFD3" MU\n",si->sp->Mstar/si->Nstar);
  }
  fprintf(file,"\n");
  fprintf(file,"Sampling properties\n\n");
  fprintf(file,"|Cr| = "OFD3" LU       Cr = ("OFD4", "OFD4", "OFD4") MU\n",si->stuff->Cr[0],si->stuff->Cr[1],si->stuff->Cr[2],si->stuff->Cr[3]);
  fprintf(file,"|Cv| = "OFD3" LU TU^-1 Cv = ("OFD4", "OFD4", "OFD4") LU TU^-1\n",si->stuff->Cv[0],si->stuff->Cv[1],si->stuff->Cv[2],si->stuff->Cv[3]);
  fprintf(file,"\n");
  fprintf(file,"Etot = "OFD4" MU LU^2 TU^-2\n",si->stuff->Etot);
  fprintf(file,"Ekin = "OFD4" MU LU^2 TU^-2\n",si->stuff->Ekin);
  fprintf(file,"Epot = "OFD4" MU LU^2 TU^-2\n",si->stuff->Epot);
  fprintf(file,"Rvir = |2*Ekin/Epot| = %g\n",fabs(2*si->stuff->Ekin/si->stuff->Epot));
  fprintf(file,"\n");
  fprintf(file,"Ntot                = "OFD3" = "OFI1"\n",(DOUBLE)si->stuff->N,si->stuff->N);
  fprintf(file,"rimp                = "OFD3" LU\n",si->rimp);
  if (si->halo_flag == 1){ 
    fprintf(file,"r1                  = "OFD3" LU\n",si->r1);
    fprintf(file,"r100                = "OFD3" LU\n",si->r100);
  }
  fprintf(file,"Mtheo               = "OFD3" MU\n",bh->mass + si->sp->M + si->sp->Mstar);
  fprintf(file,"Msamp               = "OFD3" MU\n",si->stuff->Mp);
  fprintf(file,"(Msamp-Mtheo)/Mtheo = "OFD3"\n",si->stuff->Mp/(bh->mass + si->sp->M + si->sp->Mstar)-1.0);
  fprintf(file,"Random seed         = "OFD3"\n",randomseed);
  fprintf(file,"\n");
  fprintf(file,"Times for individual steps\n\n");
  fprintf(file,"Calculation of halo properties and initialisation of grid in r: "OFD1" seconds.\n",t1-t0);
  fprintf(file,"Initialisation of grid for distribution function: "OFD1" seconds.\n",t2-t1);
  fprintf(file,"Setting particle positions: "OFD1" seconds\n",t3-t2);
  fprintf(file,"Setting particle velocities: "OFD1" seconds\n",t4-t3);
  fprintf(file,"Setting remaining particle attributes: "OFD1" seconds\n",t5-t4);
  fprintf(file,"Calculating a few things and correct center of mass: "OFD1" seconds\n",t6-t5);
  fprintf(file,"Selecting star particles: "OFD1" seconds!\n",t7-t6);
  t8 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
  fprintf(file,"Writing output: "OFD1" seconds\n",t8-t7);
  fprintf(file,"Total time: "OFD1" seconds\n",t8-t0);
  fclose(file);
   
  fprintf(stderr,"Done in "OFD1" seconds\nTotal time needed was "OFD1" seconds\n",t8-t7,t8-t0);

  free(bh);
  free(si);
  exit(0);

} /* end of main function */

