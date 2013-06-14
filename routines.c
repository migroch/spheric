/* 
** routines.c 
**
** Routines for spherIC_0.1
*/

#include <stdio.h>


/*
** Routine for initialising systems
*/


#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "definitions.h"
#include "functions.h"
#include "routines.h"

void initialise_parameters(SI *si) {
	
  int i;

  si->sp->alpha = -1;
  si->sp->beta = -1;
  si->sp->gamma = -1;
  si->sp->delta = -1;
  si->sp->M = -1;
  si->sp->Mstar = -1;
  si->sp->rs = -1;
  si->sp->rhalf = -1;
  si->sp->rcutoff = -1;
  si->sp->rho0 = -1;
  si->sp->rvir = -1;
  si->sp->rdecay = -1;
  si->sp->rt = -1;
  si->sp->rc = -1;
  si->sp->rp = -1;
  si->sp->rhern = -1;	
  si->sp->rhalfStar = -1;
  si->N0 = 0;
  si->Nstar = 0;
  si->rimp = SBI;
  si->r1 = -1;
  si->r100 = -1;
  si->soft0 = -1;
  for(i = 0; i < 3; i++){
    si->deltapos[i] = 0;
    si->deltavel[i] = 0;
  }
}

/*
** Routine for initialising black hole
*/

void initialise_black_hole(PARTICLE *p) {

  INT i;

  for (i = 0; i < 4; i++) {
    p->r[i] = 0;
    p->v[i] = 0;
  }
  p->mass = 0;
}

/*
** Routine for checking parameters of system
*/

void check_main_parameters(SI *si) {

  if (si->stars_flag == 0 && si->halo_flag == 0){
    fprintf(stderr,"Missing or bad parameter:\n");
    fprintf(stderr,"You have to set at least one of the system options ON (set -halo and/or -king/-plummer/-hernquist).\n");	
    usage();
  }
  if (si->nostarpot_flag == 1 && si->halo_flag == 0){
    fprintf(stderr,"Missing or bad parameter:\n");
    fprintf(stderr,"You have set the -nonstarpot flag ON and the -halo flag OFF. You can't exlude the stellar potenteial without having a halo system.\n");	
    usage();
  }
  if (si->halo_flag == 1)
    {	
      if (si->sp->alpha == -1) {
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You have not set a value for alpha.\n");
	usage();
      }
      if (si->sp->beta == -1) {
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You have not set a value for beta.\n");
	usage();
      }
      if (si->sp->gamma == -1) {
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You have not set a value for gamma.\n");
	usage();
      }
      if (si->sp->gamma >= 3) {
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You have chosen gamma = "OFD1".\n",si->sp->gamma);
	fprintf(stderr,"This means your cumulative mass function is diverging at the center.\n");
	fprintf(stderr,"Use a smaller value for gamma.\n");
	usage();
      }
      if (si->sp->M == -1) {	
	fprintf(stderr,"WARNING: You have not set a value for the halo mass (-Mhalo).\n Mhalo has been set to 1\n");
	si->sp->M = 1;
      }
      if (si->N0 == 0) {
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You have not set a value for the number of paricles in the halo (-Nhalo).\n");
	usage();
      }
    }
	
  if (si->stars_flag == 1)
    {	
      if (si->sp->Mstar == -1) {	
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You have not set a value for the stellar mass (-Mstar).\n");
	usage();
      }	
      if ((si->nostarpot_flag == 1 ) && (si->Nstar != 0)) {
	fprintf(stderr,"Warning: You have asked not to include the stellar potential in the calculation, hence you can not set Nstar your self. Nstar will be ignored.\n");
	si->Nstar = 0;
      }
      if ((si->nostarpot_flag == 0 ) && (si->Nstar == 0)) {
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You want to include the stellar potential in the calculation, hence you have to set the number of star particles (-Nstar) your self .\n");
	usage();
      }
      if (((si->king_flag == 1 ) && (si->plummer_flag == 1)) || ((si->king_flag == 1 ) && (si->hernquist_flag == 1)) || ((si->plummer_flag == 1 ) && (si->hernquist_flag == 1))) {		
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You can only use one of the stellar profile options (-king, -plummer, or -hernquist)\n");
	usage();
      }

      if(si->king_flag == 1){
	if (si->sp->rc == -1) {
	  fprintf(stderr,"Missing or bad parameter:\n");
	  fprintf(stderr,"You have not set a value for the core radius Rc defining the King model (-rc)\n");
	  usage();
	}
	if (si->sp->rt == -1) {
	  fprintf(stderr,"Missing or bad parameter:\n");
	  fprintf(stderr,"You have not set a value for the tidal radius Rt defining the King model (-rt)\n");
	  usage();
	}
      }
      if(si->plummer_flag == 1){
	if(si->sp->rp == -1){
	  fprintf(stderr,"Missing or bad parameter:\n");
	  fprintf(stderr,"You have not set a value for the scale radius defining the Plummer model (-rp)\n");
	  usage();
	}
      }
      if(si->hernquist_flag == 1){
	if(si->sp->rhern == -1){
	  fprintf(stderr,"Missing or bad parameter:\n");
	  fprintf(stderr,"You have not set a value for the scale radius defining the Hernquist model (-rp)\n");
	  usage();
	}
      }
    }
		
  //     if (si->soft0 == -1) {
  //         fprintf(stderr,"Missing or bad parameter:);
  //         fprintf(stderr,"You have not set a value for the softening of the particles soft0.\n");
  //         usage();
  //         }
}

/*
** Routine for calculating system parameters
*/

void calculate_parameters(SI *si) {

  DOUBLE IM, IMcutoff,IMStar;

  if (si->halo_flag == 1)
    {
      if (si->sp->beta > 3) {
	/*
	** Finite mass models
	*/
	if (si->sp->rs == -1) {
	  fprintf(stderr,"Missing or bad parameter:\n");
	  fprintf(stderr,"For finite mass models you have to set a value for the scale radius rs.\n");
	  usage();
	}
	if (si->sp->rcutoff != -1) {
	  fprintf(stderr,"Warning: ");
	  fprintf(stderr,"For finite mass models the cutoff radius rcutoff is not needed!\n");
	  fprintf(stderr,"Hence, your input for the cutoff radius rcutoff (= "OFD1" LU) was ignored.\n",si->sp->rcutoff);
	}
	si->sp->rdecay = 0;
	si->sp->delta = 0;
	IM = exp(lgamma((si->sp->beta-3)/si->sp->alpha));
	IM *= exp(lgamma((3-si->sp->gamma)/si->sp->alpha));
	IM /= (si->sp->alpha*exp(lgamma((si->sp->beta-si->sp->gamma)/si->sp->alpha)));
	si->sp->rho0 = si->sp->M/(4*M_PI*(si->sp->rs*si->sp->rs*si->sp->rs)*IM);
	si->sp->rcutoff = SBI;
      }
      else {
	/*
	** Cutoff models
	*/
	if ((si->sp->rs == -1) || (si->sp->rcutoff == -1)) {
	  fprintf(stderr,"Missing or bad parameter:\n");
	  fprintf(stderr,"Specify values for the scale radius -rs and cutoff radius -rcutoff for models with cutoff.\n");
	  usage();
	}
	si->sp->rdecay = CutoffFac*si->sp->rcutoff;
	si->sp->delta = si->sp->rcutoff/si->sp->rdecay + dlrhodlr(si->sp->rcutoff,si);
	IM = pow(1e-6,3-si->sp->gamma)/(3-si->sp->gamma); /* approximate inner integral */
	IM += integral(integrandIM,1e-6,si->sp->rcutoff/si->sp->rs,si);
	IMcutoff = 1/tau(si->sp->rcutoff,si);
	IMcutoff *= 1/(si->sp->rs*si->sp->rs*si->sp->rs);
	IMcutoff *= integral(integrandIMcutoff,si->sp->rcutoff,si->sp->rcutoff+2000*si->sp->rdecay,si);
	IM += IMcutoff;
	si->sp->rho0 = si->sp->M/(4*M_PI*(si->sp->rs*si->sp->rs*si->sp->rs)*IM);
      }		

	
      si->rinner = FACTORRINNER*si->sp->rs;
      si->router = si->sp->rs;
      while (rho(si->sp->rs,si)/rho(si->router,si) < FACTORROUTER) {
	si->router = si->router*10;
      }

      fprintf(stderr,"\nHere are some given and calculated parameters for the halo (for more look at the output file!):\n\n");
      fprintf(stderr,"Nhalo     = "OFI1"\n",si->N0);
      fprintf(stderr,"haloPmass = "OFD3" MU\n",si->sp->M/si->N0);
      fprintf(stderr,"Mhalo     = "OFD3" MU\n",si->sp->M);
      fprintf(stderr,"rho0      = "OFD3" MU LU^-3  -> Normalization constant for the halo model\n",si->sp->rho0);
      fprintf(stderr,"alpha     = "OFD1"\n",si->sp->alpha);
      fprintf(stderr,"beta      = "OFD1"\n",si->sp->beta);
      fprintf(stderr,"gamma     = "OFD1"\n",si->sp->gamma);
      fprintf(stderr,"rs        = "OFD3" LU\n",si->sp->rs);
      if (si->sp->rcutoff != SBI) {
	fprintf(stderr,"rcutoff   = "OFD3" LU\n",si->sp->rcutoff);
	fprintf(stderr,"rdecay    = "OFD3" LU\n\n",si->sp->rdecay);
      }
      else fprintf(stderr,"\n");
    }	
	
  /*
  ** Stellar Models
  */

  if (si->stars_flag == 1) {	

    if (si->king_flag == 1){
      si->rinner = FACTORRINNER*si->sp->rc;
      si->routerStar = si->sp->rt*0.9999;
    } else if (si->plummer_flag == 1){
      si->rinner = FACTORRINNER*si->sp->rp;
      si->routerStar = si->sp->rp;
      while (rhoStar(si->sp->rp,si)/rhoStar(si->routerStar,si) < FACTORROUTERSTAR) {
	si->routerStar = si->routerStar*1.1;
      }
    } else if (si->hernquist_flag == 1){
      si->rinner = FACTORRINNER*si->sp->rhern;
      if(si-> halo_flag == 1){
	si->routerStar = si->router;
      } else 	si->routerStar = si->sp->rhern*500;
      //		while (rhoStar(si->sp->rhern,si)/rhoStar(si->routerStar,si) < FACTORROUTERSTAR) {
      //		si->routerStar = si->routerStar*1.1;
      //		}
    }
    if (si-> halo_flag == 0) si->router = si->routerStar;
	
    IMStar = MencInnerStar(1e-6,si);
    IMStar += integral(integrandMencStar,1e-6,si->routerStar,si);
    si->sp->K = si->sp->Mstar/IMStar;
    fprintf(stderr,"Here are some given and calculated parameters for the stellar model (for more look at the output file!):\n\n");
    fprintf(stderr,"Nstar     = "OFI1"\n",si->Nstar);
    fprintf(stderr,"starPmass = "OFD3" MU\n",si->sp->Mstar/si->Nstar);
    fprintf(stderr,"Mstar     = "OFD3" MU\n",si->sp->Mstar);	
    fprintf(stderr,"K         = "OFD3" MU LU^-3  ->  Normaliztion constant for the stellar model. \n",si->sp->K);
    if (si->king_flag == 1) {
      fprintf(stderr,"rt        = "OFD3" LU\n",si->sp->rt);
      fprintf(stderr,"rc        = "OFD3" LU\n\n",si->sp->rc);
    }
    if (si->plummer_flag == 1) {
      fprintf(stderr,"rp        = "OFD3" LU\n\n",si->sp->rp);
    }
    if (si->hernquist_flag == 1) {
      fprintf(stderr,"rhern     = "OFD3" LU\n\n",si->sp->rhern);
    }
  }
	
  //fprintf(stderr,"rho(rs),rho(router) = %g, %g\n", rho(si->sp->rs,si),rho(si->router,si));
  //fprintf(stderr,"router = %g, routerStar = %g ,rinner = %g\n\n" ,si->router,si->routerStar,si->rinner);

}

/* 
** Routine for initialising gridr 
*/

void initialise_gridr( PARTICLE *bh, SI *si) {

  INT i,hf,sf;
  DOUBLE dlogr, logr, r, r3;
  DOUBLE MencHalor, DeltaMencHalor,MencStarr;
  DOUBLE rhoHalor,rhoStarr;
  DOUBLE rhoencHalor,rhoencStarr;
  DOUBLE Potr, Potoutr;
  GRIDR *gridr;
  SP *sp;

  gridr = si->gridr;
  sp = si->sp;
  sf = si->stars_flag;
  hf = si->halo_flag;
		
  dlogr = (log(si->router)-log(si->rinner))/(NGRIDR-1);
   	
  i = 0;
  logr = log(si->rinner);
  r = exp(logr);
  r3 = r*r*r;
  gridr->logr[i] = logr;
  si->logr[i] = logr;
  gridr->r[i] = r;
	
  if(hf == 1){
    rhoHalor = rho(r,si);
    DeltaMencHalor = 4*M_PI*rhoHalor*r3/(3-sp->gamma); /* approximate inner gridpoint by analytical calculation */
    MencHalor = DeltaMencHalor;
    rhoencHalor = MencHalor/(4*M_PI*r3/3.0);
    gridr->rhoHalo[i] = rhoHalor;
    gridr->logrhoHalo[i] = log(rhoHalor);
    gridr->rhoencHalo[i] = rhoencHalor;
    gridr->logrhoencHalo[i] = log(rhoencHalor);
    gridr->MencHalo[i] = MencHalor;
    gridr->logMencHalo[i] = log(MencHalor);
  } else {
    gridr->rhoHalo[i] = 0;
    gridr->logrhoHalo[i] = 0;
    gridr->rhoencHalo[i] = 0;
    gridr->logrhoencHalo[i] = 0;
    gridr->MencHalo[i] = 0;
    gridr->logMencHalo[i] = 0;	
  }
	   
	  
  if (sf == 1) {
    rhoStarr = sp->K*rhoStar(r,si);
    MencStarr = sp->K*MencInnerStar(r,si);          /* approximate inner gridpoint by analytical calculation */
    rhoencStarr = MencStarr/(4*M_PI*r3/3.0);
    gridr->rhoStar[i] = rhoStarr;
    gridr->logrhoStar[i] = log(rhoStarr);
    gridr->rhoencStar[i] = rhoencStarr;
    gridr->logrhoencStar[i] = log(rhoencStarr);
    gridr->MencStar[i] = MencStarr;
    gridr->logMencStar[i] = log(MencStarr);
  } else {
    gridr->rhoStar[i] = 0;
    gridr->logrhoStar[i] = 0;
    gridr->rhoencStar[i] = 0;
    gridr->logrhoencStar[i] = 0;
    gridr->MencStar[i] = 0;
    gridr->logMencStar[i] = 0;
  }	
		
  gridr->rho[i] = gridr->rhoHalo[i] + gridr->rhoStar[i];
  gridr->logrho[i] = log(gridr->rho[i]); 
  gridr->rhoenc[i] = gridr->rhoencHalo[i] + gridr->rhoencStar[i];
  gridr->logrhoenc[i] = log(gridr->rhoenc[i]);
  si->logrhoenc[i] = log(gridr->rhoenc[i]);
  gridr->Menc[i] = bh->mass + gridr->MencHalo[i] + gridr->MencStar[i]; 
  gridr->logMenc[i] = log(gridr->Menc[i]);
  si->logMenc[i] = log(gridr->Menc[i]);
  gridr->eqrvcmax[i] = gridr->Menc[i] - 4*M_PI*gridr->rho[i]*r3; 
	 
  for (i = 1; i < NGRIDR; i++) {
	
    logr = log(si->rinner) + i*dlogr;
    r = exp(logr);
    r3 = r*r*r;
    gridr->logr[i] = logr;
    si->logr[i] = logr;
    gridr->r[i] = r;
	
		
    if(hf == 1){
      rhoHalor = rho(r,si);
      DeltaMencHalor = integral(integrandMenc,gridr->r[i-1],r,si);
      MencHalor = gridr->MencHalo[i-1] + DeltaMencHalor;
      rhoencHalor = MencHalor/(4*M_PI*r3/3.0);
      gridr->rhoHalo[i] = rhoHalor;
      gridr->logrhoHalo[i] = log(rhoHalor);
      gridr->rhoencHalo[i] = rhoencHalor;
      gridr->logrhoencHalo[i] = log(rhoencHalor);
      gridr->MencHalo[i] = MencHalor;
      gridr->logMencHalo[i] = log(MencHalor);
    } else {
      gridr->rhoHalo[i] = 0;
      gridr->logrhoHalo[i] = 0;
      gridr->rhoencHalo[i] = 0;
      gridr->logrhoencHalo[i] = 0;
      gridr->MencHalo[i] = 0;
      gridr->logMencHalo[i] = 0;	
    }
	   
    if (sf == 1) {
      rhoStarr = sp->K*rhoStar(r,si);
      MencStarr = (r <= si->routerStar) ? gridr->MencStar[i-1] + IMencStar(gridr->r[i-1],r,si) : gridr->MencStar[i-1] ;
      rhoencStarr = MencStarr/(4*M_PI*r3/3.0);
      gridr->rhoStar[i] = rhoStarr;
      gridr->logrhoStar[i] = log(rhoStarr);
      gridr->rhoencStar[i] = rhoencStarr;
      gridr->logrhoencStar[i] = log(rhoencStarr);
      gridr->MencStar[i] = MencStarr;
      gridr->logMencStar[i] = log(MencStarr);
    } else {
      gridr->rhoStar[i] = 0;
      gridr->logrhoStar[i] = 0;
      gridr->rhoencStar[i] = 0;
      gridr->logrhoencStar[i] = 0;
      gridr->MencStar[i] = 0;
      gridr->logMencStar[i] = 0;
    }	
		
    gridr->rho[i] = gridr->rhoHalo[i] + gridr->rhoStar[i];
    gridr->logrho[i] = log(gridr->rho[i]); 
    gridr->rhoenc[i] = gridr->rhoencHalo[i] + gridr->rhoencStar[i];
    gridr->logrhoenc[i] = log(gridr->rhoenc[i]);
    si->logrhoenc[i] = log(gridr->rhoenc[i]);
    gridr->Menc[i] = gridr->MencHalo[i] + gridr->MencStar[i]; 
    gridr->logMenc[i] = log(gridr->Menc[i]);
    si->logMenc[i] = log(gridr->Menc[i]);
    gridr->eqrvcmax[i] = gridr->Menc[i] - 4*M_PI*gridr->rho[i]*r3; 
  }
	
  i = NGRIDR-1;
  Potoutr = 0; /* approximate outer gridpoint by 0 */
  if(si->nostarpot_flag == 1 || si->stars_flag == 0){
    Potr = (-1)*G*(gridr->MencHalo[i]/gridr->r[i] + Potoutr);
  } else if (si->halo_flag == 1){
    Potr = (-1)*G*(gridr->Menc[i]/gridr->r[i]+Potoutr);
  }else{	
    Potr = (-1)*G*(gridr->MencStar[i]/gridr->r[i] + Potoutr);
  }
  gridr->Pot[i] = Potr;
  gridr->logPot[i] = log(-Potr);
  gridr->Potoutr[i] = Potoutr;
  for (i = (NGRIDR-2); i >= 0; i--) {
    if(si->nostarpot_flag == 1 || si->stars_flag == 0){
      Potoutr = integral(integrandPotHalo,gridr->r[i],gridr->r[i+1],si) + gridr->Potoutr[i+1];
      Potr = (-1)*G*(gridr->MencHalo[i]/gridr->r[i] + Potoutr);
    } else if (si->halo_flag == 1){
      Potoutr = integral(integrandPotHalo,gridr->r[i],gridr->r[i+1],si) + integral(integrandPotStar,gridr->r[i],gridr->r[i+1],si) + gridr->Potoutr[i+1];
      Potr = (-1)*G*(gridr->Menc[i]/gridr->r[i] + Potoutr);
    } else if (si->halo_flag == 0){	
      Potoutr = integral(integrandPotStar,gridr->r[i],gridr->r[i+1],si) + gridr->Potoutr[i+1];
      Potr = (-1)*G*(gridr->MencStar[i]/gridr->r[i] + Potoutr);	
    }
    gridr->Pot[i] = Potr;
    gridr->logPot[i] = log(-Potr);
    gridr->Potoutr[i] = Potoutr;
  }
}

/*
** Routine for initialising griddf
*/

void initialise_griddf(SI *si) {

  INT k,i, j, dj,warning_flag,count;
  DOUBLE fE,fEstar,gEj,IdMdE,IdMstardE,rnegative;
  GRIDR *gridr;
  GRIDDF *griddf;

  warning_flag = 0;
  count = 0;
  rnegative = 0;
  IdMdE = 0;
  IdMstardE = 0;
  gridr = si->gridr;
  griddf = si->griddf;
  dj = (NGRIDR-1) / (NGRIDDF-1);
	
  for (i = (NGRIDDF-1), j = (NGRIDR-1) ; i >= 0 ; i--, j = j-dj) {
		
    if(i == NGRIDDF-1){
      fprintf(stderr,"\n");
      fprintf(stderr,"[");
      for(k = 0; k < 50; k++) fprintf(stderr," ");        //All this is jut to have a progress bar..
      fprintf(stderr,"]");
      fprintf(stderr,"\b");
      for(k = 0; k < 50; k++) fprintf(stderr,"\b");
    }
	
    k = (NGRIDDF-1)/50;
    if(count == k){
      count = 0;
      fprintf(stderr,".");
    }
    count++;
	
    //	rotate();
    //  fprintf(stderr,"\b\b\b\b%4d",i);                    //an alternative progress indicator..
	
	
    //Here we calculate the Distribution Function and Energy Distribution for the different modes of operation.

    //There are 4 possible modes of operation. The first one is a "halo + star" system that excludes the stellar potential. 
    //We can take care of this mode together with the "just a halo" and no stars mode.	
		
    if(si->nostarpot_flag == 1 || si->stars_flag == 0){
		
      fE = integraldf(d2rhodPhi2,j,si);
      if (fE < 0) {
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You have chosen parameters that lead to a negative and hence unphysical distribution function for the halo system.\n");
	usage();	
      }
	
      gEj = (i > 0) ? gE(gridr->Pot[j],gridr->r[j],si) : 0.0;
	
      griddf->r[i] = gridr->r[j];
      griddf->logr[i] = gridr->logr[j];
      griddf->E[i] = gridr->Pot[j];
      griddf->logE[i] = gridr->logPot[j];
      griddf->fE[i] = fE;
      griddf->logfE[i] = log(fE);
      griddf->gE[i] = gEj;
      griddf->loggE[i] = log(gEj);
	
      if (si->stars_flag == 1) {
	fEstar = integraldf(d2rhoStardPhi2,j,si);
	if (fEstar < 0) {
	  fEstar = 0.0;
	  if (warning_flag == 0){
	    rnegative = gridr->r[j]; 
	    warning_flag = 1;
	  }	
	}
	griddf->fEstar[i] = fEstar;
	griddf->logfEstar[i] = log(fEstar);
      }
      else{
	griddf->fEstar[i] = 0.0;
	griddf->logfEstar[i] = 0.0;		
      }	
    }
		
    //Else, we can have a halo + stars full potentential model.

    else if (si->halo_flag == 1){
					
      fE = integraldf(d2rhodPhi2,j,si);
      if (fE < 0) {
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You have chosen parameters that lead to a negative and hence unphysical distribution function for the halo system.\n");
	usage();	
      }
		
      fEstar = integraldf(d2rhoStardPhi2,j,si);
      if (fEstar < 0) {
	fEstar = 0.0;
	if (warning_flag == 0){
	  rnegative = gridr->r[j]; 
	  warning_flag = 1;
	}			
      }
		
      gEj = (i > 0) ? gE(gridr->Pot[j],gridr->r[j],si) : 0.0;
		
      griddf->r[i] = gridr->r[j];
      griddf->logr[i] = gridr->logr[j];
      griddf->E[i] = gridr->Pot[j];
      griddf->logE[i] = gridr->logPot[j];
      griddf->fE[i] = fE;
      griddf->logfE[i] = log(fE);
      griddf->gE[i] = gEj;
      griddf->loggE[i] = log(gEj);
      griddf->fEstar[i] = fEstar;
      griddf->logfEstar[i] = log(fEstar);
		
    }		

    //Lastly we can have just an stellar system model.

    else{
				
      fE = 0.0;			
			
      fEstar = integraldf(d2rhoStardPhi2,j,si);
      if (fEstar < 0) {
	fprintf(stderr,"Missing or bad parameter:\n");
	fprintf(stderr,"You have chosen parameters that lead to a negative and hence unphysical distribution function for the stellar system.\n");
	usage();	
      }				
							
      gEj = (i > 0) ? gE(gridr->Pot[j],gridr->r[j],si) : 0.0;
	
      griddf->r[i] = gridr->r[j];
      griddf->logr[i] = gridr->logr[j];
      griddf->E[i] = gridr->Pot[j];
      griddf->logE[i] = gridr->logPot[j];
      griddf->fE[i] = fE;
      griddf->logfE[i] = log(fE);
      griddf->gE[i] = gEj;
      griddf->loggE[i] = log(gEj);
      griddf->fEstar[i] = fEstar;
      griddf->logfEstar[i] = log(fEstar);
    }
				
	
    if (i == 0) fprintf(stderr,"\n");
		
  }	
	
  if (warning_flag == 1){
    fprintf(stderr,"WARNING: Negative and hence unphysical values of the stellar distribution function were found at one or more points on the grid. These were set to 0.\n");
    fprintf(stderr,"Negative values of the stellar distribution function were found for r <= %lf LU.\n\n",rnegative);
  }

  //Uncomment this if you want to check the mass obtained by integrating the distribution function.

  //if (si->halo_flag == 1) { 
    //IdMdE = integral(dMdE,griddf->E[0],griddf->E[NGRIDDF-1],si);
    //fprintf(stderr,"Mhalo from integrating dM = fhalo(E)g(E)dE : %5.4e\n",IdMdE);
    //}
  //if (si->stars_flag == 1) {
    // IdMstardE = integral(dMstardE,griddf->E[0],Pot(si->routerStar,si),si);
    //fprintf(stderr,"Mstar from integrating dM = fstar(E)g(E)dE : %5.4e\n",IdMstardE);
    //  }

}

/*
** Routines for initialising halo and stellar structures
*/

void initialise_haloStructure(SI *si) {
    
  GRIDR *gridr;
	
  gridr = si->gridr;
  si->N = si->N0/2;
  si->p = malloc(si->N*sizeof(PARTICLE));
  si->Mmin = gridr->MencHalo[0]; 
  si->Mmax = gridr->MencHalo[NGRIDR-1];
  si->mass = si->sp->M/(2.0*si->N);
}


void initialise_starStructure(SI *si) {
	
  GRIDR *gridr;
	
  gridr = si->gridr;
  si->Nstar = si->Nstar/2;
  si->pstar = malloc(si->Nstar*sizeof(PARTICLE));
  si->MminStar = gridr->MencStar[0];
  if(si->halo_flag == 1){
    si->MmaxStar = exp(lininterpolate(NGRIDR,gridr->logr,gridr->logMencStar,log(si->routerStar)));
  } else {
    si->MmaxStar = gridr->MencStar[NGRIDR-1];
  }	
  si->massStar = si->sp->Mstar/(2.0*si->Nstar);
}
/*
** Routine for setting positions of halo particles
*/

void set_haloPositions(SI *si) {
    
  INT i, N;
  DOUBLE Mrand, logMrand, Mmin, Mmax;
  DOUBLE rrand, logrrand;
  DOUBLE theta, phi;
  PARTICLE *p;
  GRIDR *gridr;
	
  gridr = si->gridr;
  N = si->N;
  p = si->p;
  Mmin = si->Mmin;
  Mmax = si->Mmax;
  for (i = 0; i < N; i++) {
    Mrand = Mmin + rand01()*(Mmax - Mmin);
    logMrand = log(Mrand);
    logrrand = lininterpolate(NGRIDR,gridr->logMencHalo,gridr->logr,logMrand);
    rrand = exp(logrrand);
    theta = acos(2.0*rand01() - 1.0);
    phi = rand01()*2.0*M_PI;
    p[i].r[0] = rrand;
    p[i].r[1] = rrand*sin(theta)*cos(phi);
    p[i].r[2] = rrand*sin(theta)*sin(phi);
    p[i].r[3] = rrand*cos(theta);
  }
}

/*
** Routine for setting positions of star particles
*/

void set_starPositions(SI *si) {
    
  INT i, N;
  DOUBLE Mrand, logMrand, Mmin, Mmax;
  DOUBLE rrand, logrrand;
  DOUBLE theta, phi;
  PARTICLE *p;
  GRIDR *gridr;
	
  gridr = si->gridr;
  N = si->Nstar;
  p = si->pstar;
  Mmin = si->MminStar;
  Mmax = si->MmaxStar;
  for (i = 0; i < N; i++) {
    Mrand = Mmin + rand01()*(Mmax - Mmin);
    logMrand = log(Mrand);
    logrrand = lininterpolate(NGRIDR,gridr->logMencStar,gridr->logr,logMrand);
    rrand = exp(logrrand);
    theta = acos(2.0*rand01() - 1.0);
    phi = rand01()*2.0*M_PI;
    p[i].r[0] = rrand;
    p[i].r[1] = rrand*sin(theta)*cos(phi);
    p[i].r[2] = rrand*sin(theta)*sin(phi);
    p[i].r[3] = rrand*cos(theta);
  }
}
/*
** Routine for setting velocities of halo particles 
*/

void set_haloVelocities(SI *si) {

  INT i, N;
  DOUBLE r, Erand, Potr;
  DOUBLE fEmax, fErand, fEcheck;
  DOUBLE vesc, vrand;
  DOUBLE theta, phi;
  GRIDR *gridr;
  GRIDDF *griddf;
  PARTICLE *p;

  gridr = si->gridr;
  griddf = si->griddf;
  N = si->N;
  p = si->p;
  for (i = 0; i < N; i++) {
    r = p[i].r[0];
    Potr = Pot(r,si);
    vesc = vescape(r,si);
    fEmax = f2(r,si);
    vrand = 0;
    Erand = 0;
    fErand = 0;
    fEcheck = 1;
    while (fEcheck > fErand) {
      vrand = pow(rand01(),1.0/3.0)*vesc;
      Erand = 0.5*vrand*vrand + Potr;
      fErand = f1(Erand,si);
      fEcheck = rand01()*fEmax;
    }
    theta = acos(2.0*rand01() - 1.0);
    phi = rand01()*2.0*M_PI;
    p[i].v[0] = vrand;
    p[i].v[1] = vrand*sin(theta)*cos(phi);
    p[i].v[2] = vrand*sin(theta)*sin(phi);
    p[i].v[3] = vrand*cos(theta);
    p[i].E = Erand;
  }
}
	
/*
** Routine for setting velocities of star particles 
*/

void set_starVelocities(SI *si) {

  INT i, N;
  DOUBLE r, Erand, Potr; // ,vr;
  DOUBLE fEmax, fErand, fEcheck;
  DOUBLE vesc, vrand;
  DOUBLE theta, phi;
  GRIDR *gridr;
  GRIDDF *griddf;
  PARTICLE *p;

  gridr = si->gridr;
  griddf = si->griddf;
  N = si->Nstar;
  p = si->pstar;
  for (i = 0; i < N; i++) {
    r = p[i].r[0];
    Potr = Pot(r,si);
    vesc = vescapeStar(r,si);
    fEmax = f2star(r,si);
    /*		if( r < si->routerStar &&  fEmax == 0 ){
		vrand = integral(integrandJeansStar,r,si->routerStar,si)/(si->sp->K*rhoStar(r,si));  // this is Vr^2 from the jeans equation. We are forcing the velocities to be what we need to stay in equilibrium.
		vr = sqrt(vrand);
		p[i].v[1] = rand01()*(vr*r/p[i].r[1]);
		p[i].v[2] = rand01()*(vr*r - p[i].v[1]*p[i].r[1])/p[i].r[2];		
		p[i].v[3] = rand01()*(vr*r - p[i].v[1]*p[i].r[1] - p[i].v[2]*p[i].r[2])/p[i].r[3];
		p[i].v[0] = sqrt(p[i].v[1]*p[i].v[1] + p[i].v[2]*p[i].v[2] + p[i].v[3]*p[i].v[3]);
		Erand = 0.5*p[i].v[0]*p[i].v[0]  + Potr;
		p[i].E = Erand;
		}
		else{*/
    vrand = 0;
    Erand = 0;
    fErand = 0;
    fEcheck = 1;
    while (fEcheck > fErand) {
      vrand = pow(rand01(),1.0/3.0)*vesc;
      Erand = 0.5*vrand*vrand + Potr;
      fErand = f1star(Erand,si);
      fEcheck = rand01()*fEmax;
    }
    theta = acos(2.0*rand01() - 1.0);
    phi = rand01()*2.0*M_PI;
    p[i].v[0] = vrand;
    p[i].v[1] = vrand*sin(theta)*cos(phi);
    p[i].v[2] = vrand*sin(theta)*sin(phi);
    p[i].v[3] = vrand*cos(theta);
    p[i].E = Erand;
    //		}
  }
}	
	
	
/*
** Routine for doubling halo particles with mirror halo
*/

void double_haloParticles(SI *si) {

  INT i, N;
  PARTICLE *p;

  N = si->N;
  p = si->p;
  p = realloc(p,2*N*sizeof(PARTICLE));
  if (N > 0) {
    assert(p != NULL);
  }
  for (i = 0; i < N; i++) {
    p[N+i].r[0] = p[i].r[0];
    p[N+i].r[1] = -p[i].r[1];
    p[N+i].r[2] = -p[i].r[2];
    p[N+i].r[3] = -p[i].r[3];
    p[N+i].v[0] =  p[i].v[0];
    p[N+i].v[1] = -p[i].v[1];
    p[N+i].v[2] = -p[i].v[2];
    p[N+i].v[3] = -p[i].v[3];
    p[N+i].mass = p[i].mass;
    p[N+i].E = p[i].E;
    p[i].star_flag = 0;
    p[N+i].star_flag = p[i].star_flag;
		 
  }
  si->N = 2*N;
  si->p = p;
}

/*
** Routine for doubling star particles with mirror stellar system
*/

void double_starParticles(SI *si) {

  INT i, N;
  PARTICLE *p;

  N = si->Nstar;
  p = si->pstar;
  p = realloc(p,2*N*sizeof(PARTICLE));
  if (N > 0) {
    assert(p != NULL);
  }
  for (i = 0; i < N; i++) {
    p[N+i].r[0] = p[i].r[0];
    p[N+i].r[1] = -p[i].r[1];
    p[N+i].r[2] = -p[i].r[2];
    p[N+i].r[3] = -p[i].r[3];
    p[N+i].v[0] =  p[i].v[0];
    p[N+i].v[1] = -p[i].v[1];
    p[N+i].v[2] = -p[i].v[2];
    p[N+i].v[3] = -p[i].v[3];
    p[N+i].mass = p[i].mass;
    p[N+i].E = p[i].E;
    p[i].star_flag = 1;
    p[N+i].star_flag = p[i].star_flag;
		 
  }
  si->Nstar = 2*N;
  si->pstar = p;
}

void set_index(SI *si){

  INT i, N,Nstar;
  PARTICLE *p;
  PARTICLE *pstar;
    
  if(si->halo_flag == 1){
    N = si->N;
    p = si->p;
    for (i = 0; i < N; i++) {
      p[i].index = i + 1;
    }
  }
  else N = 0;	


  if(si->stars_flag == 1 && si->nostarpot_flag == 0){
    Nstar = si->Nstar;
    pstar = si->pstar;
    for (i = 0; i < Nstar; i++) {
      pstar[i].index =  i + N + 1 ;
    }
  }



}

/*
** Routine for assigning light to mass (f(E)star/f(E)dark) raios to halo particles for when the stellar potential is
** not included in the calculation.
*/

void set_lighttomass(SI *si){

  INT i,j,N;
  DOUBLE fs_f;
  GRIDDF *griddf;
  PARTICLE *p;
  N = si->N;
  p = si->p;
  griddf = si->griddf;
		
  for( i = 0; i < NGRIDDF-1 ; i++ ){
    fs_f = (griddf->fEstar[i])/(griddf->fE[i]);
    if (fs_f > 0){
      for(j = 0; j < N; j++){		
	if ((griddf->E[i]) <= p[j].E && (p[j].E < griddf->E[i+1])){
	  p[j].star_flag = fs_f;
	}
      }
    }
  }
}


/*
** Routine for tagging halo particles as stars when the stellar potential is
** not included in the calculation.(This routine is not used in this version)
*/

void select_stars(SI *si){

  INT i,j,starCount[NGRIDDF-1],Nstar,N,warning_flag ;
  DOUBLE nstar,fs_f,M;
  GRIDDF *griddf;
  PARTICLE *p;
  N = si->N;
  M = si->sp->M;
  p = si->p;
  griddf = si->griddf;
	
	
  warning_flag = 0;
  Nstar = 0;	
  for( i = 0; i < NGRIDDF-1 ; i++ ){
    starCount[i] = 0;
    nstar = N*integral(dMstardE,griddf->E[i],griddf->E[i+1],si)/M;
    fs_f = (griddf->fEstar[i])/(griddf->fE[i]);
    if (fs_f > 1) warning_flag = 1;
    if (nstar >= 0.5){
      for(j = 0; starCount[i] < nstar && j < N; j++){		
	if (p[j].star_flag == 0){
	  if ( (griddf->E[i]) <= p[j].E && (p[j].E < griddf->E[i+1])){
	    p[j].star_flag = 1;
	    starCount[i]++;
	  }
	}
      }
    }
    Nstar += starCount[i];
  }
	
  if (warning_flag == 1){
    fprintf(stderr," WARNING: f(E)/fStars(E) was found to be less than 1 at one or more points of the grid. You are loosing some stellar mass! This can be solved with a smaller Mstar.\n"); 
  }
  si->Nstar = Nstar;
  fprintf(stderr,"Total number of star particles = %d\n", si->Nstar);
}


/*
** Routine for setting the remaining attributes
*/

void set_attributes(SI *si) {

  INT i, N,Nstar;
  DOUBLE mass,massStar;
  PARTICLE *p;
  PARTICLE *pstar;
    
  if(si->halo_flag == 1){
    N = si->N;
    p = si->p;
    mass = si->mass;
    for (i = 0; i < N; i++) {
      p[i].mass = mass;
    }
  }
	
  if(si->stars_flag == 1 && si->nostarpot_flag == 0){
    Nstar = si->Nstar;
    pstar = si->pstar;
    massStar = si->massStar;
    for (i = 0; i < Nstar; i++) {
      pstar[i].mass = massStar;
    }
  }
}


/*
** Routine for calculating center of mass position and velocity, and angular momentum
*/

void calculate_stuff( PARTICLE *bh, SI *si) {

  INT i, j, N,Nstar;
  DOUBLE mass, x, y, z, sumrvir;
  STUFF *stuff;
  PARTICLE *p;
  PARTICLE *pstar;
	
  sumrvir = 0;
  stuff = si->stuff;
  stuff->N = 0;
  stuff->Mp = 0;
  stuff->Ekin = 0;
  stuff->Epot = 0;
  for(i = 0; i < 4; i++) {
    stuff->Cr[i] = 0;
    stuff->Cv[i] = 0;
  }
  N = si->N;
  Nstar = si->Nstar; 
  p = si->p;
  pstar = si->pstar;
	
	
  if (bh->mass > 0) {
    stuff->N++;
    stuff->Mp += bh->mass;
  }
	
  stuff->N += N + Nstar;
	
  for (i = 0; i < N; i++) {
    if(p[i].r[0] < si->rimp) {
      si->rimp = p[i].r[0];
    }
    mass = p[i].mass;
    stuff->Mp += mass;
    stuff->Cr[1] += mass*p[i].r[1];
    stuff->Cr[2] += mass*p[i].r[2];
    stuff->Cr[3] += mass*p[i].r[3];
    stuff->Cv[1] += mass*p[i].v[1];
    stuff->Cv[2] += mass*p[i].v[2];
    stuff->Cv[3] += mass*p[i].v[3];
    stuff->Ekin += mass*p[i].v[0]*p[i].v[0]/2.0;
    stuff->Epot += mass*Pot(p[i].r[0],si);
    if (si->dorvirexact == 1) {
      for (j = i+1; j < N; j++) {
	x = p[i].r[1] - p[j].r[1];
	y = p[i].r[2] - p[j].r[2];
	z = p[i].r[3] - p[j].r[3];
	sumrvir += p[i].mass*p[j].mass/sqrt(x*x + y*y +z*z);
      }
    }
	
  }
    
	
  if(si->stars_flag == 1 && si->nostarpot_flag == 0){
		
    for (i = 0; i < Nstar; i++) {
      if(pstar[i].r[0] < si->rimp) {
	si->rimp = pstar[i].r[0];
      }
      mass = pstar[i].mass;
      stuff->Mp += mass;
      stuff->Cr[1] += mass*pstar[i].r[1];
      stuff->Cr[2] += mass*pstar[i].r[2];
      stuff->Cr[3] += mass*pstar[i].r[3];
      stuff->Cv[1] += mass*pstar[i].v[1];
      stuff->Cv[2] += mass*pstar[i].v[2];
      stuff->Cv[3] += mass*pstar[i].v[3];
      stuff->Ekin += mass*pstar[i].v[0]*pstar[i].v[0]/2.0;
      stuff->Epot += mass*Pot(pstar[i].r[0],si);
      if (si->dorvirexact == 1) {
	for (j = i+1; j < Nstar; j++) {
	  x = pstar[i].r[1] - pstar[j].r[1];
	  y = pstar[i].r[2] - pstar[j].r[2];
	  z = pstar[i].r[3] - pstar[j].r[3];
	  sumrvir += pstar[i].mass*pstar[j].mass/sqrt(x*x + y*y +z*z);
	}
      }
    }	
  }	
		
  for(i = 1; i < 4; i++) {
    stuff->Cr[i] = stuff->Cr[i]/stuff->Mp;
    stuff->Cv[i] = stuff->Cv[i]/stuff->Mp;
  }
  stuff->Cr[0] = sqrt(stuff->Cr[1]*stuff->Cr[1]+stuff->Cr[2]*stuff->Cr[2]+stuff->Cr[3]*stuff->Cr[3]);
  stuff->Cv[0] = sqrt(stuff->Cv[1]*stuff->Cv[1]+stuff->Cv[2]*stuff->Cv[2]+stuff->Cv[3]*stuff->Cv[3]);
  stuff->Epot = stuff->Epot/2.0;
  stuff->Etot = stuff->Ekin + stuff->Epot;
  if(si->halo_flag == 1){
    si->sp->rhalf = exp(lininterpolate(NGRIDR,si->gridr->MencHalo,si->gridr->logr,si->sp->M/2.0));
    si->r1 = pow(((3.0-si->sp->gamma)*si->mass)/(4*M_PI*si->sp->rho0*pow(si->sp->rs,si->sp->gamma)),1.0/(3.0-si->sp->gamma));
    si->r100 = pow(((3.0-si->sp->gamma)*100*si->mass)/(4*M_PI*si->sp->rho0*pow(si->sp->rs,si->sp->gamma)),1.0/(3.0-si->sp->gamma));
  }
  if(si->stars_flag == 1){
    si->sp->rhalfStar = exp(lininterpolate(NGRIDR,si->gridr->MencStar,si->gridr->logr,si->sp->Mstar/2.0));
  }
  if (si->dorvirexact == 1) {
    si->sp->rvir = (stuff->Mp*stuff->Mp)/(2.0*sumrvir);
  }
  else {
    sumrvir = (-1)*stuff->Epot/G;
    si->sp->rvir = (stuff->Mp*stuff->Mp)/(2.0*sumrvir);
  }
}


/*
** Routine for converting velocities from internal units to the sytems of units determined by "velConvert" **(definitions.h). Offsests are set here as well.
*/

void displace(SI *si){

  PARTICLE *p;
  PARTICLE *pstar;
  INT Nstar,N,i,j;

  p = si->p;
  pstar = si->pstar;
  N = si->N;
  Nstar = si->Nstar;

  for(i = 0; i < N; i++) {
    for(j = 0; j < 3 ; j++){
      p[i].r[j+1] = p[i].r[j+1] + si->deltapos[j];
      p[i].v[j+1] = p[i].v[j+1]*velConvert + si->deltavel[j];	
    }
  }

  for(i = 0; i < Nstar; i++) {
    for(j = 0; j < 3 ; j++){
      pstar[i].r[j+1] = pstar[i].r[j+1] + si->deltapos[j];
      pstar[i].v[j+1] = pstar[i].v[j+1]*velConvert + si->deltavel[j];	
    }
  }
}

/*
** Routine for transfering particles to tipsy structure
*/

void transfer_particles(const PARTICLE *bh,const SI *si, TIPSY_STRUCTURE *ts) {

  INT i, k, Ndark,Nstar,Ntotal;
  /* PARTICLE *p;*/

  Ndark = 0;
  if (bh->mass != 0) {
    Ndark++;
  }
  Ndark = Ndark + si->N; 
  Nstar =  si->Nstar;
  Ntotal = (si->nostarpot_flag == 0) ? Nstar + Ndark : Ndark;
  /*
  ** Initialise tipsy structure
  */
  ts->th = malloc(sizeof(TIPSY_HEADER));
  ts->gp = NULL;
  ts->dp = malloc(Ndark*sizeof(DARK_PARTICLE));
  ts->sp = malloc(Nstar*sizeof(STAR_PARTICLE));
  /*
  ** Initialise tipsy header
  */
  ts->th->time = 0;
  ts->th->ntotal = Ntotal;
  ts->th->ndim = 3;
  ts->th->ngas = 0;
  ts->th->ndark = Ndark;
  ts->th->nstar = Nstar;
  /*
  ** Transfer particles to tipsy structure
  */
    
      
  if (bh->mass != 0) {
    for (k = 0; k < 3; k++) {
      ts->dp->pos[k] = bh->r[k+1];
      ts->dp->vel[k] = bh->v[k+1];
    }
    ts->dp->mass = bh->mass;
    /**ts->dp->eps = bh->soft;
       ts->dp->phi = bh->Epot;*/ 
    ts->dp++;
  }
   

  for (i = 0; i < Ndark; i++, ts->dp++) {
    for (k = 0; k < 3; k++) {
      ts->dp->pos[k] = si->p[i].r[k+1];
      ts->dp->vel[k] = si->p[i].v[k+1];
    }
      
    ts->dp->mass = si->p[i].mass ;
    ts->dp->eps = si->soft0;
    ts->dp->phi = Pot(si->p[i].r[0],si);
  }
	
  for (i = 0; i < Nstar; i++, ts->sp++) {
    for (k = 0; k < 3; k++) {
      ts->sp->pos[k] = si->pstar[i].r[k+1];
      ts->sp->vel[k] = si->pstar[i].v[k+1];
    }
      
    ts->sp->mass = si->pstar[i].mass ;
    ts->sp->eps = si->soft0;
    ts->sp->phi = Pot(si->pstar[i].r[0],si);
  }	

  ts->dp = ts->dp - Ndark;
  ts->sp = ts->sp - Nstar;
    
}

/*
** Usage description
*/

void usage() {

  fprintf(stderr,"\n");
  fprintf(stderr,"You can specify the following arguments.\n\n");
  fprintf(stderr,"-halo               : generate a dark matter halo\n");
  fprintf(stderr,"-Nhalo <value>      : number of halo particles\n");
  fprintf(stderr,"-Mhalo <value>      : total mass of the halo (default Mhalo = 1)\n");
  fprintf(stderr,"-a <value>          : alpha parameter in the halo density profile\n");
  fprintf(stderr,"-b <value>          : beta parameter in the halo density profile\n");
  fprintf(stderr,"-c <value>          : gamma parameter in the halo density profile\n");
  fprintf(stderr,"-rs <value>         : scale radius (default rs = 1)\n");
  fprintf(stderr,"-rcutoff <value>    : cutoff radius for cutoff halo models (i.e. beta >= 3)\n");
  //    fprintf(stderr,"-soft0 <value>      : softening of particles (optional for writing on TIPSY file only)\n");
  fprintf(stderr,"-king               : generate a stellar component with a King profile\n");
  fprintf(stderr,"-hernquist          : generate a stellar component with a Hernquist profile\n");
  fprintf(stderr,"-plummer            : generate a stellar component with a Plummer profile\n");
  fprintf(stderr,"-Nstar <value>      : number of star particles\n");
  fprintf(stderr,"-Mstar <value>      : total stellar mass \n");
  fprintf(stderr,"-rc <value>         : king core radius\n");
  fprintf(stderr,"-rt <value>         : king tidal radius\n");
  fprintf(stderr,"-rhern <value>      : scale radius for Hernquist profile\n");
  fprintf(stderr,"-rp <value>         : scale radius for Plummer profile\n");
  fprintf(stderr,"-MBH <value>        : mass of black hole (default MBH = 0)\n");
  fprintf(stderr,"-name <value>       : name of the output file\n");
  fprintf(stderr,"-dx/dy/dz <value>   : position offset for the initial conditions (default All = 0)\n");
  fprintf(stderr,"-dvx/dvy/dvz <value>: velocity offset for the initial conditions (default All = 0)\n");
  fprintf(stderr,"-ogr                : set this flag for outputting grid in r in an ASCII file\n");
  fprintf(stderr,"-ogdf               : set this flag for outputting grid for distribution function in an ASCII file\n");
  fprintf(stderr,"-ogb                : set this flag for generating a GADGET2 initial conditions binary file\n");
  fprintf(stderr,"-otb                : set this flag for generating a TIPSY initial conditions binary file\n");
  fprintf(stderr,"-oift               : set this flag to write positions for IFRIT binary file \n");
  fprintf(stderr,"-opfs               : set this flag to write a table of density profiles in an ASCII file\n");
  fprintf(stderr,"-nostarpot          : set this flag for excluding the stellar potential \n");
  fprintf(stderr,"-randomseed <value> : set this flag for setting a value for a random seed (default: random value)\n");
  fprintf(stderr,"-dorvirexact        : set this flag for calculating rvir exactly via N^2 sum - Warning: time consuming for large N!\n");
  fprintf(stderr,"\n");
  exit(1);
}

/*
** This routine could be used for a fancy progress indicator. Not used at the moment.
*/

void rotate()
{
  static short index=0;
  fprintf(stderr," ");
  switch(index)
    {
    case 0 : fprintf(stderr,"\b\\");
      index++;
      break;
    case 1 : fprintf(stderr,"\b|");
      index++;
      break;
    case 2 : fprintf(stderr,"\b/");
      index++;
      break;
    case 3 : fprintf(stderr,"\b-");
      index = 0;
      break;
    }
}
