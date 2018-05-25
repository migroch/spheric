/* 
** io.c
**
** IO routines for spherIC_0.1
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "definitions.h"
#include "io.h"


/*
** Routine for writing grid in r to an ASCII file.
*/

void write_gridr(FILE *file,const SI *si) {

  INT i,halo_flag,stars_flag;
  GRIDR *gridr;
	
  gridr = si->gridr;
  halo_flag = si->halo_flag;
  stars_flag = si->stars_flag;
	
  if(halo_flag == 1 && stars_flag == 1){
    fprintf(file,"#     0             1             2             3             4             5             6             7\n");
    fprintf(file,"#     r          rhoHalo       rhoStar         rho         MencHalo      MencStar        Menc          Pot\n");
    for (i = 0; i < NGRIDR; i++) {
      fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
	      gridr->r[i],gridr->rhoHalo[i],gridr->rhoStar[i],gridr->rho[i],gridr->MencHalo[i],gridr->MencStar[i],gridr->Menc[i],gridr->Pot[i]);
    }
  }
	
  if(halo_flag == 1 && stars_flag == 0){
    fprintf(file,"#     0             1             2             3\n");
    fprintf(file,"#     r          rhoHalo       MencHalo        Pot\n");
    for (i = 0; i < NGRIDR; i++) {
      fprintf(file,OFD2" "OFD2" "OFD2" "OFD2"\n",
	      gridr->r[i],gridr->rhoHalo[i],gridr->MencHalo[i],gridr->Pot[i]);
    }
  }
	
  if(halo_flag == 0 && stars_flag == 1){
    fprintf(file,"#     0             1             2             3\n");
    fprintf(file,"#     r          rhoStar       MencStar        Pot\n");
    for (i = 0; i < NGRIDR; i++) {
      fprintf(file,OFD2" "OFD2" "OFD2" "OFD2"\n",
	      gridr->r[i],gridr->rhoStar[i],gridr->MencStar[i],gridr->Pot[i]);
    }
  }
	
}	

/*
** Routine for writing grid in f (distribution function) to an ASCII file.
*/

void write_griddf(FILE *file, const SI *si) {

  INT i,halo_flag,stars_flag;
	
  halo_flag = si->halo_flag;
  stars_flag = si->stars_flag;
  GRIDDF *griddf;
	
	
  if(halo_flag == 1 && stars_flag == 1){
    fprintf(file,"#     0             1             2             3             4\n");
    fprintf(file,"#     r             E             f           fstar          gE\n");
    griddf = si->griddf;
    for (i = 0; i < NGRIDDF; i++) {
      fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
	      griddf->r[i],griddf->E[i],griddf->fE[i],griddf->fEstar[i],griddf->gE[i]);
    }
  }

  if(halo_flag == 1 && stars_flag == 0){
    fprintf(file,"#     0             1             2             3\n");
    fprintf(file,"#     r             E           fhalo          gE\n");
    griddf = si->griddf;
    for (i = 0; i < NGRIDDF; i++) {
      fprintf(file,OFD2" "OFD2" "OFD2" "OFD2"\n",
	      griddf->r[i],griddf->E[i],griddf->fE[i],griddf->gE[i]);
    }
  }
	
  if(halo_flag == 0 && stars_flag == 1){
    fprintf(file,"#     0             1             2             3\n");
    fprintf(file,"#     r             E           fstar          gE\n");
    griddf = si->griddf;
    for (i = 0; i < NGRIDDF; i++) {
      fprintf(file,OFD2" "OFD2" "OFD2" "OFD2"\n",
	      griddf->r[i],griddf->E[i],griddf->fEstar[i],griddf->gE[i]);
    }
  }
}


/*
** Routine for writing ifrit binary with positions and star_flag.
*/

void write_ifrit(FILE *F,const SI *si){

  int dummy,n,i,Nstar,N;
  float xl, xh,temp; 
  xh = si->sp->rs*10;
  xl = -xh;
   
  N = si->N;
  Nstar = si->Nstar;
   
  if(si->nostarpot_flag == 1 || si->stars_flag == 0){   
    n = N;
    Nstar = 0;
  } else {
    n = Nstar + N;
  }
   
  dummy = sizeof(int);
  assert(fwrite(&dummy,sizeof(int),1,F) == 1);
  assert(fwrite(&n,sizeof(int),1,F) == 1);
  assert(fwrite(&dummy,sizeof(int),1,F) == 1);
   
  dummy = 6*sizeof(float);
  assert(fwrite(&dummy,sizeof(int),1,F) == 1);
  assert(fwrite(&xl,sizeof(float),1,F) == 1);
  assert(fwrite(&xl,sizeof(float),1,F) == 1);
  assert(fwrite(&xl,sizeof(float),1,F) == 1);
  assert(fwrite(&xh,sizeof(float),1,F) == 1);
  assert(fwrite(&xh,sizeof(float),1,F) == 1);
  assert(fwrite(&xh,sizeof(float),1,F) == 1);
  assert(fwrite(&dummy,sizeof(int),1,F) == 1);
   
  dummy = n*sizeof(float);
  assert(fwrite(&dummy,sizeof(int),1,F) == 1) ;
  for(i = 0; i < N; i++) {
    temp = si->p[i].r[1];
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }   
  for(i = 0; i < Nstar; i++) {
    temp = si->pstar[i].r[1];
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }
  assert(fwrite(&dummy,sizeof(int),1,F) == 1);
   
  assert(fwrite(&dummy,sizeof(int),1,F) == 1) ;
  for(i = 0; i < N; i++) {
    temp = si->p[i].r[2];
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }
  for(i = 0; i < Nstar; i++) {
    temp = si->pstar[i].r[2];
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }
  assert(fwrite(&dummy,sizeof(int),1,F) == 1);
   
  assert(fwrite(&dummy,sizeof(int),1,F) == 1) ;
  for(i = 0; i < N; i++) {
    temp = si->p[i].r[3];
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }
  for(i = 0; i < Nstar; i++) {
    temp = si->pstar[i].r[3];
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }
  assert(fwrite(&dummy,sizeof(int),1,F) == 1);
   
  assert(fwrite(&dummy,sizeof(int),1,F) == 1) ;
  for(i = 0; i < N; i++) {
    temp = si->p[i].r[0];
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }
  for(i = 0; i < Nstar; i++) {
    temp = si->pstar[i].r[0];
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }
  assert(fwrite(&dummy,sizeof(int),1,F) == 1);
   
  assert(fwrite(&dummy,sizeof(int),1,F) == 1) ;
  for(i = 0; i < N; i++) {
    temp = si->p[i].star_flag;
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }
  for(i = 0; i < Nstar; i++) {
    temp = si->pstar[i].star_flag;
    assert(fwrite(&temp,sizeof(float),1,F) == 1);
  }
  assert(fwrite(&dummy,sizeof(int),1,F) == 1);
   
}

/*
** Routine for writing density and velocity dispersion profiles for analysis.
*/

void write_profiles(FILE *F, const SI *si,int Ngrid, double rlow, double rup){

  int i,j;
  float dr;
  double rgrid[Ngrid],rho[Ngrid-1],rhoStar[Ngrid-1],vel2[Ngrid-1];
  double vel2star[Ngrid-1],M[Ngrid-1],Mstar[Ngrid-1],r,Msum,MstarSum,vr;
  PARTICLE *p;
  PARTICLE *pstar;

  p = si->p;
  pstar = si->pstar;
  dr = (rup - rlow)/(Ngrid-1.0);
  Msum = 0;
  MstarSum = 0;

  for (j = 0 ; j < Ngrid ; j++){
    if (j < Ngrid -1){
      rho[j] = 0 ;
      rhoStar[j] = 0;
      vel2[j] = 0;
      vel2star[j] = 0;
      M[j] = 0;
      Mstar[j] = 0;
    }
    rgrid[j] = rlow + j*dr;
  }

  if(si->nostarpot_flag == 1){
    for(i = 0 ; i < si->N ; i++){
      for (j = 0 ; j < Ngrid-1 ; j++){
	if ( (p[i].r[0] >= rgrid[j]) && (p[i].r[0] < rgrid[j+1]) ){ 
	  rho[j] += 1;
	  vr = (p[i].v[1]*p[i].r[1] + p[i].v[2]*p[i].r[2] + p[i].v[3]*p[i].r[3])/p[i].r[0];
	  vel2[j] += vr*vr;
	  if (p[i].star_flag > 0) rhoStar[j] += p[i].star_flag;
	}
      }
    }
  }

  if(si->nostarpot_flag == 0){

    for(i = 0 ; i < si->N ; i++){
      for (j = 0 ; j < Ngrid-1 ; j++){
	if ( (p[i].r[0] >= rgrid[j]) && (p[i].r[0] < rgrid[j+1]) ){ 
	  rho[j] += 1;
	  vr = (p[i].v[1]*p[i].r[1] + p[i].v[2]*p[i].r[2] + p[i].v[3]*p[i].r[3])/p[i].r[0];
	  vel2[j] +=  vr*vr;
	}
      }
    }
	
    for(i = 0 ; i < si->Nstar ; i++){
      for (j = 0 ; j < Ngrid-1 ; j++){
	if ( (pstar[i].r[0] >= rgrid[j]) && (pstar[i].r[0] < rgrid[j+1]) ){ 
	  rhoStar[j] += 1;
	  vr = (pstar[i].v[1]*pstar[i].r[1] + pstar[i].v[2]*pstar[i].r[2] + pstar[i].v[3]*pstar[i].r[3])/pstar[i].r[0];
	  vel2star[j] += vr*vr;
	}
      }
    }

  }

  fprintf(F,"#     0             1             2             3             4             5             6 \n");
  fprintf(F,"#     r          rhoHalo       rhoStar       MencHalo      MencStar       vel^2       velStar^2\n");
  for (j = 0 ; j < Ngrid - 1 ; j++){
    r = (rgrid[j] + rgrid[j+1])/2;
    vel2[j] = vel2[j]/rho[j];
    vel2star[j] = vel2star[j]/rhoStar[j];
    Msum += rho[j];
    MstarSum += rhoStar[j];
    M[j] = si->sp->M*Msum/si->N; 
    rho[j] = si->sp->M*rho[j]/(si->N*4*M_PI*r*r*dr);
    if(si->nostarpot_flag == 1){
      rhoStar[j] = si->sp->M*rhoStar[j]/(si->N*4*M_PI*r*r*dr);
      Mstar[j] = si->sp->M*MstarSum/si->N; 
    }else{
      rhoStar[j] = si->sp->Mstar*rhoStar[j]/(si->Nstar*4*M_PI*r*r*dr);
      Mstar[j] = si->sp->Mstar*MstarSum/si->Nstar; 
    }
    fprintf(F, OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n", r,rho[j],rhoStar[j],M[j],Mstar[j],vel2[j],vel2star[j]);
  }

}

/*
** Routine for writing Gadget binary with initial conditions.
*/


void write_gadget(FILE *fp, const PARTICLE *bh, const SI *si){

  int i, j, dummy, Ntotal,Ndark,Nstar, temp_i;
  float temp;
  PARTICLE *p;
  PARTICLE *pstar;
  GH gh;
    
  p = si->p;
  pstar = si->pstar; 
		
  /*
  **  Assign number of particles.
  */
   
  Nstar = si->Nstar;	
  Ndark = si->N; 
  Ntotal = (si->nostarpot_flag == 0) ? Nstar + Ndark : Ndark;	
  if (bh->mass != 0) {
    Ntotal++;
  }

  /*
  ** Initialise header
  */
  gh.npart[0] = 0;
  gh.npart[1] = Ndark;
  gh.npart[2] = 0;
  gh.npart[3] = 0;
  gh.npart[4] = Nstar;
  gh.npart[5] = 0;
  gh.mass[0] = 0;
  gh.mass[1] = si->mass;
  gh.mass[2] = 0;
  gh.mass[3] = 0;
  gh.mass[4] = si->massStar;
  gh.mass[5] = 0;
  gh.time = 0;
  gh.redshift = 0;
  gh.flag_sfr = 0;
  gh.flag_feedback = 0;
  gh.npartTotal[0] = 0;
  gh.npartTotal[1] = Ndark;
  gh.npartTotal[2] = 0;
  gh.npartTotal[3] = 0;
  gh.npartTotal[4] = Nstar;
  gh.npartTotal[5] = 0;
  gh.flag_cooling = 0;
  gh.num_files = 1;
  gh.BoxSize = 0;
  gh.Omega0 = 0;
  gh.OmegaLambda = 0;
  gh.HubbleParam = 0;
  gh.flag_stellarage = 0;
  gh.flag_metals = 0;
  for (i = 0; i < 6; i++) {
    gh.npartTotalHighWord[i] = 0;
  }
  gh.flag_entropy_instead_u = 0;
  /*
  ** Write out header
  */
  dummy = sizeof(GH);
  assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
  assert(fwrite(&gh,sizeof(GH),1,fp) == 1);
  assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
  /*
  ** Write out positions
  */
  dummy = 3*Ntotal*sizeof(float);
  assert(fwrite(&dummy,sizeof(int),1,fp) == 1); 
    
  for(i = 0; i < Ndark; i++) {
    if (p[i].star_flag == 0) {
      for(j = 1; j < 4; j++) {
	temp = p[i].r[j];
	assert(fwrite(&temp,sizeof(float),1,fp) == 1);
      }
    }
  }
    
  // 	if(si->nostarpot_flag == 1){
  // 		for(i = 0; i < Ntotal; i++) {
  // 		if (p[i].star_flag == 1) {
  // 			for(j = 1; j < 4; j++) {
  // 				temp = p[i].r[j];
  // 				assert(fwrite(&temp,sizeof(float),1,fp) == 1);
  // 				}
  // 			}
  // 		}
  // 	}else{
  for(i = 0; i < Nstar; i++) {
    for(j = 1; j < 4; j++) {
      temp = pstar[i].r[j];
      assert(fwrite(&temp,sizeof(float),1,fp) == 1);
    }
		
  }
  //	}		
	
  assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    
  /*
  ** Write out velocities
  */
	
  dummy = 3*Ntotal*sizeof(float);
  assert(fwrite(&dummy,sizeof(int),1,fp) == 1); 
    		
  for(i = 0; i < Ndark; i++) {
    if (p[i].star_flag == 0) {
      for(j = 1; j < 4; j++) {
	temp = p[i].v[j];
	assert(fwrite(&temp,sizeof(float),1,fp) == 1);
      }
    }
  }
	
  //     if(si->nostarpot_flag == 1){
  // 		for(i = 0; i < Ntotal; i++) {
  // 		if (p[i].star_flag == 1) {
  // 			for(j = 1; j < 4; j++) {
  // 				temp = p[i].v[j];
  // 				assert(fwrite(&temp,sizeof(float),1,fp) == 1);
  // 				}
  // 			}
  // 		}
  // 	}else{
  for(i = 0; i < Nstar; i++) {
    for(j = 1; j < 4; j++) {
      temp = pstar[i].v[j];
      assert(fwrite(&temp,sizeof(float),1,fp) == 1);
    }		
  }
  //	}
			
  assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    
  /*
  ** Write out ids
  */
  dummy = Ntotal*sizeof(int);
  assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    
  for(i = 0; i < Ndark; i++) {
    temp_i = p[i].index;
    assert(fwrite(&temp_i,sizeof(int),1,fp) == 1);
  }
  for(i = 0; i < Nstar; i++) {
    temp_i = pstar[i].index;
    assert(fwrite(&temp_i,sizeof(int),1,fp) == 1);
  }
  assert(fwrite(&dummy,sizeof(int),1,fp) == 1);

}

/*
** Routine for writing ACII files with initial conditions.
*/

void write_ics(FILE *F, const SI *si,int flag){

  int i,j;

  //     if (bh->mass != 0) {
  // 	assert(fprintf(file,OFI1" ",index) > 0);
  //         assert(fprintf(file,OFD5" ",bh->mass) > 0);
  // 	for (j = 0; j < 3; j++) {
  // 	    assert(fprintf(file,OFD6" ",bh->r[j+1]) > 0);
  //             }
  //         for (j = 0; j < 3; j++) {
  //             assert(fprintf(file,OFD6" ",bh->v[j+1]) > 0);
  //             }
  // 	assert(fprintf(file,"\n") > 0);
  // 	index++;
  //         }
	
  if (flag == 1) {
    fprintf(F,"#0       1              2              3              4              5              6              7              8\n");
    fprintf(F,"#n   haloPmass        haloX          haloY          haloZ          haloVx         haloVy        haloVz      f(E)star/f(E)\n");
    for (i = 0; i < si->N; i++) {
      assert(fprintf(F,OFI1" ",si->p[i].index) > 0);
      assert(fprintf(F,OFD3" ",si->p[i].mass) > 0);
      for (j = 0; j < 3; j++) {
	assert(fprintf(F,OFD3" ",si->p[i].r[j+1]) > 0);
      }
      for (j = 0; j < 3; j++) {
	assert(fprintf(F,OFD3" ",si->p[i].v[j+1]) > 0);
      }
      assert(fprintf(F,OFD3" ",si->p[i].star_flag) > 0);
      assert(fprintf(F,"\n") > 0);
    }
  }
	
  if (flag == 2) {
    fprintf(F,"#0       1             2              3              4              5              6              7\n");
    fprintf(F,"#n   starPmass       starX          starY          starZ          starVx         starVy         starVz\n");		
    for (i = 0; i < si->Nstar; i++) {
      assert(fprintf(F,OFI1" ",si->pstar[i].index) > 0);
      assert(fprintf(F,OFD3" ",si->pstar[i].mass) > 0);
      for (j = 0; j < 3; j++) {
	assert(fprintf(F,OFD3" ",si->pstar[i].r[j+1]) > 0);
      }
      for (j = 0; j < 3; j++) {
	assert(fprintf(F,OFD3" ",si->pstar[i].v[j+1]) > 0);
      }	
      assert(fprintf(F,"\n") > 0);
		
    }
  }
}

/*
** Routine for writing typsy binary with initial conditions.
*/

void write_tipsy(FILE *fp, const TIPSY_STRUCTURE *ts) {

  TIPSY_HEADER *th;
  GAS_PARTICLE *gp;
  DARK_PARTICLE *dp;
  STAR_PARTICLE *sp;

  th = ts->th;
  gp = ts->gp;
  dp = ts->dp;
  sp = ts->sp;
  /*
  ** Write out header
  */
  assert(fwrite(th,sizeof(TIPSY_HEADER),1,fp) == 1);
  /*
  ** Write out particles
  */
  assert(fwrite(gp,sizeof(GAS_PARTICLE),th->ngas,fp) == th->ngas);
  assert(fwrite(dp,sizeof(DARK_PARTICLE),th->ndark,fp) == th->ndark);
  assert(fwrite(sp,sizeof(STAR_PARTICLE),th->nstar,fp) == th->nstar);
}
