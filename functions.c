/* 
** functions.c
**
** Functions for spherIC_0.1
*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "definitions.h"
#include "functions.h"

/*
** Functions for integration via Simpson's rule
** Iteration until precision (TOL) is reached 
*/

DOUBLE integral(DOUBLE (*function)(DOUBLE,const SI (*)), DOUBLE a, DOUBLE b, const SI *si) {

  INT i;
  DOUBLE T[3], SA, SB;
  T[0] = trapez(function,a,b,0,si);
  T[1] = trapez(function,a,b,1,si);
  T[2] = trapez(function,a,b,2,si);
  SA = (4.0*T[1]-T[0])/3.0;
  SB = (4.0*T[2]-T[1])/3.0;
  for(i = 3; i < (NINTMAX+1); i++) {
    T[0] = T[1];
    T[1] = T[2];
    T[2] = trapez(function,a,b,i,si);
    SA = SB;
    SB = (4.0*T[2]-T[1])/3.0;
    if(i > NINTMIN) {
      if((fabs(SB-SA) < TOL*fabs(SB)) || (SA == 0.0 && SB == 0.0)) {
	return SB;
      }
    }
  }
  fprintf(stderr,"Warning!\n");
  fprintf(stderr,"Too many steps in function integral!\n");
  fprintf(stderr,"a = "OFD3"\n",a);
  fprintf(stderr,"b = "OFD3"\n",b);
  fprintf(stderr,"sum = "OFD3"\n",SB);
  fprintf(stderr,"Sum not converged within tolerance of %e\n",TOL);
  fprintf(stderr,"Abort tolerance was %e\n",fabs((SA-SB)/SB));
  return SB;
}

DOUBLE trapez(DOUBLE (*function)(DOUBLE, const SI(*)), DOUBLE a, DOUBLE b, INT n, const SI *si) {

  INT i, j, N[n+1];
  DOUBLE deltax, sum, sumN;

  if (n == 0) {
    return (0.5*(b-a)*((*function)(a,si)+(*function)(b,si)));
  }
  else {
    for (i = 0; i < (n+1); i++) {
      N[i] = pow(2,i);
    }
    sum = 0.5*((*function)(a,si)+(*function)(b,si));
    for(i = 1; i < (n+1); i++) {
      deltax = (b-a)/N[i];
      sumN = 0;
      for (j = 1; j < (N[i-1]+1); j++) {
	sumN += (*function)(a+(2*j-1)*deltax,si);
      }
      sum = sum + sumN;
    }
    return (sum*(b-a)/N[n]);
  }
}

DOUBLE integraldf(DOUBLE (*d2anyrhodPhi2)(DOUBLE,const SI (*)),INT j, const SI *si) {
    
  INT i;
  DOUBLE T[3], SA, SB;

  T[0] = trapezdf(d2anyrhodPhi2,j,0,si);
  T[1] = trapezdf(d2anyrhodPhi2,j,1,si);
  T[2] = trapezdf(d2anyrhodPhi2,j,2,si);
  SA = (4.0*T[1]-T[0])/3.0;
  SB = (4.0*T[2]-T[1])/3.0;
  for(i = 3; i < (NINTMAXDF+1); i++) {
    T[0] = T[1];
    T[1] = T[2];
    T[2] = trapezdf(d2anyrhodPhi2,j,i,si);
    SA = SB;
    SB = (4.0*T[2]-T[1])/3.0;
    if(i > NINTMINDF) {
      if((fabs(SB-SA) < TOLDF*fabs(SB)) || (SA == 0.0 && SB == 0.0)) {
	return SB;
      }
    }
  }
  fprintf(stderr,"Warning!\n");
  fprintf(stderr,"Too many steps in function integraldf!\n");
  fprintf(stderr,"E = "OFD3" LU^2 Gyr^-2\n",si->gridr->Pot[j]);
  fprintf(stderr,"r = "OFD3" LU\n",si->gridr->r[j]);
  fprintf(stderr,"f(E) = "OFD3" MU Gyr^3 LU^-6\n",SB);
  fprintf(stderr,"f(E) not converged within tolerance of %e\n",TOLDF);
  fprintf(stderr,"Abort tolerance was %e\n",fabs((SA-SB)/SB));
  if (SB != SB) {
    fprintf(stderr,"f(E) value set to %e MU Gyr^3 LU^-6\n",DFFAILUREMAX);
    return DFFAILUREMAX;
  }
  else {
    return SB;
  }
}

DOUBLE trapezdf(DOUBLE (*d2anyrhodPhi2)(DOUBLE,const SI (*)),INT k, INT n,const SI *si) {

    
  INT i, j, N[n+1];
  DOUBLE deltax, sum, sumN;
  DOUBLE x, xlower, xupper;
  DOUBLE r, rlower, rupper;
  DOUBLE Phi, Philower, Phiupper;
  GRIDR *gridr;

  gridr = si->gridr;
  Philower = gridr->Pot[NGRIDR-1];
  Phiupper = gridr->Pot[k];
  rupper = gridr->r[k];
  xlower = asin(sqrt(Philower/Phiupper));
  xupper = M_PI/2;  
  rlower = gridr->r[NGRIDR-1];
  if (n == 0) {
    return ((sqrt(-Phiupper/2)/(M_PI*M_PI))*0.5*(xupper-xlower)*((*d2anyrhodPhi2)(rlower,si)*sin(xlower)+(*d2anyrhodPhi2)(rupper,si)*sin(xupper)));
  }
  else {
    for (i = 0; i < (n+1); i++) {
      N[i] = pow(2,i);
    }
    sum = 0.5*((*d2anyrhodPhi2)(rlower,si)*sin(xlower)+(*d2anyrhodPhi2)(rupper,si)*sin(xupper));
    for(i = 1; i < (n+1); i++) {
      deltax = (xupper-xlower)/N[i];
      sumN = 0;
      for (j = 1; j < (N[i-1]+1); j++) {
	x = xlower+(2*j-1)*deltax;
	Phi =  Phiupper*sin(x)*sin(x);
	r = exp(lininterpolate(NGRIDR,gridr->logPot,gridr->logr,log(-Phi))); 
	sumN += (*d2anyrhodPhi2)(r,si)*sin(x);
      }
      sum = sum + sumN;
    }
    return ((sqrt(-Phiupper/2)/(M_PI*M_PI))*sum*(xupper-xlower)/N[n]);
  }
}

/*
** Function for locating x on a monotonic increasing or decreasing grid 
*/

INT locate(INT n, const DOUBLE *grid, DOUBLE x) {

  INT jl, jm, ju, ascend;

  jl = -1;
  ju = n;
  ascend = (grid[n-1] >= grid[0]);
  while ((ju-jl) > 1) {
    jm = (ju+jl)/2;
    if ((x >= grid[jm]) == ascend) {
      jl = jm;
    }
    else {
      ju = jm;
    }
  }
  if (x == grid[0]) {
    return 0;
  }
  else if (x == grid[n-1]) {
    return (n-2);
  }
  else {
    return jl;
  }
}

/*
** Function for linear interpolation between two grid points 
*/

DOUBLE lininterpolate(INT n, const DOUBLE *gridx, const DOUBLE *gridy, DOUBLE x) {

  INT i; 
  DOUBLE dgx, dgy, dx, m, y, RE;

  i = locate(n,gridx,x);
  if (i < 0) {
    fprintf(stderr,"Warning!\n");
    fprintf(stderr,"x = "OFD3" was below range of array! Index i = "OFI1" / n = "OFI1".\n",x,i,n);
    fprintf(stderr,"Array between "OFD3" and "OFD3".\n",gridx[0],gridx[n-1]);
    RE = fabs((x-gridx[0])/gridx[0]);
    if (RE < TOLLININT) {
      fprintf(stderr,"Relative error (= "OFD1") at lower boundary was within tolerance of %e\n",RE,TOLLININT);
      fprintf(stderr,"Index set to i = 0\n");
      i = 0;
    }
  }
  else if (i > n-2) {
    RE = fabs((x-gridx[n-1])/gridx[n-1]);
    if (RE > TOLLININT) {
      fprintf(stderr,"ERROR!\n");
      fprintf(stderr,"x = "OFD3" was above range of array! Index i = "OFI1" / n = "OFI1".\n",x,i,n);
      fprintf(stderr,"Array between "OFD3" and "OFD3".\n",gridx[0],gridx[n-1]);		
      fprintf(stderr,"Relative error (= "OFD1") at upper boundary was greater than the tolerance value of %e\n",RE,TOLLININT);
    }
    else i = n-2;
  }
  dgy = gridy[i+1] - gridy[i];
  dgx = gridx[i+1] - gridx[i];
  m = dgy/dgx;
  dx = x - gridx[i];
  y = gridy[i] + dx*m;

  return y;
}

/*
** Function for generating random numbers between [0,1) 
*/

DOUBLE rand01() {

  return ( ((DOUBLE) rand()) / ((DOUBLE) (RAND_MAX + 1.0)) );
}

/* 
** alpha-beta-gamma density function with exponential cutoff 
** except for finite mass models 
*/

DOUBLE rho(DOUBLE r, const SI *si) {

  DOUBLE fac1, fac2;
  SP *sp;

  sp = si->sp;
  if (sp->beta > 3) {
    /*
    ** Finite mass models
    */
    return (sp->rho0/tau(r,si));
  }
  else {
    /*
    ** Cutoff models
    */
    if (r <= sp->rcutoff) {
      return (sp->rho0/tau(r,si));
    }
    else {
      fac1 = pow((r/sp->rcutoff),sp->delta);
      fac2 = exp(-(r-sp->rcutoff)/sp->rdecay);
      return (sp->rho0/tau(sp->rcutoff,si)*fac1*fac2);
    }
  }
}

/* 
** Derivative of density drho/dr 
*/

DOUBLE drhodr(DOUBLE r, const SI *si) {

  SP *sp;

  sp = si->sp;
  if (sp->beta > 3) {
    /*
    ** Finite mass models
    */
    return (-rho(r,si)*eta(r,si));
  }
  else {
    /*
    ** Cutoff models
    */
    if (r <= sp->rcutoff) {
      return (-rho(r,si)*eta(r,si));
    }
    else {
      return (rho(r,si)*(sp->delta/r-1/sp->rdecay));
    }
  }
}

/* 
** Derivative of density d^2rho/dr^2 
*/

DOUBLE d2rhodr2(DOUBLE r, const SI *si) {

  SP *sp;

  sp = si->sp;
  if (sp->beta > 3) {
    /*
    ** Finite mass models
    */
    return (rho(r,si)*(pow(eta(r,si),2)-detadr(r,si)));
  }
  else {
    /*
    ** Cutoff models
    */
    if (r <= sp->rcutoff) {
      return (rho(r,si)*(pow(eta(r,si),2)-detadr(r,si)));
    }
    else {
      return (rho(r,si)*(pow((sp->delta/r-1/sp->rdecay),2)-sp->delta/(r*r)));
    }
  }
}

/* 
** Derivative of density d^2rho/dPhi^2 
*/

DOUBLE d2rhodPhi2(DOUBLE r, const SI *si) {

  DOUBLE Mencr,rhor;
  DOUBLE fac1, fac2;
  SP *sp;
  sp = si->sp;
	
  if(si->nostarpot_flag == 1 || si->king_flag == 0){
    Mencr = MencHalo(r,si);
    rhor = rho(r,si);
  }	
  else{
    Mencr = Menc(r,si);
    rhor = rho(r,si) + si->sp->K*rhoStar(r,si);
  }
		
		
  fac1 = r*r/(G*G*Mencr*Mencr);
  fac2 = 2*r-4*M_PI*r*r*r*r*rhor/Mencr;
  return (fac1*(r*r*d2rhodr2(r,si)+fac2*drhodr(r,si)));
}

/* 
** Derivative of density dlrho/dlr 
*/

DOUBLE dlrhodlr(DOUBLE r, const SI *si) {

  return (r/rho(r,si)*drhodr(r,si));
}

/*
** Auxiliary functions 
*/

DOUBLE eta(DOUBLE r, const SI *si) {

  DOUBLE fac1, fac2, fac3;
  SP *sp;

  sp = si->sp;
  fac1 = (sp->beta-sp->gamma)/sp->rs;
  fac2 = pow((r/sp->rs),(sp->alpha-1));
  fac3 = 1+pow((r/sp->rs),sp->alpha);
  return ((sp->gamma/r)+(fac1*fac2/fac3));
}

DOUBLE detadr(DOUBLE r, const SI *si) {

  DOUBLE fac1, fac2, fac3, fac4;
  SP *sp;

  sp = si->sp;
  fac1 = (sp->beta-sp->gamma)/(sp->rs*sp->rs);
  fac2 = pow((r/sp->rs),(sp->alpha-2));
  fac3 = (sp->alpha-1)*(1+pow((r/sp->rs),sp->alpha))-sp->alpha*pow((r/sp->rs),sp->alpha);
  fac4 = pow((1+pow((r/sp->rs),sp->alpha)),2);
  return (-sp->gamma/(r*r)+fac1*fac2*fac3/fac4);
}

DOUBLE tau(DOUBLE r, const SI *si) {
    
  DOUBLE exp1, exp2, exp3;
  DOUBLE fac1, fac2, fac3;
  SP *sp;

  sp = si->sp;
  exp1 = sp->gamma;
  exp2 = sp->alpha;
  exp3 = (sp->beta-sp->gamma)/sp->alpha;
  fac1 = pow(r/sp->rs,exp1);
  fac2 = 1+pow(r/sp->rs,exp2);
  fac3 = pow(fac2,exp3);
  return (fac1*fac3);
}

/*
** Integrand of IM integral 
*/

DOUBLE integrandIM(DOUBLE r, const SI *si) {
    
  DOUBLE exp1, exp2, exp3;
  DOUBLE fac1, fac2, fac3;
  SP *sp;

  sp = si->sp;
  exp1 = 2-sp->gamma;
  exp2 = sp->alpha;
  exp3 = (sp->beta-sp->gamma)/sp->alpha;
  fac1 = pow(r,exp1);
  fac2 = 1+pow(r,exp2);
  fac3 = pow(fac2,exp3);
  return (fac1/fac3);
}

/*
** Integrand of IMcutoff integral
*/

DOUBLE integrandIMcutoff(DOUBLE r, const SI *si) {

  DOUBLE fac1, fac2, fac3;
  SP *sp;

  sp = si->sp;
  fac1 = r*r;
  fac2 = pow(r/sp->rcutoff,sp->delta);
  fac3 = exp(-(r-sp->rcutoff)/sp->rdecay);
  return (fac1*fac2*fac3);
}

/*
** Integrand of enclosed mass integral 
*/

DOUBLE integrandMenc(DOUBLE r, const SI *si) {

  return (4*M_PI*r*r*(rho(r,si)));
}

/*
** Integrand of outer potential integral 
*/

DOUBLE integrandPotHalo(DOUBLE r,const SI *si) {
	
  return (4*M_PI*r*rho(r,si));
}

DOUBLE integrandPotStar(DOUBLE r,const SI *si) {
	
  return (4*M_PI*r*si->sp->K*rhoStar(r,si));
}




/* 
** Function for calculating total enclosed mass within radius r 
*/

//Total Mass
DOUBLE Menc(DOUBLE r, const SI *si) {

  return (exp(lininterpolate(NGRIDR,si->gridr->logr,si->gridr->logMenc,log(r))));
}
//Just Halo
DOUBLE MencHalo(DOUBLE r, const SI *si) {

  return (exp(lininterpolate(NGRIDR,si->gridr->logr,si->gridr->logMencHalo,log(r))));
}
//Just Stars
DOUBLE MencStar(DOUBLE r, const SI *si) {

  return (exp(lininterpolate(NGRIDR,si->gridr->logr,si->gridr->logMencStar,log(r))));
}

/*
** Function for calculating potential at radius r 
*/

DOUBLE Pot(DOUBLE r, const SI *si) {

  return (-exp(lininterpolate(NGRIDR,si->gridr->logr,si->gridr->logPot,log(r))));
}

/*
** Function for calculating the escape velocity at position r 
*/

DOUBLE vescape(DOUBLE r, const SI *si) {

  return (sqrt(2.0*fabs(Pot(r,si)-si->gridr->Pot[NGRIDR-1])));
}

DOUBLE vescapeStar(DOUBLE r,const SI *si) {

  return (sqrt(2.0*fabs(Pot(r,si)-Pot(si->routerStar,si))));
}

/*
** Function for calculating the dynamical time at position r 
*/

DOUBLE Tdyn(DOUBLE r, const SI *si) {

  return (2*M_PI*sqrt((r*r*r)/(G*Menc(r,si))));
}

/*
** Function for calculating the value of the distribution function at energy E 
*/

DOUBLE f1(DOUBLE E, const SI *si) {

  return (lininterpolate(NGRIDDF,si->griddf->logE,si->griddf->fE,log(-E)));
}

DOUBLE f1star(DOUBLE E, const SI *si) {

  return (lininterpolate(NGRIDDF,si->griddf->logE,si->griddf->fEstar,log(-E)));
}

/*
** Function for calculating the value of the distribution function at position r 
*/

DOUBLE f2(DOUBLE r, const SI *si) {

  return (lininterpolate(NGRIDDF,si->griddf->logr,si->griddf->fE,log(r)));
}

DOUBLE f2star(DOUBLE r, const SI *si) {

  return (lininterpolate(NGRIDDF,si->griddf->logr,si->griddf->fEstar,log(r)));
}


/*
** Function for calculating the value of dM/dE 
*/

DOUBLE dMdE(DOUBLE E,const SI *si){

  return (f1(E,si)*gEi(E,si));
}

/*
** Function for calculating the value of the dMstar/dE 
*/

DOUBLE dMstardE(DOUBLE E,const SI *si){

  return (f1star(E,si)*gEi(E,si));
}

/////////Some functions for calculating the value of the differential energy distribution g(E).///////////

/*
** Function for calculating the value of differential energy distribution g(E) at energy E by interpolation. 
*/

DOUBLE gEi(DOUBLE E, const SI *si) {
	
  return (lininterpolate(NGRIDDF,si->griddf->logE,si->griddf->gE,log(-E)));
}


/*
** function for calculating g(E).
*/

DOUBLE gE(DOUBLE E,DOUBLE re,const SI *si){

  while ((E-Pot(re,si)) < 0){
    re *= 0.99999;
  }
  return (integralgE(integrandgE,E,si->rinner,re,si));
}


/*
** Integrand for the differential energy distribution g(E).
*/

DOUBLE integrandgE(DOUBLE r,DOUBLE E,const SI *si){ 

  DOUBLE fac1,fac2;

  fac1 = 4*M_PI*4*M_PI*r*r;
  fac2 = sqrt(2*(E-Pot(r,si)));

  return (fac1*fac2);
}

/*
** Integral for g(E)
*/

DOUBLE integralgE(DOUBLE (*function)(DOUBLE,DOUBLE,const SI (*)),DOUBLE E, DOUBLE a, DOUBLE b, const SI *si) {

  INT i;
  DOUBLE T[3], SA, SB;

  T[0] = trapezgE(function,E,a,b,0,si);
  T[1] = trapezgE(function,E,a,b,1,si);
  T[2] = trapezgE(function,E,a,b,2,si);
  SA = (4.0*T[1]-T[0])/3.0;
  SB = (4.0*T[2]-T[1])/3.0;
  for(i = 3; i < (NINTMAX+1); i++) {
    T[0] = T[1];
    T[1] = T[2];
    T[2] = trapezgE(function,E,a,b,i,si);
    SA = SB;
    SB = (4.0*T[2]-T[1])/3.0;
    if(i > NINTMINDF) {
      if((fabs(SB-SA) < TOLDF*fabs(SB)) || (SA == 0.0 && SB == 0.0)) {
	return SB;
      }
    }
  }
	
  fprintf(stderr,"Warning!\n");
  fprintf(stderr,"Too many steps in function integralgE!\n");
  fprintf(stderr,"a = "OFD3"\n",a);
  fprintf(stderr,"b = "OFD3"\n",b);
  fprintf(stderr,"sum = "OFD3"\n",SB);
  fprintf(stderr,"Sum not converged within tolerance of %e\n",TOL);
  fprintf(stderr,"Abort tolerance was %e\n",fabs((SA-SB)/SB));
  return SB;
}

DOUBLE trapezgE(DOUBLE (*function)(DOUBLE, DOUBLE,const SI (*) ),DOUBLE E ,DOUBLE a, DOUBLE b, INT n,const SI *si) {

  INT i, j, N[n+1];
  DOUBLE deltax, sum, sumN;

  if (n == 0) {
    return (0.5*(b-a)*((*function)(a,E,si)+(*function)(b,E,si)));
  }
  else {
    for (i = 0; i < (n+1); i++) {
      N[i] = pow(2,i);
    }
    sum = 0.5*((*function)(a,E,si)+(*function)(b,E,si));
    for(i = 1; i < (n+1); i++) {
      deltax = (b-a)/N[i];
      sumN = 0;
      for (j = 1; j < (N[i-1]+1); j++) {
	sumN += (*function)(a+(2*j-1)*deltax,E,si);
      }
      sum = sum + sumN;
    }
    return (sum*(b-a)/N[n]);
  }
}


///////////////////////////functions for the Stellar models////////////////////

/*
** King density profile
*/

DOUBLE xKing (DOUBLE r, const SI *si){


  DOUBLE rc,rt;
  DOUBLE fac1,fac2;
  SP *sp;

  sp = si->sp;
  rc = sp->rc;
  rt = sp->rt;

  fac1 = 1 + (r/rc)*(r/rc);
  fac2 = 1 + (rt/rc)*(rt/rc);
  return (sqrt(fac1/fac2));

}


DOUBLE rhoKing(DOUBLE r, const SI *si) {

  DOUBLE rc,rt;
  SP *sp;

  sp = si->sp;
  rc=sp->rc;
  rt=sp->rt;	
	
  if(r >= rt ){
    return ( 0.0 );
  }
  else {
	
    DOUBLE fac1, fac2;
    fac1 = acos(xKing(r,si))/xKing(r,si);
    fac2 = sqrt(1 - xKing(r,si)*xKing(r,si));
    return ((fac1-fac2)/(xKing(r,si)*xKing(r,si)));	
  }
}

/*
**Approximate mass enclosed at small r for the King Model
*/


DOUBLE MencInnerKing(DOUBLE r, const SI *si){

  DOUBLE rc,rt,ck;
  DOUBLE fac1,fac2,fac3,fac4,fac5;
  SP *sp;

  sp = si->sp;
  rc=sp->rc;
  rt=sp->rt;
  ck=rt/rc;
			
  fac1 = ck/(1+ck*ck);
  fac2 = sqrt(1/(1+ck*ck));
  fac3 = acos(fac2);
  fac4 = (fac3 - fac1)/(fac2*fac2*fac2);
  fac5 = ck*ck*(fac2 - (3.0*fac3/(fac2*ck)))/(2.0*fac1*rc*rc);
	
	
  return (4*M_PI*(fac4*r*r*r/3.0 + fac5*r*r*r*r*r/5.0));
}

/*
**Derivative with respect to r for the x factor in the King density profile
*/

DOUBLE dxKingdr(DOUBLE r, const SI *si) {

  DOUBLE rc,rt,ck;
  SP *sp;

  sp = si->sp;
  rc = sp->rc;
  rt = sp->rt;
  ck = rt/rc;

  return (r/((1 + ck*ck)*xKing(r,si)*rc*rc));
}

/*
**Second derivative with respect to r for the x factor in the King density profile
*/

DOUBLE d2xKingdr2(DOUBLE r, const SI *si) {

  DOUBLE rc,rt,ck;
  DOUBLE fac1,fac2;
  SP *sp;

  sp = si->sp;
  rc = sp->rc;
  rt = sp->rt;
  ck = rt/rc;

  fac1 =  1/((1 + ck*ck)*xKing(r,si)*rc*rc);
  fac2 =  (r*r)/((1 + ck*ck)*(1 + ck*ck)*pow(xKing(r,si),3)*pow(rc,4));

  return (fac1-fac2);
}


/*
**Derivative with respect to x for the King density profile
*/

DOUBLE drhoKingdxKing(DOUBLE r, const SI *si) {

  DOUBLE fac1,fac2;

  fac1 = 3*acos(xKing(r,si));
  fac2 = xKing(r,si)*sqrt(1-xKing(r,si)*xKing(r,si));
	

  return ((-fac1+fac2)/pow(xKing(r,si),4));
}

/*
**Second derivative with respect to x for the King density profile
*/

DOUBLE d2rhoKingdxKing2(DOUBLE r, const SI *si) {

  DOUBLE fac1,fac2,fac3;

  fac1 = 6*acos(xKing(r,si));
  fac2 = pow(xKing(r,si),3)/sqrt(1-xKing(r,si)*xKing(r,si));
  fac3 = 2.0/pow(xKing(r,si),5);

  return ((fac1+fac2)*fac3);
}

/*
**Derivative with respect to r for the King density profile
*/

DOUBLE drhoKingdr(DOUBLE r, const SI *si){
	
  SP *sp;
  DOUBLE rt;
	
  sp = si->sp;
  rt=sp->rt;	
	
  if(r >= rt ){
    return ( 0.0 );
  }else{
    return (drhoKingdxKing(r,si)*dxKingdr(r,si));
  }
}


/*
**Second derivative with respect to x for the King density profile
*/

DOUBLE d2rhoKingdr2(DOUBLE r, const SI *si){
	
  SP *sp;
  DOUBLE rt;
	
  sp = si->sp;
  rt=sp->rt;	
	
  if(r >= rt ){
    return ( 0.0 );
  }else{
    return (d2rhoKingdxKing2(r,si)*dxKingdr(r,si)*dxKingdr(r,si) + drhoKingdxKing(r,si)*d2xKingdr2(r,si));
  }
}

/*
** Plummer density profile
*/


DOUBLE rhoPlummer(DOUBLE r, const SI *si){

  DOUBLE rp;
  DOUBLE fac1;
  SP *sp;


  sp = si->sp;
  rp=sp->rp;

  fac1 = 1 + (r*r)/(rp*rp);

  return (1/pow(fac1,5.0/2.0));

}

/*
**Approximate mass enclosed at small r for the Plummer Model
*/

DOUBLE MencInnerPlummer(DOUBLE r, const SI *si){

  DOUBLE rp;
  DOUBLE fac1,fac2;
  SP *sp;

  sp = si->sp;
  rp=sp->rp;

  fac1 = r*r*r/3.0;
  fac2 = r*r*r*r*r/(2.0*rp*rp);

  return ( 4*M_PI*(fac1-fac2));
}


/*
** Derivative with respect to r for the Plummer density profile
*/


DOUBLE drhoPlummerdr(DOUBLE r, const SI *si){

  DOUBLE rp;
  DOUBLE fac1,fac2;
  SP *sp;


  sp = si->sp;
  rp=sp->rp;

  fac1 = 1 + (r*r)/(rp*rp);
  fac2 = (5.0*r)/(rp*rp);
	
  return (-fac2/pow(fac1,7.0/2.0));

}

/*
**Second derivative with respect to r for the Plummer density profile
*/


DOUBLE d2rhoPlummerdr2(DOUBLE r, const SI *si){

  DOUBLE rp;
  DOUBLE fac1,fac2,fac3;
  SP *sp;


  sp = si->sp;
  rp=sp->rp;

  fac1 = 1 + (r*r)/(rp*rp);
  fac2 = 5.0/(rp*rp);
  fac3 = (35.0*r*r)/(rp*rp*rp*rp);	

  return (fac3/pow(fac1,9.0/2.0) - fac2/pow(fac1,7.0/2.0));

}


/*
** Hernquist density profile
*/


DOUBLE rhoHernquist(DOUBLE r, const SI *si){

  DOUBLE rhern;
  DOUBLE fac1;
  SP *sp;


  sp = si->sp;
  rhern=sp->rhern;

  fac1 = 1 + r/rhern;

  return (rhern/(r*fac1*fac1*fac1));

}

/*
**Approximate mass enclosed at small r for the Hernquist model
*/

DOUBLE MencInnerHernquist(DOUBLE r, const SI *si){

  DOUBLE rhern;
  DOUBLE fac1,fac2,fac3,fac4;
  SP *sp;

  sp = si->sp;
  rhern=sp->rhern;

  fac1 = rhern*r*r/2.0;
  fac2 = r*r*r;
  fac3 = 3.0*r*r*r*r/(2.0*rhern);
  fac4 = 2.0*r*r*r*r*r/(rhern*rhern);

  return ( 4*M_PI*(fac1-fac2 + fac3 -fac4));
}


/*
** Derivative with respect to r for the Hernquist density profile
*/


DOUBLE drhoHernquistdr(DOUBLE r, const SI *si){

  DOUBLE rhern;
  DOUBLE fac1,fac2,fac3;
  SP *sp;


  sp = si->sp;
  rhern=sp->rhern;

  fac1 = 1 + r/rhern;
  fac2 = 3.0/(r*fac1*fac1*fac1*fac1);
  fac3 = rhern/(r*r*fac1*fac1*fac1);
	
  return (-fac2-fac3);

}

/*
**Second derivative with respect to r for the Hernquist density profile
*/


DOUBLE d2rhoHernquistdr2(DOUBLE r, const SI *si){

  DOUBLE rhern;
  DOUBLE fac1,fac2,fac3,fac4;
  SP *sp;


  sp = si->sp;
  rhern=sp->rhern;

  fac1 = 1 + r/rhern;
  fac2 = 12.0/(rhern*r*fac1*fac1*fac1*fac1*fac1);
  fac3 = 6.0/(r*r*fac1*fac1*fac1*fac1);
  fac4 = 2.0*rhern/(r*r*r*fac1*fac1*fac1);	

  return (fac2+fac3+fac4);

}

/*
**Integrand for stellar enclosed mass integral
*/

DOUBLE integrandMencStar(DOUBLE r, const SI *si) {
	
  return (4.0*M_PI*r*r*(rhoStar(r,si)));
}


/*
**Mass enclosed for the stellar system
*/

DOUBLE IMencStar(DOUBLE rin,DOUBLE r,const SI *si){

  return (si->sp->K*integral(integrandMencStar,rin,r,si));

}

/*
**Approximate mass enclosed at small r for the Stellar Model
*/

DOUBLE MencInnerStar(DOUBLE r, const SI *si){
		
  if (si->king_flag == 1){
    return (MencInnerKing(r,si));
  }else if (si->plummer_flag == 1){
    return  (MencInnerPlummer(r,si));
  }else if (si->hernquist_flag == 1){
    return  (MencInnerHernquist(r,si));
  }  else { 
    fprintf(stderr,"WARNING: MencInnerStar returned 0 \n");
    return (0.0);
  }
}


/*
**The following 3 functions select the stellar density profile and its derivatives.
*/


DOUBLE rhoStar(DOUBLE r, const SI *si){

  if (si->king_flag == 1){
    return (rhoKing(r,si));
  }else if (si->plummer_flag == 1){
    return (rhoPlummer(r,si));
  }else if (si->hernquist_flag == 1){
    return (rhoHernquist(r,si));
  } else { 
    fprintf(stderr,"WARNING: rhoStar returned 0 \n");
    return (0.0);
  }
}

DOUBLE drhoStardr(DOUBLE r, const SI *si){

  if (si->king_flag == 1){
    return (drhoKingdr(r,si));
  }else if (si->plummer_flag == 1){
    return  (drhoPlummerdr(r,si));
  }else if (si->hernquist_flag == 1){
    return  (drhoHernquistdr(r,si));
  } else { 
    fprintf(stderr,"WARNING: drhoStardr returned 0 \n");
    return (0.0);
  }
}

DOUBLE d2rhoStardr2(DOUBLE r, const SI *si){

  if (si->king_flag == 1){
    return (d2rhoKingdr2(r,si));
  }else if (si->plummer_flag == 1){
    return  (d2rhoPlummerdr2(r,si));
  }else if (si->hernquist_flag == 1){
    return  (d2rhoHernquistdr2(r,si));
  } else { 
    fprintf(stderr,"WARNING: d2rhoStardr2 returned 0 \n");
    return (0.0);
  }
}

/*
**Second derivative with respect to phi for the stellar density profile
*/

DOUBLE d2rhoStardPhi2(DOUBLE r, const SI *si) {
	
  DOUBLE rt,K;
  DOUBLE Mencr,rhor;
  DOUBLE fac1, fac2, result;
  SP *sp;

  sp = si->sp;
  rt=sp->rt;
  K = sp->K;

  if (si->king_flag == 1){
    if (r >= rt*0.9999) return(0.0);
  }

  if(si->nostarpot_flag == 1){
    Mencr = MencHalo(r,si);
    rhor = rho(r,si);
  }	
  else if(si->halo_flag == 1){
    Mencr = Menc(r,si);
    rhor = rho(r,si) + K*rhoStar(r,si);
  }
  else {
    Mencr = MencStar(r,si);
    rhor = K*rhoStar(r,si);
  }	
				
  fac1 = r*r/(G*G*Mencr*Mencr);
  fac2 = 2*r-4*M_PI*r*r*r*r*rhor/Mencr;
  result = fac1*(r*r*K*d2rhoStardr2(r,si)+fac2*K*drhoStardr(r,si));
  return (result);	
}


DOUBLE integrandJeansStar(DOUBLE r, const SI *si){

  DOUBLE K,Menc;
  SP *sp;

  sp = si->sp;
  K = sp->K;


  Menc = exp(lininterpolate(NGRIDR,si->logr,si->logMenc,log(r)));

  return (K*rhoStar(r,si)*Menc/(r*r));

}
