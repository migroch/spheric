/* 
** functions.h
**
** Header file for functions.c
*/


DOUBLE Ecosmo(const SI (*));
DOUBLE integral(DOUBLE (*)(DOUBLE, const SI(*)), DOUBLE, DOUBLE, const SI (*));
DOUBLE trapez(DOUBLE (*)(DOUBLE, const SI(*)), DOUBLE, DOUBLE, INT, const SI (*));
DOUBLE integraldf(DOUBLE (*)(DOUBLE,const SI (*)),INT, const SI (*));
DOUBLE trapezdf(DOUBLE (*)(DOUBLE, const SI (*)),INT, INT, const SI (*));
INT locate(INT, const DOUBLE (*), DOUBLE);
DOUBLE lininterpolate(INT, const DOUBLE (*), const DOUBLE (*), DOUBLE);
DOUBLE rand01();
DOUBLE rho(DOUBLE, const SI (*));
DOUBLE drhodr(DOUBLE, const SI (*));
DOUBLE d2rhodr2(DOUBLE, const SI (*));
DOUBLE d2rhodPhi2(DOUBLE, const SI (*));
DOUBLE dlrhodlr(DOUBLE, const SI (*));
DOUBLE eta(DOUBLE, const SI (*));
DOUBLE detadr(DOUBLE, const SI (*));
DOUBLE tau(DOUBLE, const SI (*));
DOUBLE integrandIM(DOUBLE, const SI (*));
DOUBLE integrandIMcutoff(DOUBLE, const SI (*));
DOUBLE integrandMenc(DOUBLE, const SI (*));
DOUBLE integrandPotHalo(DOUBLE,const SI (*));
DOUBLE integrandPotStar(DOUBLE,const SI (*));
INT split(INT, DOUBLE, const SI (*));
DOUBLE Menc(DOUBLE, const SI (*));
DOUBLE MencHalo(DOUBLE, const SI (*));
DOUBLE MencStar(DOUBLE, const SI (*));
DOUBLE Pot(DOUBLE, const SI (*));
DOUBLE vescape(DOUBLE, const SI (*));
DOUBLE vescapeStar(DOUBLE r,const SI (*));
DOUBLE Tdyn(DOUBLE, const SI (*));
DOUBLE f1(DOUBLE, const SI (*));
DOUBLE f2(DOUBLE, const SI (*));
DOUBLE f1star(DOUBLE, const SI (*));
DOUBLE f2star(DOUBLE, const SI (*));
DOUBLE xKing (DOUBLE, const SI (*));
DOUBLE rhoKing(DOUBLE, const SI (*));
DOUBLE MencInnerKing(DOUBLE, const SI (*));
DOUBLE dxKingdr(DOUBLE, const SI (*));
DOUBLE d2xKingdr2(DOUBLE, const SI (*));
DOUBLE drhoKingdxKing(DOUBLE, const SI (*));
DOUBLE d2rhoKingdxKing2(DOUBLE, const SI (*));
DOUBLE drhoKingdr(DOUBLE, const SI (*));
DOUBLE d2rhoKingdr2(DOUBLE, const SI (*));
DOUBLE rhoPlummer(DOUBLE, const SI (*));
DOUBLE MencInnerPlummer(DOUBLE, const SI (*));
DOUBLE drhoPlummerdr(DOUBLE, const SI (*));
DOUBLE d2rhoPlummerdr2(DOUBLE, const SI (*));
DOUBLE rhoHernquist(DOUBLE, const SI (*));
DOUBLE MencInnerHernquist(DOUBLE, const SI (*));
DOUBLE drhoHernquistdr(DOUBLE, const SI (*));
DOUBLE d2rhoHernquistdr2(DOUBLE, const SI (*));
DOUBLE rhoStar(DOUBLE, const SI (*));
DOUBLE drhoStardr(DOUBLE, const SI (*));
DOUBLE d2rhoStardr2(DOUBLE, const SI (*));
DOUBLE integrandMencStar(DOUBLE, const SI (*));
DOUBLE IMencStar(DOUBLE,DOUBLE,const SI (*));
DOUBLE MencInnerStar(DOUBLE, const SI (*));
DOUBLE d2rhoStardPhi2(DOUBLE, const SI (*));
DOUBLE gE(DOUBLE,DOUBLE,const SI (*));
DOUBLE integrandgE(DOUBLE,DOUBLE,const SI (*));
DOUBLE integralgE(DOUBLE (*)(DOUBLE,DOUBLE,const SI (*)),DOUBLE,DOUBLE,DOUBLE,const SI (*));
DOUBLE trapezgE(DOUBLE (*)(DOUBLE, DOUBLE,const SI (*)),DOUBLE,DOUBLE, DOUBLE, INT,const SI (*));
DOUBLE gEi(DOUBLE , const SI (*));
DOUBLE dMdE(DOUBLE,const SI (*));
DOUBLE fstar1(DOUBLE, const SI (*));
DOUBLE fstar2(DOUBLE, const SI (*));
DOUBLE dMstardE(DOUBLE,const SI (*));
DOUBLE integrandJeansStar(DOUBLE,const SI (*));
