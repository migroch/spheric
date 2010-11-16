/* 
** definitions.h
**
** Definitions for spherIC_0.1
*/

/*
 *These are the number of grid points, and the inner and outer radius for the
 *output density profile.
 */
#define ngrid_profile 30
#define rin_profile 0.1
#define rout_profile 5
/*
 *Number of grid points for the grid in r and the grid in the distribution
 *function.
 */
#define NGRIDR 2001
#define NGRIDDF 1001

/*This factor converts velocities from internal units (G = 1) to km/s. velConvert = sqrt(G)*/
#define velConvert 2.07382e2  

/*
 *Stuff
 */
#define NINTMIN 5
#define NINTMAX 28
#define NINTMINDF 8
#define NINTMAXDF 28
#define TOL 1e-10
#define TOLDF 1e-4
#define TOLLININT 1e-10
#define DFFAILUREMAX 1e20
#define FACTORRINNER 1e-3
#define FACTORROUTER 1e10
#define FACTORROUTERSTAR 1e10	
#define SBI 1e100
#define CutoffFac 0.3
#define G 1
#define STRINGSIZE 50
#define INT int
#define FLOAT float
#define DOUBLE double
#define CHAR char
#define OFI1 "%d"
#define OFI2 "%-3d"
#define OFI3 "%-10d"
#define OFD1 "%g"
#define OFD2 "%8.7e"
#define OFD3 "%14.7e"
#define OFD4 "%+8.7e"
#define OFD5 "%16.15e"
#define OFD6 "%+16.15e"

typedef struct gridr {

  DOUBLE r[NGRIDR];
  DOUBLE logr[NGRIDR];
  DOUBLE rho[NGRIDR];
  DOUBLE logrho[NGRIDR];
  DOUBLE rhoHalo[NGRIDR];
  DOUBLE logrhoHalo[NGRIDR];
  DOUBLE rhoenc[NGRIDR];
  DOUBLE logrhoenc[NGRIDR];
  DOUBLE rhoencHalo[NGRIDR];
  DOUBLE logrhoencHalo[NGRIDR];
  DOUBLE rhoStar[NGRIDR];
  DOUBLE logrhoStar[NGRIDR];
  DOUBLE rhoencStar[NGRIDR];
  DOUBLE logrhoencStar[NGRIDR];
  DOUBLE MencStar[NGRIDR];
  DOUBLE logMencStar[NGRIDR];
  DOUBLE Menc[NGRIDR];
  DOUBLE logMenc[NGRIDR];
  DOUBLE MencHalo[NGRIDR];
  DOUBLE logMencHalo[NGRIDR];
  DOUBLE Pot[NGRIDR];
  DOUBLE logPot[NGRIDR];
  DOUBLE Potoutr[NGRIDR];
  DOUBLE eqrvcmax[NGRIDR];
} GRIDR;

typedef struct griddf {

  DOUBLE r[NGRIDDF];
  DOUBLE logr[NGRIDDF];
  DOUBLE E[NGRIDDF];
  DOUBLE logE[NGRIDDF];
  DOUBLE fE[NGRIDDF];
  DOUBLE logfE[NGRIDDF];
  DOUBLE fEstar[NGRIDDF];
  DOUBLE logfEstar[NGRIDDF];
  DOUBLE gE[NGRIDDF];
  DOUBLE loggE[NGRIDDF];
} GRIDDF;

typedef struct systemparameters {

  //halo system parameters

  DOUBLE alpha;
  DOUBLE beta;
  DOUBLE gamma;
  DOUBLE delta;
  DOUBLE M;
  DOUBLE rho0;
  DOUBLE rs;
  DOUBLE rvir;
  DOUBLE rcutoff;
  DOUBLE rdecay;
  DOUBLE rhalf;

  //stellar system parameters
		
  DOUBLE rt;
  DOUBLE rc;	
  DOUBLE rp;
  DOUBLE rhern;
  DOUBLE K;
  DOUBLE Mstar;
  DOUBLE rhalfStar;			
} SP;
    
typedef struct particle {

  DOUBLE r[4];
  DOUBLE v[4];
  DOUBLE mass;
  DOUBLE E;
  DOUBLE star_flag;
  unsigned int index;	
} PARTICLE;

typedef struct stuff {

  INT N;
  DOUBLE Mp;
  DOUBLE Ekin;
  DOUBLE Epot;
  DOUBLE Etot;
  DOUBLE Cr[4];
  DOUBLE Cv[4];
} STUFF;

typedef struct systeminfo {

  INT N0;
  INT N;
  INT Nstar;	
  DOUBLE soft0;
  DOUBLE Mmin;
  DOUBLE Mmax;
  DOUBLE MminStar;
  DOUBLE MmaxStar;
  DOUBLE mass;
  DOUBLE massStar;
  DOUBLE rimp;
  DOUBLE r1;
  DOUBLE r100;
  DOUBLE logr[NGRIDR];
  DOUBLE logMenc[NGRIDR];
  DOUBLE logrhoenc[NGRIDR];
  DOUBLE rinner;
  DOUBLE router;
  DOUBLE routerStar;
  DOUBLE deltapos[3];
  DOUBLE deltavel[3];
  //flags
  INT halo_flag;
  INT king_flag;
  INT plummer_flag;
  INT hernquist_flag;
  INT stars_flag;
  INT nostarpot_flag;
  INT dorvirexact;
  //Structure
  SP *sp;
  PARTICLE *p;
  PARTICLE *pstar;
  GRIDDF *griddf;
  STUFF *stuff;
  GRIDR *gridr;
  //name
  CHAR systemname[STRINGSIZE];
} SI;


typedef struct gadget_header {
    
  unsigned int npart[6];    // number of particles of each type in this file.
  DOUBLE mass[6];           // mass of particles of each type. If 0, then the masses are explicitly stored in the-
                            // mass-block of the snapshot file, otherwise they are omitted.
  DOUBLE time;              // time of snapshot file.
  DOUBLE redshift;          // redshift of snapshot file.
  INT flag_sfr;             // flags whether the simulation was including star formation 
  INT flag_feedback;        // flags whether feedback was included (obsolete) 
  INT npartTotal[6];        // total number of particles of each type in this snapshot This can be different from npart 
                            // if one is dealing with a multi-file snapshot.
  INT flag_cooling;         // flags whether cooling was included.
  INT num_files;            // number of files in multi-file snapshot.
  DOUBLE BoxSize;           // box-size of simulation in case periodic boundaries were used.
  DOUBLE Omega0;            // matter density in units of critical density.
  DOUBLE OmegaLambda;       // cosmological constant parameter.
  DOUBLE HubbleParam;       //Hubble parameter in units of 100 km/sec/Mpc.
  INT flag_stellarage;      //flags whether the file contains formation times of star particles.
  INT flag_metals;          //flags whether the file contains metallicity values for gas and star particles.
  INT npartTotalHighWord[6];  //High word of the total number of particles of each type.
  INT  flag_entropy_instead_u;    // flags that IC-file contains entropy instead of u.
  CHAR fill[60];                       // fills to 256 Bytes.
} GH;

/*
** Tipsy definitions
*/

typedef struct tipsy {

  double time;
  int ntotal;
  int ndim;
  int ngas;
  int ndark;
  int nstar;
} TIPSY_HEADER;

typedef struct gas_particle {

  float mass;
#ifdef DPP
  double pos[3];
#else
  float pos[3];
#endif
  float vel[3];
  float rho;
  float temp;
  float hsmooth;
  float metals;
  float phi;
} GAS_PARTICLE;

typedef struct dark_particle {

  float mass;
#ifdef DPP
  double pos[3];
#else
  float pos[3];
#endif
  float vel[3];
  float eps;
  float phi;
} DARK_PARTICLE;

typedef struct star_particle {

  float mass;
#ifdef DPP
  double pos[3];
#else
  float pos[3];
#endif
  float vel[3];
  float metals;
  float tform;
  float eps;
  float phi;
} STAR_PARTICLE;

typedef struct tipsy_structure {

  TIPSY_HEADER *th;
  GAS_PARTICLE *gp;
  DARK_PARTICLE *dp;
  STAR_PARTICLE *sp;
} TIPSY_STRUCTURE;
