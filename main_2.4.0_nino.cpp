/*####################################################################
 ###  TROLL
 ###  Individual-based forest dynamics simulator
 ###    Version 1: Jerome Chave
 ###    Version 2.1 & 2.2: Isabelle Marechaux & Jerome Chave
 ###    Version 2.3 onwards: Isabelle Marechaux, Fabian Fischer, Jerome Chave
 ###
 ###  History:
 ###    version 0.1 --- JC - 22 Sep 97
 ###    version 0.2 --- JC - 06 Oct 97
 ###    version 0.3 --- JC - 11-14 Nov 97
 ###    version 1.0 --- JC - stable version Chave, Ecological Modelling (1999)
 ###    version 1.1 --- JC - 02-30 Sep 98
 ###    version 1.2 --- JC - 22 Jan 00
 ###    version 1.3 --- JC - 28 Sep 01 stable version Chave, American Naturalist (2001)
 ###
 ###    version 2.0 --- JC - 23 Mar 11 (physiology-based version, translation of comments into English)
 ###    version 2.01 --- IM - oct-nov 13
 ###    version 2.02 --- IM - apr-may 2015
 ###    version 2.03 --- IM - jul 2015
 ###    version 2.04 --- IM - jul 2015 (monthly timestep)
 ###    version 2.1  --- IM - dec 2015 (undef/defined alternative versions)
 ###    version 2.11 --- IM - jan 2016 timestep better used
 ###    version 2.12 --- JC - jan 2016 porting to GitHub for social coding, check of the MPI routines and update, new header for code, trivia: reindentation (orphan lines removed)
 ###    version 2.2  --- IM - may 2016 core changes in: daily coupling with environment; respiration; treefall module
 ###    version 2.3  --- FF - oct-nov 2016: bug fixing (including UpdateSeed() bug), general reworking of code, changes in CalcLAI(), initialisation from data, toolbox with alternative fine flux calculation (cf. end of script)
 ###    version 2.3.0 --- IM & JC - janv-feb 2017: new tree size threshold to maturity (dbhmature), changes in input file structure, corrections of temperature dependencies in GPPleaf. Code acceleration (use of lookup tables instead of functions in the calculations of the Fluxh() and GPPleaf() routines; faster whole-tree GPP calculation. This results in an increase in speed of a factor 4.
 ###    version 2.3.1 --- IM & JC - feb-mar 2017: fast seed dispersal introduces the dispersal cell, or dcell concept (square subplot of linear size length_dcell)
 ###    version 2.4.0 --- FF - feb-mar 2018: intraspecific variation, several changes to remove strong discretisation effects (crown expansion, flux averaging, proportional GPP allocation), treefall functions slightly rewritten, relative height-diameter allometry, crown shape parameter to adjust crown volume
 
 ###    code acceleration FASTGPP now uses dailyGPPcrown() instead of dailyGPPleaf()
 ###    DCELL/INTRASPECIFIC: GNU scientific library is needed -- on osx, type "brew install gsl"
 ###    DCELL/INTRASPECIFIC: Compile command (osx/linux):
 ###                    g++ -O3 -Wall -o troll main_xx.cpp -lgsl -lgslcblas -lm
 ###    DCELL/INTRASPECIFIC: Code profiling: g++ -O3 -Wall -o troll main_xx.cpp -lgsl -lgslcblas -lm -g -pg
 ###    INTRASPECIFIC_covariance: a recent version of the GNU scientific library is needed for the gsl_linalg.h header fie
 ###
 ####################################################################*/


/*
 Glossary: MPI = Message Passing Interface. Software for sharing information     f
 across processors in parallel computers. If global variable MPI is not defined,
 TROLL functions on one processor only.
 */

#undef MPI              /* if flag MPI defined, parallel routines (MPI software) are switched on */
#undef easyMPI          /* if flag easyMPI defined, parallel routine for fast search of parameter space are switched on */
#undef toolbox          /* never to be defined! Toolbox is an assortment of alternative formulations of TROLL procedures, attached to the code */
#undef DCELL            /* this explores the need for an intermediate grid cell size */

#define INTRASPECIFIC               /* new in v.2.4.0: introduces variation around all species-specific traits. All traits are distributed lognormally, except for wood density which is assumed to have a normal distribution (distributions for FG inferred from Bridge data set (Baraloto et al. 2012)) */

#define INTRASPECIFIC_covariance    /* new in v.2.4.0: only works with INTRASPECIFIC. Introduces covariance for variation around traits such as N, P and LMA, and crown radius and height. Possible to deactivate, but requires change in input sheet (could be done differently). ToDo: if covariance is zero, option has to be deactivated, as Cholesky matrix is not pos. definite and produces an error - can be prevented by introducing a condition/catching the error */

#undef FLUX_AVG                    /* new in v.2.4.0: replaces calculations of PPFD, VPD, and T at top of layer by average across a layer, avoids runaway due to unconstrained top layers, developed only for 1m height resolution in v.2.4.0, but easily extended to other resolution */

#undef GPP_PROPORTIONAL            /* new in v.2.4.0: if defined, GPP is weighted by separately for top and bottom layers of tree crown which are usually smaller than the other layers */

#undef CROWN_EXPANSION                /* new in v.2.4.0: replaces strictly symmetric growth (with large jumps in crown area) by gradual filling of crown inside out, lookup table is currently hardcoded to a maximum crown diameter of 51, should probably be set higher to accomodate for extreme parameter combinations */

#undef CROWN_SHAPE                /* new in v.2.4.0: only works with CROWN_EXPANSION. Crown shape can be varied: rotated trapezoid above and inverse rotated trapezoid below, where the trapezoid above has a steeper slope to account for top-heaviness of trees. The rotated trapezoid below quickly approaches a cone. The slopes are governed by a global shape parameter which in the limit 1 returns the cylinder. Main idea: Empirically measured crown depth and crown diameter are usually the maximum extensions of a crown, so cylindric volumes will considerably overestimate the space taken up in the canopy. Simple shape parameter allows for assessment of overestimation. Alternative could be: effective cylinder volume (rescale according to a rescale parameter) or more realistic shapes (ellipsoid), but: the former misrepresents the crown radius and depth, the latter comes with additional computational effort */

#define ALLOM_relative              /* new in v.2.4.0: if activated, the Michaelis Menten parameters are calculated according to the following formula: h_rel = hlim * dbh_rel / (dbh_rel + ah) where dbh_rel = dbh/dmax and h_rel = h/h_dmax. This has to be rescaled for input in TROLL by multiplying with h_dmax so that we get hmax = hlim*h_dmax and h = hmax * dbh_rel / (dbh_rel + ah). Nota bene: h_dmax is the height at dmax (and not the maximum height of the species) and both dmax and hmax are directly taken from the data. In future versions, dmax and h_dmax should be estimated as well and included as parameters in the hierarchical model. We could take dmax as the 99th percentile of a lognormal distribution and knowing that at dbh_relmax = 1, h = 1 we infer 1 = hlim / (1 + ah), translating into: hlim = 1 + ah, and finally: h = h_dmax * (1 + ah) * dbh_rel / (dbh_rel + ah) */

#undef LL_limit                     /* new in v.2.4.0: lower limit to leaflifespan (currently set to 3 months), quick fix, if leafdem_resolution = 1, but could also be justified from rarity of low leaflifespans in French Guiana */


/* Libraries */
# include <cstdio>
# include <iostream>
# include <fstream>
# include <cstdlib>
# include <string>
# include <limits>
# include <ctime>
# include <cmath>
# ifdef MPI
# include "mpi.h"
# endif

#if defined (INTRASPECIFIC) || defined(DCELL)
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

#ifdef INTRASPECIFIC_covariance
#include <gsl/gsl_linalg.h>
#endif
#endif

// NINO
#include <iomanip>
#include <array>
#include <sstream>
#include <stdio.h>

using namespace std;

/* Global constants (e.g. PI and various derivatives...) */
# define PI 3.141592654
# define twoPi 6.2831853071
# define Pis2 1.570796327
# define iPi 0.3183099
char buffer[256], inputfile[256], inputfile_climate[256], inputfile_sylviculture[256], outputinfo[256], inputfile_data[256], *bufi(0), *bufi_climate(0), *buf(0), *bufi_data(0), *bufi_sylviculture(0);

/* random number generators */
double genrand2(void);
void sgenrand2(unsigned long);
unsigned long genrand2i(void);
void sgenrand2i(unsigned long);

/* file output streams */

fstream out,out2, out3; /* NINO : added out3 */
fstream output[40]; /* NINO : should this size be modified ? */

/****************/
/* User control */
/****************/

/* options can be turned on (1) or off (0). This comes, however, at computational costs. Where routines have to be called frequently, if-conditioning should be done as far outside the loop as possible (e.g. for DAILYLIGHT outside voxel loops) */
/* currenly, options are set below, but inclusion in parameter sheet needed (for control from R) */

bool
_NONRANDOM=1,           /* If _NONRANDOM == 1, the seeds for the random number generators will be kept fixed at 1, for bug fixing */
_FASTGPP=0,             /* This defines an option to compute only GPP from the topmost value of PPFD and GPP, instead of looping within the crown. Much faster and more accurate [but lateral competition?] */
_BASICTREEFALL=1,       /* if defined: treefall is a source of tree death (and if TREEFALL not defined, this is modeled through simple comparison between tree height and a threshold t_Ct, if not defined, treefall is not represented as a separated and independent source of death, but instead, all tree death are due to the deathrate value) */
_TREEFALL=0,            /* computation of the force field if TREEFALL is defined, neighboring trees contribute to fell each tree */
_DAILYLIGHT=1,          /* if defined: the rate of carbon assimilation integrates an average daily fluctuation of light (thanks to GPPDaily). Should be defined to ensure an appropriate use of Farquhar model */
_SEEDTRADEOFF=0,        /* if defined: the number of seeds produced by each tree is determined by the tree NPP allocated to reproduction and the species seed mass, otherwise the number of seeds is fixed; besides, seedling recruitment in one site is not made by randomly and 'equiprobably' picking one species among the seeds present at that site but the probability of recruitment among the present seeds is proportional to the number of seeds (in s_Seed[site]) time the seed mass of each species */
_NDD=0,                 /* if defined, negative density dependant processes affect both the probability of seedling recruitment and the local tree death rate. The term of density-dependance is computed as the sum of conspecific tree basal area divided by their distance to the focal tree within a neighbourhood (circle of radius 15m) */
_OUTPUT_reduced=1,      /* reduced set of ouput files */
_OUTPUT_last100=0,      /* output that tracks the last 100 years of the simulation for the 1whole grid (2D) */
_FromData=1,            /* if defined, an additional input file can be provided to start simulations from an existing data set or a simulated data set (5 parameters are needed: x and y coordinates, dbh, species_label, species */

/* NINO : ADDED THIS */
_INPUT_fullFinal=1,	   /* works only with _FromData=1. initializes simulation from a full final pattern with all tree attributes */
_OUTPUT_fullFinal=1,	/* output of full final pattern with all tree attributes */


_DISTURBANCE=0,			/* if defined: implementation of a basic perturbance module at a given iteration step in waiting for a more sofisticated sylviculture module */
_LOGGING=0;			/* if defined: implementation of a selective logging module imulating sylviculture as it is happening in french guiana */

/********************************/
/* Parameters of the simulation */
/********************************/

int sites,          /* number of sites */
cols,               /* nb of columns */
rows,               /* nb of lines */
numesp,             /* nb of species */
iterperyear,        /* nb of iter in a year (=12 if monthly timestep, =365 if daily timestep) */
nbiter,             /* total nb of timesteps */
iter,               /* current timestep */
nbout,              /* nb of outputs */
freqout;            /* frequency HDF outputs */

#if defined(INTRASPECIFIC) || defined(DCELL)
gsl_rng *gslrand;
#ifdef INTRASPECIFIC_covariance
gsl_matrix *mcov_N_P_LMA;  /* covariance matrices for crown radius - tree height covariance and leaf_properties, respectively */
gsl_vector *mu_N_P_LMA, *variation_N_P_LMA; /* the mean values of the distributions as well as the result vector for the multivariate draw */
#endif
#endif

#ifdef DCELL
int length_dcell,   /* v2.3.1 linear size of a dcell */
linear_nb_dcells,    /* linear number of dcells note that nbdcells = linear_nb_dcells^2 */
sites_per_dcell,    /* number of sites per dcell */
nbdcells;           /* total number of dcells */
int **MAP_DCELL(0);   /* list of sites per dcell (for fast fillin) */
double *prior_DCELL(0); /* prior for picking the dispersal sites within the dcell -- vector of size sites_per_dcell and with all entries equal to 1/sites_per_dcell (because the dispersal is equiprobable) */
unsigned int *post_DCELL(0); /* number of seeds per site within the dcell */
double *prior_GERM(0); /* prior for picking the germination event within a site -- vector of size numesp and with entries equal to the number of seeds multiplied by seedsize */
unsigned int *post_GERM(0); /* vector with only one non-null entry (the successful germination event) */
#endif

int HEIGHT,         /* max height (in m) */
dbhmaxincm,         /* max DBH times 100 (ie dbh in cm *100 = in meters) */
RMAX,               /* max crown radius */
SBORD,              /* RMAX*cols */
leafdem_resolution; /* resolution for leaf demography model */

float NV,           /* nb cells per m (vertical) */
NH,                 /* nb cells per m (horizontal) */
LV,                 /* LV = 1.0/NV */
LH,                 /* LH = 1.0/NH */
timestep;           /* duration of one timestep (in years)=1/iterperyear */

float p_nonvert,    /* ratio of non-vertical incident light */
Cseedrain,          /* constant used to scale total seed rain per hectare across species */
nbs0,               /* nb of seeds produced and dispersed by each mature tree when SEEDTRADEOFF is not defined */
Cair,               /* atmosphericCO2 concentration, if we aim at making CO2 vary (scenarios), CO2 will have to have the same status as other climatic variables  */
iCair;              /* inverse of Cair */

#ifdef CROWN_SHAPE
float shape_crown;  /* this corresponds to the ratio of radius at the crown top to the radius at the crown center (assumed to be the maximum) */
#endif


/* new version 2.2 */
float daily_light[24];    /* normalized (ie between 0 and 1) daily light variation (used if DAILYLIGHT defined) */
float daily_vpd[24];      /* normalized (ie between 0 and 1) daily vpd variation (used if DAILYLIGHT defined) */
float daily_T[24];        /* normalized (ie between 0 and 1) daily T variation (used if DAILYLIGHTdefined) */

/*parameter for NDD*/
float R,                    /* distance beyond which NDD effect is not accounted anymore*/
deltaR,                     /* NDD strength parameter in recruitment */
deltaD,                     /* NDD strength parameter in deathrate */
BAtot;

/* NINO : re-added sylvain's disturbance module */
int disturb_iter;			/* iteration step where the disturbation occure */
float disturb_intensity;	/* intensity of disturbance in percent of BA */

/* NINO : re-added sylvain's logging module */
int numespharvestable,		/* number of harvestable species */
designated_volume,			/* volume designated for harvesting in m3/ha */
harvested_volume;       	/* volume harvested in m3/ha */

/*********************************************/
/* Environmental variables of the simulation */
/*********************************************/

/* Climate input data  */
/* these climate input data are given in the input file, its structure depends on the timestep and scenario used for the simulation */
/* new version 2.2 */
float *Temperature(0);                      /* in degree Celsius */
float *DailyMaxTemperature(0);              /* in degree Celsius */
float *NightTemperature(0);                 /* in degree Celsius */
float *Rainfall(0);                         /* in mm */
float *WindSpeed(0);                        /* in m/s */
float *MaxIrradiance(0);                    /* in W/m2 */
float *MeanIrradiance(0);                   /* in W/m2 */
float *SaturatedVapourPressure(0);          /* in kPa */
float *VapourPressure(0);                   /* in kPa */
float *VapourPressureDeficit(0);            /* in kPa */
float *DailyVapourPressureDeficit(0);       /* in kPa */
float *DailyMaxVapourPressureDeficit(0);    /* in kPa */

/* New in v.2.3.0 */
int nbTbins;                    /*nb of bins for the temperature lookup tables */
float iTaccuracy;                /* inverse of accuracy of a temperature bin (e.g. Taccuracy us 0.1 or 0.5 °C, so iTaccuracy is 10.0 or 2.0, resp) */
float *LookUp_KmT(0);                   /* lookup table for fast comput of Farquhar as a function of T */
/* !! leaf temperature must be comprised between 0°C and 60°C
 (T_leaf is stored every 0.5°C, so 120 values in total */
float *LookUp_GammaT(0);                   /* lookup table for fast comput of Farquhar */
float *LookUp_tempRday(0);                 /* lookup table for fast comput of Farquhar */
float *LookUp_VcmaxT(0);                   /* lookup table for fast comput of Farquhar */
float *LookUp_JmaxT(0);                   /* lookup table for fast comput of Farquhar */
float *LookUp_flux(0);                 /* lookup table for faster computation of PPFD / new in v.2.4: averaging instead of top value*/
float *LookUp_VPD(0);                  /* lookup table for faster computation of VPD / new in v.2.4:averaging instead of top value */
float *LookUp_T(0);                    /* lookup table for faster computation of T / new in v.2.4: averaging instead of top value */
float *LookUp_Rstem(0);                /* lookup table for faster computation of Rstem */
float *LookUp_Rnight(0);                /* lookup table for faster computation of Rstem */

#ifdef CROWN_EXPANSION
int LookUp_Crown_site[2601];               /* new in v.2.4: lookup table to fill crown cylinder sequentially from inside to outside, allowing for smooth crown growth */
#endif

/***** Environmental variables, changed at each timestep *****/
float temp,         /* Temperature */
tmax,               /* Daily max temperature */
tnight,             /* Night mean temperature */
precip,             /* Rainfall  */
WS,                 /* WindSpeed */
Wmax,               /* Daily max irradiance (average for timestep) (in micromol PAR photon/m^2/s)*/
/* used in the photosynthesis part. see if it would not be better to have value in the right unit in the input file, however W/m2 is the common unit of meteo station */
/* below: new version 2.2 */
Wmean,              /* mean irradiance (in W/m2)*/
e_s,                /* SaturatedVapourPressure */
e_a,                /* VapourPressure*/
VPDbasic,           /* VapourPressureDeficit */
VPDday,             /* DailyVapourPressureDeficit */
VPDmax;             /* DailyMaxVapourPressureDeficit */


/****************************************/
/*  Common variables for the species    */
/* (simplifies initial version 170199)  */
/****************************************/

//int Cm;           /* Treefall threshold */
float klight,       /* light absorption rate or extinction cefficient used in Beer-Lambert law to compute light within the canopy */
phi,                /* apparent quantum yield (in micromol C/micromol photon). phi is the quantum yield multiplied by leaf absorbance (in the literature, the quantum yield is often given per absorbed light flux, so one should multiply incident PPFD by leaf absorbance, see Poorter et al American Journal of Botany (assumed to be 0.91 for tropical tree species). Even though holding phi constant across species is widely assumed, in reality, phi varies across species and environmental conditions (see eg Domingues et al 2014 Plant Ecology & Diversity) */
theta,              /* parameter of the Farquhar model set to 0.7 in this version. For some authors, it should be species-dependent or temperature dependent, but these options are not implemented here */
g1,                 /* g1 parameter of Medlyn et al's model of stomatal conductance. v230: defined as a global parameter shared by species, instead of a Species class's variable */
alpha,              /* apparent quantum yield to electron transport in mol e-/mol photons, equal to the true quantum yield multiplied by the leaf absorbance -- v.2.2 */
vC,                 /* variance of treefall threshold */
H0,                 /* initial height (in m) */
DBH0,               /* initial DBH (in m) */
de0,                /* initial crown Crown_Depth (in m) */
de1,                /* Crown_Depth/height slope */
/* fallocwood and falloccanopy new -- v.2.2 */
fallocwood,         /* fraction of biomass allocated to above ground wood (branches+stem) */
falloccanopy,       /* fraction of biomass allocated to canopy (leaves + reproductive organs + twigs) */
dens,               /* initial leaf density (in m^2/m^3) */
ra1,                /* crown radius - dbh slope */
ra0,                /* initial crown radius (in m) */

#ifdef INTRASPECIFIC
/* sigmas for intraspecific variation, currently assumed to be the same for all species */
sigma_height, sigma_CR, sigma_CD, sigma_P, sigma_N, sigma_LMA, sigma_wsg, sigma_dmax,
#ifdef INTRASPECIFIC_covariance
/* correlation / covariance for selected traits */
corr_CR_height, corr_N_P, corr_N_LMA, corr_P_LMA,
cov_N_P, cov_N_LMA, cov_P_LMA,
#endif
#endif

p_tfsecondary,      /* new in v.2.4.0: probability that a death that results from a treefall is a treefall itself */
hurt_decay,         /* new in v.2.4.0: factor by which t_hurt is decaying each timestep */
m,                  /* basal death rate */
m1;                 /* deathrate-wsg slope -- new v.2.2 */


#ifdef INTRASPECIFIC
/* LookUp_tables for intraspecific variation (but maybe drawing from distribution is equally quick?) */
float d_intraspecific_height[100000],
d_intraspecific_CR[100000],
d_intraspecific_CD[100000],
d_intraspecific_P[100000],
d_intraspecific_N[100000],
d_intraspecific_LMA[100000],
d_intraspecific_wsg[100000],
d_intraspecific_dmax[100000];
#endif

float **LAI3D(0);   /* leaf density (per volume unit) */
unsigned short *Thurt[3];            /* Treefall field */
unsigned short *Tlogging[3];            /* Logging tree field (0 for logged, 1 for tracks, 2 for gap damages) NINO*/

int    *SPECIES_GERM (0);
float  *PROB_S (0); /* _SEEDTRADEOFF */
float tempRday;     /* temporary variable used for the computation of leaf day respiration -- new v.2.2 */


/***************/
/* Diagnostics */
/***************/

int nblivetrees,            /* nb live trees for each timestep  */
nbtrees_n10,            /* nb trees dbh > 10 cm, computed at beginning of each timestep */
nbdead_n1,              /* nb deaths dbh > 1 cm, computed at each timestep */
nbdead_n10,             /* nb deaths dbh > 10 cm, computed at each timestep */
nbdead_n30,             /* nb deaths dbh > 30 cm, computed at each timestep */
nbTreefall1,            /* nb treefalls for each timestep (dbh > 1cm), _BASICTREEFALL */
nbTreefall10,           /* nb treefalls for each timestep (dbh > 10 cm), _BASICTREEFALL */
nbTreefall30;           /* nb treefalls for each timestep (dbh > 30 cm), _BASICTREEFALL */
/* NINO : A diagnostics for the volume or number of commercial trees is to be added */

//long int *persist;	/* persistence histogram */
int *nbdbh(0);          /* dbh size distribution */
float *layer(0);        /* vertical LAI histogram */

/**************/
/* Processors */
/**************/

int mpi_rank,mpi_size;
int easympi_rank;


/******************/
/* MPI procedures */
/******************/

#ifdef MPI
unsigned short **LAIc[2];
void
MPI_ShareSeed(unsigned char **,int),
MPI_ShareField(unsigned short **,unsigned short ***,int),
MPI_ShareTreefall(unsigned short **,int);
#endif


/**********************/
/* Simulator routines */
/**********************/

void
Initialise(void),
InitialiseFromData(void),
InitialiseLogging(fstream&),	// _LOGGING /*NINO*/
AllocMem(void),
BirthInit(void),
Evolution(void),
UpdateField(void),
TriggerTreefall(void),              // TROLL v 2.4.0: rewritten UpdateTreefall(), new name for clarity's sake
TriggerTreefallSecondary(void),     // TROLL v 2.4.0: new function, triggering secondary treefalls at the beginning of each round from previous damage
UpdateTree(void),
Average(void),
OutputField(void),
FreeMem(void),
/* NINO */
Disturbance(void),					// _DISTURBANCE
SelectiveLogging(void);				// _LOGGING

/******************************************/
/* Selective logging routines (_LOGGING) */
/****************************************/

void
Designate(),
Select(void),
Rot(void),
Fell(void),
MainTracks(void),
SecondaryTracks(void),
GapDamages(void);

/**********************/
/** Output routines ***/
/**********************/

void
OutputSnapshot(fstream& output),
OutputSnapshotDetail(fstream& output),
OutputSnapshot10cm(fstream& output),
OutputSpeciesParameters(fstream& output),
OutputLAI(fstream& output_CHM, fstream& output_LAI),
OutputLAIFull(fstream& output_LAI3D),
OutputSnapshotFullFinal(fstream& output); /* NINO : fullfinal that needs to be readapted */

/****************************/
/* Various inline functions */
/****************************/


inline float flor(float f) {
    if(f>0.) return f;
    else return 0.;
}
inline float florif(int i) {
    if(i>0) return float(i);
    else return 0.;
}
inline float maxf(float f1, float f2) {
    if(f1>f2) return f1;
    else return f2;
}
inline float minf(float f1, float f2) {
    if(f1<f2) return f1;
    else return f2;
}
inline int min(int i1, int i2) {
    if(i1<i2) return i1;
    else return i2;
}
inline int max(int i1, int i2) {
    if(i1>i2) return i1;
    else return i2;
}
inline int sgn(float f) {
    if(f>0.0) return 1;
    else return -1;
}


/*############################################
 ############################################
 ###########     Species  class   ###########
 ############################################
 ############################################*/

class Species {
    
public:
    int    s_nbind,			/* nb of individuals per species */
    s_dormDuration,         /* seed dormancy duration -- not used in v.2.2 */
    s_nbext;                /* total number of incoming seeds in the simulated plot at each timestep (seed rain) -- v.2.2 */
    char	s_name[256];	/* species name */
    float  s_LCP,			/* light compensation point  (in micromol photon/m^2/s) */
    s_Rdark,                /* dark respiration rate (at PPFD = 0) in micromol C/m^2/s) */
    s_ds,                   /* average dispersal distance */
    //  de1,                /* (crown depth) - height slope deprecated v.2.1 */
    s_dmax,                 /* maximal dbh (given in m) */
    s_hmax,                 /* maximal height (given in m) */
    s_dbh0,                 /* Initial dbh at recruitement, computed for each species with the species-specific allometric relationship at h=H0=1m -- in v230 */
    s_Vcmax,                /* maximal rate of carboxylation, on an area basis (in micromolC/m^2/s) */
    s_Vcmaxm,               /* maximal rate of carboxylation, on a mass basis (in micromolC/g-1/s) */
    s_Jmax,                 /* maximal rate of electron transport, on an area basis (in micromol/m^2/s) */
    s_Jmaxm,                /* maximal rate of electron transport, on a mass basis (in micromol/g-1/s) */
    s_fci,                  /* fraction of CO2 partial pressure in intercellular spaces divided by ambiant CO2 partial pressure (both in microbar, or ppm = micromol/mol) */
    s_Gamma,                /* compensation point for the carboxylation rate, here NORMALIZED by atm CO2 concentration (Cair) */
    s_Km,                   /* apparent kinetic constant for the rubiscco = Kc*(1+[O]/Ko), here normalized by atm CO2 concentration (Cair) */
    //s_d13C,                 /* isotopic carbon discrimination NOW normalized at zero height -- deprecated v.2.2 */
    s_LMA,                  /* leaf mass per area (in g/m^2) */
    s_Nmass,                /* leaf nitrogen concentration (in g/g) v.2.01 */
    s_Pmass,                /* leaf phosphorous concentration (in g/g) v.2.01 */
    s_wsg,                  /* wood specific gravity (in g/cm^3) */
    s_ah,                   /* parameter for allometric height-dbh equation */
    s_seedmass,             /* seed mass, in g (from Baraloto & Forget 2007 dataset, in classes) v.2.3: seeminlgy deprecated in v.2.2, but still necessary for SEEDTRADEOFF */
    s_iseedmass,            /* inverse of seed mass, v.2.3: seeminlgy deprecated in v.2.2, but still necessary for SEEDTRADEOFF */
    //s_factord13Ca,        /* deprecated v.2.2 -- factor used for a previous version of ci/ca ratio computation, from d13C value */
    //s_factord13Cb,        /* deprecated v.2.2 -- factor used for a previous version of ci/ca ratio computation, from d13C value */
    /* Below: new in v.2.2 */
    s_leaflifespan,         /* average leaf lifespan, in month */
    s_time_young,           /* leaf resident time in the young leaf class */
    s_time_mature,          /* leaf resident time in the mature leaf class */
    s_time_old,             /* leaf resident time in the old leaf class */
    s_output_field[24];         /* scalar output fields per species (<24) */

#ifdef DCELL
    int *s_DCELL;	/* number of seeds from the species in each dcell */
    int *s_Seed;	/* presence/absence of seeds at each site; if def SEEDTRADEOFF, the number of seeds */
#else
    int *s_Seed;	/* presence/absence of seeds; if def SEEDTRADEOFF, the number of seeds */
#endif
    
#ifdef MPI
    unsigned char *s_Gc[4]; /* MPI: seeds on neighboring procs */
#endif

	    /* species logging parameters */
	bool s_harvestable;		/*is the species in the list of harvestable species*/
	float s_dbhmin,			/*minimum harvestable diameter*/
	s_dbhmax;				/*maximum harvestable diameter*/
	int s_interest;				/*species interest to filter harvested trees in designated trees*/
    
    
    Species() {
        s_nbind=0;
        s_Seed=0;
#ifdef DCELL
        s_DCELL=0;
#endif
    };                              /* constructor */
    
    virtual ~Species() {
        delete [] s_Seed;
#ifdef DCELL
        delete [] s_DCELL;
#endif
    };                              /* destructor */
    
    void Init(int,fstream&);        /* init Species class */
#ifdef DCELL
    void FillSeed(int,int,int);         /* assigns the produced seeds to s_DCELL */
#else
    void FillSeed(int,int);         /* fills s_Seed field (and s_Gc (MPI)) */
#endif
    void UpdateSeed(void);       /* Updates s_Seed field */
#ifdef MPI
    void AddSeed(void);             /* MPI: adds fields s_Gc  to field s_Seed */
#endif
    
    inline float DeathRateNDD(float, float, int, float); /* _NDD, overloading with function in following line */
    inline float DeathRate(float, float, int);  /* actual death rate -- new argument int v.2.2 */
    inline float GPPleaf(float, float, float);    /* Computation of the light-limited leaf-level NPP per m^2 (in micromol/m^2/s) -- two new arguments float v.2.2 */
    /* Farquhar von Caemmerer Berry model */
    inline float dailyGPPleaf(float, float, float);    /* computation of the daily average assimilation rate, taking into account the daily variation in light, VPD and temperature two new arguments float v.2.2, _DAILYLIGHT */
    inline float dailyGPPcrown(float, float, float, float, float);  /* v.2.4.0: proposition of reinsertion of fastdailyGPPleaf() as dailyGPPcrown(), main reason: despite sharing some code with dailyGPPleaf, structurally very different */
};


/*############################
 ###  Initialize Species  ###
 ###    Species::Init     ###
 ############################*/

void Species::Init(int nesp,fstream& is) {
    
    int site;
    float regionalfreq;     // "regional" relative abundance of the species -- new name v.2.2
    float SLA;              // specific leaf area = 1/LMA
    
    /*** Read parameters ***/
    
    //new input file -- in v230
    is  >> s_name >> s_LMA >> s_Nmass >>  s_Pmass  >> s_wsg >> s_dmax >> s_hmax  >> s_ah >> s_seedmass >> regionalfreq;
    
    // instead of seedmass we are given seedvolume
    // from this we assume a conversion factor of 1 to wet mass (~density of water, makes seeds float)
    // to convert to drymass we use a conversion factor of 0.4 (~40% of the seed are water)
    
    s_seedmass *= 0.4;
    s_iseedmass=1.0/s_seedmass;
    s_ds=40.0;
    
    s_dbh0=s_ah*H0/(s_hmax-H0);
    
#ifdef DCELL
    /* NEW in 2.3.1: input of seeds per timestep and per dcell is assumed proportional to seed size (large seeds are less numerous) this is compensated for by a lower recruitment probability for small seeds compared to large ones. Thus, in essence, we still have the seed number regeneration tradeoff but we allow more stochasticity in this process. Also, at least one seed arrives in the dcell from each species at each timestep (this is the +1 term) */
#endif
    
    // uniform composition of the seed rain -- in v230
    regionalfreq=1.0/float(numesp);
    
    if(_SEEDTRADEOFF){
        s_nbext = (int(regionalfreq*Cseedrain*s_iseedmass)+1);
    }
    else {
        s_nbext = int(regionalfreq*Cseedrain*(sites*LH*LH/10000));
    }
    
    SLA=10000.0/s_LMA;    // computation of specific leaf area in cm^2/g for use in Domingues et al 2010 model of photosynthetic capacities
    s_leaflifespan = pow(10,(2.040816*(2.579713-log10(SLA))));    //this is the expression from Reich et al. 1997 PNAS (provides probably more realistic estimates for species with high LMA).
    //s_leaflifespan=1.5+pow(10,(7.18+3.03*log10(s_LMA*0.0001)));           //this is the expression from Reich et al 1991 Oecologia (San Carlos Rio Negro).
    //s_leaflifespan=0.5+pow(10,(-2.509+1.71*log10(s_LMA)));    //this is the expression from Wright et al 2004 Nature (leaf economics spectrum).

#ifdef LL_limit
    s_leaflifespan = maxf(s_leaflifespan,3.0);
#endif
    
    s_time_young=1.0;
    s_time_mature=s_leaflifespan/3.0;
    s_time_old=s_leaflifespan-s_time_mature-s_time_young;
    
    //cout << "Initialise species: " << s_name << "(" << nesp << ") with leaflifespan (in months): " << s_leaflifespan << " and turnover: " << s_time_young << " (young) | " << s_time_mature << " (mature) | " << s_time_old << " (old) " << endl;
    
    /*** Normalization of the parameters ***/
    /* vertical (NV) and horizontal (NH) scales */
    
    s_ah   *= NV*LH;
    s_ds   *= NH;
    s_hmax *= NV;
    s_dmax *= NH;
    s_dbh0 *= NH;
    
    s_nbind=0;
    s_fci = 0.0;
    
    s_Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(s_Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(s_Pmass*1000.0)+0.25*log10(SLA))));
    // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
    s_Jmaxm=pow(10.0, minf((-1.50+0.41*log10(s_Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(s_Pmass*1000.0)+0.32*log10(SLA))));
    // added as a Species member variable 14-04-2015; this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
    
    s_Vcmax=s_Vcmaxm*s_LMA;
    s_Jmax=s_Jmaxm*s_LMA;
    
    s_Rdark=s_LMA*(8.5341-130.6*s_Nmass-567.0*s_Pmass-0.0137*s_LMA+11.1*s_Vcmaxm+187600.0*s_Nmass*s_Pmass)*0.001;
    
    //s_Rdark corresponds to leaf maintenance respiration. From Table 6 in Atkin et al 2015 New phytologist v.2.0 */
    
    //s_Rdark=(82.36*(s_LMA*1e-3)-0.1561)*(s_LMA*1e-3);                 /* from Domingues et al 2007 */
    //s_Rdark=0.01*s_Vcmax;                                             /* parameterization of Rdark commonly used in vegetation models */
    //s_Rdark=0.02*s_Vcmax-0.01;                                        /* parameterization of von Caemmerer 2000 Table 2.3 page 45 */
    
    s_Gamma = 38.0*iCair;
    
    // s_Gamma at 25°C computed according to von Caemmerer 2000 formula: gamma=Kc*O*0.25/(2*Ko), with Kc=260 microbar, Ko=179mbar and O=210 mbar (the last value is from Farquhar et al 1980, the first two one are from von Caemmerer 2000 table 2.3 page 45). gamma is set to 36.9 on Atkin et al 2015. Could be a global variable. v.2.0
    
    s_LCP = s_Rdark/phi;    /* Computation of the light compensation point from dark respiration and the quantum yield phi. By definition, Rdark is in micromolC/m^2/s and it is used in the Species::NPP() routine */
    
      	/*** Sylviculture parameters NINO***/  
    s_harvestable=0;

#ifdef DCELL
    if (NULL==(s_DCELL = new int[nbdcells])) cerr<<"!!! Mem_Alloc s_DCELLn";
    for(int dcell=0;dcell<nbdcells;dcell++) s_DCELL[dcell]=0;
    /***  field initialization ***/
    if (NULL==(s_Seed = new int[sites])) cerr<<"!!! Mem_Alloc\n";
    for(site=0;site<sites;site++) s_Seed[site]=0;
#else
    /***  field initialization ***/
    if (NULL==(s_Seed = new int[sites])) cerr<<"!!! Mem_Alloc\n";
    for(site=0;site<sites;site++) s_Seed[site]=0;
#endif
    
#ifdef MPI
    for(int i=0;i<4;i++) {
        if (NULL==(s_Gc[i] = new unsigned char[sites])) cerr<<"!!! Mem_Alloc\n";
        for(site=0;site<sites;site++) s_Gc[i][site]=0;
    }
#endif
}


/*############################
 ###     Species Seeds     ###
 ###   Species::FillSeed   ###
 ###  Species::UpdateSeed  ###
 ###   Species::AddSeed    ###
 #############################*/

/*###  Species::FillSeed  ###*/
/* creates one seed at point (col,row) */

#ifdef DCELL
/* in the new approach with a mean field seed flux (DCELL), the function FillSeed has a new
 role: it fills the vector s_DCELL that stores the number of produced seeds per timestep and per dcell */

void Species::FillSeed(int dcol, int drow, int nbs) {
    s_DCELL[dcol+linear_nb_dcells*drow]+=nbs;
}

#else

void Species::FillSeed(int col, int row) {
    int site;
    if(col < cols) {
        if((row >=0) && (row < rows)) {
            site=col+cols*row;
            if(_SEEDTRADEOFF){
                s_Seed[site]++;                         /* ifdef SEEDTRADEOFF, s_Seed[site] is the number of seeds of this species at that site */
            }
            else{
                if(s_Seed[site]!=1) s_Seed[site]=1;     /* If s_Seed[site] = 0, site is not occupied, if s_Seed[site] > 1, s_Seed[site] is the age of the youngest seed  */
            }
        }
        
#ifdef MPI                                       /* on each processor a stripe of forest is simulated.
Nearest neighboring stripes are shared. Rque, this version is not valid ifdef SEEDTRADEOFF */
        else if((row+rows >=0) && (row < 0)) {
            site=col+cols*(row+rows);
            if(s_Gc[0][site]!=1) s_Gc[0][site]=1;
        }
        else if((row >=rows) && (row < 2*rows)) {
            site=col+cols*(row-rows);
            if(s_Gc[1][site]!=1) s_Gc[1][site]=1;
        }
#endif
    }
}
#endif



/*### Species::UpdateSeed ###*/
/* updates s_Seed field */
/* new in v.2.3: not called within loop over sites, instead includes loop --> less function calling, BUT: no check of site occupation anymore, cf. below */

#ifdef DCELL
/* in the new approach with a mean field seed flux (DCELL), the function UpdateSeed
 has a new role: it uses the vector s_DCELL to fill the s_Seed local seed bank */

void Species::UpdateSeed() {
    int site;
    for(site=0;site<sites;site++) s_Seed[site]=0;
    for(int dcell=0;dcell<nbdcells;dcell++){ // loop over dcells
        int localseeds=min(s_DCELL[dcell],sites_per_dcell);
        // store number of seeds in this dcell
        // localseeds is capped by the max number of available sites in the dcell
        s_DCELL[dcell]=0;
        gsl_ran_multinomial(gslrand,sites_per_dcell,localseeds,prior_DCELL,post_DCELL);
        //cerr << "Localseeds in dcell\t" << dcell << " : " << localseeds << endl;
        /*float sumprior=0.0,sumpost=0;
         for(int i=0;i<sites_per_dcell;i++){
         sumprior+=prior_DCELL[i];sumpost+=post_DCELL[i];
         }*/
        //cerr <<"localseeds\t"<< localseeds << "\tsumprior\t"<< sumprior << "\tsumpost\t"<< sumpost << "\n";
        // sample equiprobably all the sites in the dcell
        for(int i=0;i<sites_per_dcell;i++){ // update the s_Seed site
            site=MAP_DCELL[dcell][i];
            s_Seed[site] = post_DCELL[i];
            //cerr << "site\t" << site << "\tdcell\t" << dcell << "\tlocal_site\t" << i << "\tpost_DCELL[i]\t" << post_DCELL[i] << "\ts_Seed[site]\t" << s_Seed[site] << endl;
        }
    }
    /*int summ=0;
     for(site=0;site<sites;site++) summ=summ+s_Seed[site];
     cerr << "Localseeds of species \t" << s_name  << "\t: " << summ << endl; */
    
    
}
#else
void Species::UpdateSeed() {
    
    /* should probably be modified, since as implemented now seeds are erased every timestep (i.e. month in default mode)--> to be discussed */
    
    if(_SEEDTRADEOFF){
        for(int site=0;site<sites;site++){
# ifdef MPI
            s_Gc[0][site]=s_Gc[1][site]=s_Gc[2][site]=s_Gc[3][site]=0;
#endif
            s_Seed[site]=0;
        }
    }
    
    else{
        /* new in v.2.3: version 2.2 checked whether site was occupied by tree: T[site].t_age>0) s_Seed[site]=0;     */
        /* v.2.3 does not do this within UpdateSeed() anymore. Instead, it sets all occupied sites to zero directly within UpdateTree() */
        for(int site=0;site<sites;site++){
            
# ifdef MPI
            s_Gc[0][site]=s_Gc[1][site]=s_Gc[2][site]=s_Gc[3][site]=0;
#endif
            /* seed bank ages or disappears */
            if(s_Seed[site]==s_dormDuration) s_Seed[site]=0;
            else if(s_Seed[site]!=0) s_Seed[site]++;            // v.2.3: bug fix: before, procedure was not restricted to existing seeds, therefore creation of seeds
        }
    }
}
#endif

#ifdef MPI
/* deprecated in v.2 -- needs a new concept for spatial parallelization -- hopefully soon */
/*########################################
 ###  Calculation of shared fields s_Gc ###
 ########################################*/
void Species::AddSeed() {
    /* Stripes shared by several processors are redefined */
    for(int site=0;site<sites;site++) {
        if(p_rank){
            if(!s_Seed[site]) s_Seed[site] = s_Gc[2][site];
            if(s_Seed[site]>1)
                if(s_Gc[2][site]) s_Seed[site] = min(s_Seed[site],s_Gc[2][site]);
        }
        if(p_rank<size-1){
            if(!s_Seed[site]) s_Seed[site] = s_Gc[3][site];
            if(s_Seed[site]>1)
                if(s_Gc[3][site]) s_Seed[site] = min(s_Seed[site],s_Gc[3][site]);
        }
    }
}
#endif

/*############################
 ###   Species::DeathRate  ###
 #############################*/

/* Here we assume a light-dependent version of the mortality.
 basal is the minimal species death rate, depending on the species wood density.
 When PPFD is smaller than the light compensation point, mortality risk is increased.
 When NDD is defined, this death rate is increased by density-dependence effect that impact survival of trees with a dbh<10cm . */

/* v.2.2 Simplify function Species::DeathRate -- JChave */

/* Changed v.2.2, _NDD */
/* v.2.3.0 function has been renamed to avoid possible confusion downstream */

inline float Species::DeathRateNDD(float PPFD, float dbh, int nppneg, float ndd) {
    
    float dr=0;
    float basal=m*(1-s_wsg);
    float dd=deltaD*ndd*(1-2*dbh/s_dmax);
    
    dr=basal;
    if (nppneg > s_leaflifespan) {
        dr+=1.0/timestep;
    }
    if (dd > 0) {
        dr+=dd;
    }
    
    return dr*timestep;
}

inline float Species::DeathRate(float PPFD, float dbh, int nppneg) {
    
    float dr=0;
    float basal=maxf(m-m1*s_wsg,0.0);

    dr=basal;
    
    if (nppneg > s_leaflifespan) dr+=1.0/timestep;
    
    if (iter == int(nbiter-1))
        output[26]<< s_wsg << "\t" << basal << "\t"  << dbh << "\t"  << dr   <<  "\n";
    
    return dr*timestep;
}


/*#############################################
 ###   Farquhar von Caemmerer Berry model  ###
 ###           Species:: NPP               ###
 #############################################*/

/* This function returns the leaf-level carbon assimilation rate in micromoles C/m^2/s according to Farquhar-von Caemmerer-Berry model */
/* The function Species::dailyGPPleaf returns the primary productivity per unit leaf area, i.e. in micromoles C/m^2/s.
 It is converted into gC per m^2 of leaf per timestep by "*189.3*timestep" accounting only for the light hours (12 hours instead of 24): 189.3=12*3600*365.25*12/1000000
 BEWARE: 12 is the molar mass of C, and also the number of light hours in a day
 BEWARE: timestep is given as fraction of a year, so what is computed is actually the full assimilation per year which, in turn, is multiplied by the fraction per year that is under consideration
 BEWARE: slight inconsistency through use of 365.25 when daily timestep is likely to be given as 365, but not incorrect
 Commented version below was in use prior to version 2.3.0 -- use of lookup tables for acceleration of T dependence. cf. Bernacchi et al 2003 PCE; von Caemmerer 2000
 derivation of s_fci (ci/ca) according to Medlyn et al 2011, see also Prentice et al 2014 Ecology Letters and Lin et al 2015 Nature Climate Change --- initial version: s_fci = minf(-0.04*s_d13C-0.3*(log(PPFD)-s_factord13Cb)*s_factord13Ca-0.57, 1.0);
 from d13C (see cernusak et al 2013) without explicit model of stomatal conductance; min added in order to prevent ci:ca bigger than 1 (even though Ehleringer et al 1986 reported some values above 1 (Fig3) */

inline float Species::GPPleaf(float PPFD, float VPD, float T) {
    
    /* v.2.3.0: theta defined as a global variable */
    //theta=0.7;   // this is the fixed value of theta used by von Caemmerer 2000
    
    //float theta=0.76+0.018*T-0.00037*T*T;         // theta, but temperature dependent cf. Bernacchi et al 2003 PCE
    
    /* Parameters for Farquhar model, with temperature dependencies */
    int convT= int(iTaccuracy*T); // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---
    float KmT    =LookUp_KmT[convT];
    float GammaT =LookUp_GammaT[convT];
    tempRday    +=s_Rdark*LookUp_tempRday[convT];
    float VcmaxT =s_Vcmax*LookUp_VcmaxT[convT];
    float JmaxT  =s_Jmax*LookUp_JmaxT[convT];
    
    s_fci=g1/(g1+sqrt(VPD));
    
    /* FvCB model */
    
    float I=alpha*PPFD;
    float J = (I+JmaxT-sqrt((JmaxT+I)*(JmaxT+I)-4.0*theta*JmaxT*I))*0.5/theta;
    float A = minf(VcmaxT/(s_fci+KmT),0.25*J/(s_fci+2.0*GammaT))*(s_fci-GammaT);
    
    return A;
}

/* dailyGPPleaf returns the assimilation rate (computed from Species::GPPleaf) averaged across the daily fluctuations in climatic conditions (light, VPD and T), in micromoles C/m^2/s */

/* used only by _DAILYLIGHT */


inline float Species::dailyGPPleaf(float PPFD, float VPD, float T) {
    float ppfde,dailyA=0.0;
    
    for(int i=0; i<24; i++) {
        ppfde=PPFD*daily_light[i];
        if(ppfde > 0.1)
            // new v.2.3.0: compute GPP only if enough light is available threshold is arbitrary, but set to be low: in full sunlight ppfd is aroung 700 W/m2, and even at dawn, it is ca 3% of the max value, or 20 W/m2. The minimum threshold is set to 0.1 W/m2
            // Future update: compute slightly more efficiently, using 3-hourly values? This will have to be aligned with climate forcing layers (e.g. NCAR)
            dailyA+=GPPleaf(ppfde, VPD*daily_vpd[i], T*daily_T[i]);
        //the 6 lines in comment below corresponds to a finer version in which the multiplier is computed and used every 48 half hour, ie. with the corresponding environment instead of assuming a constant multiplier correponding the one at maximum incoming irradiance
        //float hhA=0;
        //hhA=GPPleaf(PPFD*daily_light[i], VPD*daily_vpd[i], T*daily_T[i]);
        //float alpha=phi*PPFD*daily_light[i]/hhA;
        //float D=klight*dens*CD;
        //hhA*=alpha/(D*(alpha-1))*log(alpha/(1+(alpha-1)*exp(-D)));
        //dailyA+=hhA;
    }
    //daily_light is the averaged (across one year, meteo station Nouragues DZ) and normalized (from 0 to 1) daily fluctuation of light, with half-hour time step, during the day time (from 7am to 7pm, ie 12 hours in total), same for daily_vpd and daily_T. Taking into account these daily variation is necessary considering the non-linearity of FvCB model
    
    dailyA*=0.0417;                                 // 0.0417=1/24 (24=12*2 = number of half hours in the 12 hours of daily light)
    tempRday*=0.0417;
    return dailyA;
}

/* v.2.4.0: formerly fastdailyGPPleaf(), since 2.3.2 directly included into dailyGPPleaf() */
/* but: structurally very different functions, hence proposition: renaming into dailyGPPcrown() and again separation*/
/* also: many function calls of dailyGPPleaf would result in many if(_FASTGPP) checks */

inline float Species::dailyGPPcrown(float PPFD, float VPD, float T, float dens, float CD) {
    float ppfde,dailyA=0.0;
    
    for(int i=0; i<24; i++) {
        ppfde=PPFD*daily_light[i];
        if(ppfde > 0.1)
            // new v.2.3.0: compute GPP only if enough light is available threshold is arbitrary, but set to be low: in full sunlight ppfd is aroung 700 W/m2, and even at dawn, it is ca 3% of the max value, or 20 W/m2. The minimum threshold is set to 0.1 W/m2
            // Future update: compute slightly more efficiently, using 3-hourly values? This will have to be aligned with climate forcing layers (e.g. NCAR)
            dailyA+=GPPleaf(ppfde, VPD*daily_vpd[i], T*daily_T[i]);
        //the 6 lines in comment below corresponds to a finer version in which the multiplier is computed and used every 48 half hour, ie. with the corresponding environment instead of assuming a constant multiplier correponding the one at maximum incoming irradiance
        //float hhA=0;
        //hhA=GPPleaf(PPFD*daily_light[i], VPD*daily_vpd[i], T*daily_T[i]);
        //float alpha=phi*PPFD*daily_light[i]/hhA;
        //float D=klight*dens*CD;
        //hhA*=alpha/(D*(alpha-1))*log(alpha/(1+(alpha-1)*exp(-D)));
        //dailyA+=hhA;
    }
    //daily_light is the averaged (across one year, meteo station Nouragues DZ) and normalized (from 0 to 1) daily fluctuation of light, with half-hour time step, during the day time (from 7am to 7pm, ie 12 hours in total), same for daily_vpd and daily_T. Taking into account these daily variation is necessary considering the non-linearity of FvCB model
    
    float alpha=phi*PPFD/GPPleaf(PPFD, VPD, T);             //alpha is a non-dimensional figure used to compute the multiplier below
    float D=klight*dens*CD;                                 //D is a non-dimensional figure used to compute the multiplier below
    dailyA*=alpha/(D*(alpha-1))*log(alpha/(1+(alpha-1)*exp(-D)));  // the FvCB assimilation rate computed at the top of the tree crown is multiplied by a multiplier<1, to account for the lower rate at lower light level within the crown depth. This multiplier is computed assuming that change in photosynthetic assimilation rate within a tree crown is mainly due to light decrease due to self-shading following a Michealis-menten relationship (ie. we assume that 1/ the change is not due to changes in VPD or temperature, which are supposed homogeneous at the intra-crown scale, and 2/ that other tree contributions to light decrease is neglected).
    
    dailyA*=0.0417;                                 // 0.0417=1/24 (24=12*2 = number of half hours in the 12 hours of daily light)
    tempRday*=0.0417;
    return dailyA;
}


/*############################################
 ############################################
 ############     Tree  class    ############
 ############################################
 ############################################*/

class Tree {

/* in v.2.4.0: t_C is not a tree variable anymore, but only calculated locally inside function */
//private:
//    float t_C;                    /* flexural force intensity, _TREEFALL, float? */
    
public:
    int   t_site,           /* location */
    t_NPPneg;               /* diagnostic variable: number of consecutive timesteps with NPP<0 -- V.2.2 */
    
    float t_dbh_thresh,       /* dbh threshold */
    t_hmax,                 /* allometric parameter, not real maximum */
    t_ah,                   /* new in v.2.4: allometric parameter, for consistency with t_hmax also an individual parameter */
    t_dbhmature,            /* reproductive size threshold IM jan2017 -- v230 */
    t_dbh,                  /* diameter at breast height (in m, but used in number of horizontal cells throughout all the code) */
    t_Tree_Height,          /* total tree height (in m, but used in number of vertical cells throughout all the code) */
    t_Crown_Depth,          /* crown depth (in m, but used in number of vertical cells throughout all the code) */
    t_Crown_Radius,         /* crown radius (in m, but used in number of horizontal cells throughout all the code)*/
    t_Ct,                   /* flexural force threshold, _BASICTREEFALL & t_TREEFALL */
    t_GPP,                  /* tree gross primary productivity */
    t_NPP,                  /* tree net primary productivity (in gC/timestep) */
    t_Rday,                 /* leaf respiration during day */
    t_Rnight,               /* leaf respiration during night */
    t_Rstem,                /* stem respiration */
    t_PPFD,                 /* light intensity received by the tree (computed by Tree::Flux, depending of the LAI at the tree height) */
    t_VPD,                    /* VPD at tree height -- v.2.2 */
    t_T,                    /* Temperature at tree height -- v.2.2 */
    t_ddbh,                 /* increment of dbh per timestep */
    t_age,                  /* tree age */
    t_youngLA,              /* total young leaf area, in m2 -- v.2.2  */
    t_matureLA,             /* total mature leaf area, in m2 -- v.2.2  */
    t_oldLA,                /* total old leaf area, in m2 -- v.2.2  */
    t_leafarea,             /* total crown leaf area in m2 -- v.2.2  */
    t_dens,                 /* tree crown average leaf density in m2/m2 -- v.2.2  */
    t_litter;               /* tree litterfall at each timestep, in g (of dry mass) -- v.2.2  */
    float *t_NDDfield;      /* _NDD */
 
#ifdef INTRASPECIFIC
    /* new in v.2.4.0: lognormal variation within species, hence multipliers */
    float
    t_intraspecific_multiplier_height,
    t_intraspecific_multiplier_CR,
    t_intraspecific_multiplier_P,
    t_intraspecific_multiplier_N,
    t_intraspecific_multiplier_LMA,
    t_intraspecific_multiplier_CD,
    t_intraspecific_multiplier_dmax,
    
    /* new in v.2.4.0: normal variation within species, hence absolute deviation */
    t_intraspecific_deviation_wsg,
    
    /* new in v.2.4.0: traits defined at the individual level */
    t_Pmass,                        /* phosphorus content, defined for individual instead of species */
    t_Nmass,                        /* nitrogen content, defined for individual instead of species */
    t_LMA,                          /* LMA, defined for individual instead of species  */
    t_wsg,                          /* wood specific gravity, defined for individual instead of species */
    t_dmax,                         /* maximum diameter, estimated from field data */

    t_Rdark,                        /* dark respiration rate (at PPFD = 0) in micromol C/m^2/s),  */
    t_Vcmax,                        /* maximal rate of carboxylation, on an area basis (in micromolC/m^2/s) */
    t_Vcmaxm,                       /* maximal rate of carboxylation, on an mass basis (in micromolC/g-1/s) */
    t_Jmax,                         /* maximal rate of electron transport, on an area basis (in micromol/m^2/s) */
    t_Jmaxm,                        /* maximal rate of electron transport, on a mass basis (in micromol/g-1/s) */
    t_fci,                          /* fraction of CO2 partial pressure in intercellular spaces divided by ambiant CO2 partial pressure (both in microbar, or ppm = micromol/mol) */
    t_Gamma,                /* compensation point for the carboxylation rate, here NORMALIZED by atm CO2 concentration (Cair) */
    t_Km,                   /* apparent kinetic constant for the rubiscco = Kc*(1+[O]/Ko), here normalized by atm CO2 concentration (Cair) NINO : useless ?*/
    
    t_leaflifespan,         /* average leaf lifespan, in month */
    t_time_young,           /* leaf resident time in the young leaf class */ // IF not intra ???
    t_time_mature,          /* leaf resident time in the mature leaf class */
    t_time_old;             /* leaf resident time in the old leaf class */
#endif
    
#ifdef CROWN_SHAPE
    float
    t_Crown_Slope_Top,
    t_Crown_Slope_Bottom,
    t_Crown_Volume,
    t_Crown_Volume_layer;
#endif
    
    Species *t_s;
    
    unsigned short
    t_from_Data,            /* indicator of whether tree is born through initialisation or through simulation routine */
    t_sp_lab,               /* species label */
    t_hurt;                 /* treefall index */
    
    Tree(){                 /* constructor */
        t_from_Data = 0;
        t_sp_lab = 0;
        t_age = 0;
        t_hurt = 0;
        t_NPP=t_GPP=t_Rday=t_Rnight=t_Rstem=t_PPFD=t_VPD=t_T=0.0; /* new v.2.2 */
        t_dbh = t_Tree_Height = t_Crown_Radius = t_Crown_Depth= 0.0;
        
        if(_BASICTREEFALL || _TREEFALL) t_Ct = 0.0;
        
    };
    
    virtual ~Tree() {
        delete [] t_NDDfield;   /* _NDD */
    };	/* destructor */
    
    void Birth(Species*,int,int);	/* tree birth */
    void BirthFromData(Species *S, int nume, int site0, float dbh_measured); /* tree initialisation from field data */
    void BirthFromSimulation(Species *S, int nume, int site0, float dbh_measured, float parameters[]); /* tre initialization from FullFinal*/
    void Death();                   /* tree death */
    void Growth();                  /* tree growth */
    void Fluxh(int h);             /* compute mean light flux received by the tree crown layer at height h new version (PPFD symmetrical with T and VPD) -- v230 */
    void CalcRespGPP();             /* new in v.2.4: calculates respiration and GPP across crown */
    void UpdateLeafDynamics();      /* computes within-canopy leaf dynamics (IM since v 2.1, as a standalone function in v.2.3.0) */
    void UpdateTreeBiometry();      /* compute biometric relations, including allometry -- contains a large part of empirical functions */

    void DisperseSeed();            /* update Seed field */
    void Couple(float &c_forceflex, float &angle); /* force exerted by other trees, _TREEFALL, changed in v.2.4.0: c_forceflex and angle are now defined locally and passed to Couple() by reference for clarity's sake */
    void Treefall(float angle);     /* formerly FallTree(), tree falling routine, _BASICTREEFALL & _TREEFALL, changed in v.2.4.0, only computes the actual treefall, but not the probabilities of it happening any more. Takes angle as argument and can now be used universally, be it for classic treefalls, secondary treefalls, treefelling, or other disturbances (cyclones, etc.) */
    void FellTree();				/* tree felling routine, _LOGGING NINO : To BE ADAPTED */

    void Update();                  /* tree evolution at each timestep */
    void Average();                 /* local computation of the averages */
    void CalcLAI();                 /* updating the LAI3D field */
    void histdbh();                 /* computation of dbh histograms */
    void OutputTreeStandard();  /* creates standard output for trees, writes directly to cout stream */
    void OutputTreeStandard(fstream& output);         /* overloading of function, creates standard output for trees and allows user to add 5 variables in floating point format (either global variables or locally defined ones), then writes to file specified in argument */
#ifdef INTRASPECIFIC
    /* GPP functions are now calculated at tree level */
    float DeathRateNDD(float, float, int, float); /* _NDD, overloading with function in following line */
    float DeathRate(float, float, int);  /* actual death rate -- new argument int v.2.2 */
    float GPPleaf(float, float, float);    /* Computation of the light-limited leaf-level NPP per m^2 (in micromol/m^2/s) -- two new arguments float v.2.2 */
    /* Farquhar von Caemmerer Berry model */
    float dailyGPPleaf(float, float, float);    /* computation of the daily average assimilation rate, taking into account the daily variation in light, VPD and temperature two new arguments float v.2.2, _DAILYLIGHT */
    float dailyGPPcrown(float, float, float, float, float);  /* v.2.4.0: proposition of reinsertion of fastdailyGPPleaf() as dailyGPPcrown(), main reason: despite sharing some code with dailyGPPleaf, structurally very different */
#endif
};


/*##############################################
 ####	            Tree birth              ####
 ####  called by BirthInit and UpdateTree   ####
 ##############################################*/


void Tree::Birth(Species *S, int nume, int site0) {
    t_site = site0;
    t_sp_lab = nume;            /* t_sp_lab is the species label of a site. Can be defined even if the site is empty (cf. persistence function defined in Chave, Am Nat. 2001) */
    t_NPPneg=0.0;
    t_s = S+t_sp_lab;
    t_age = 1;
    t_hurt = 0;
    t_Tree_Height = H0;
    t_hmax = (t_s->s_hmax);
    t_ah = (t_s->s_ah);

#ifdef INTRASPECIFIC
    int dev_rand = genrand2i()%100000;
    t_intraspecific_multiplier_height = d_intraspecific_height[dev_rand];
    t_intraspecific_multiplier_CR = d_intraspecific_CR[dev_rand];
    t_intraspecific_multiplier_N = d_intraspecific_N[dev_rand];
    t_intraspecific_multiplier_P = d_intraspecific_P[dev_rand];
    t_intraspecific_multiplier_LMA = d_intraspecific_LMA[dev_rand];
    t_intraspecific_multiplier_CD = d_intraspecific_CD[dev_rand];
    t_intraspecific_deviation_wsg = d_intraspecific_wsg[dev_rand];
    t_intraspecific_multiplier_dmax = d_intraspecific_dmax[dev_rand];

    t_Pmass = (t_s->s_Pmass) * t_intraspecific_multiplier_P;
    t_Nmass = (t_s->s_Nmass) * t_intraspecific_multiplier_N;
    t_LMA = (t_s->s_LMA) * t_intraspecific_multiplier_LMA;
    t_wsg = maxf((t_s->s_wsg) + t_intraspecific_deviation_wsg,0.1); // cutoff for normal distribution (could also be set at 0.0)
    t_dmax = (t_s->s_dmax) * t_intraspecific_multiplier_dmax;
#ifdef ALLOM_relative
    t_ah *= t_dmax; /* if the Michaelis Menten parameters are inferred from a relative diameter (t_dbh_rel = t_dbh/t_dmax) instead of the unscaled parameter t_dbh, then the equation h_rel = t_hmax * t_dbh_rel/ (t_dbh_rel + t_ah) is equivalent to t_hmax * t_dbh / (t_dbh + t_ah_rescaled), where ah_rescaled = t_ah * t_dmax. This means that, apart from rescaling of t_ah, everything else can be calculated as usually */
#endif
    t_dbh=t_ah*H0/(t_hmax*t_intraspecific_multiplier_height-H0);    // we make sure that all trees start at the same height, this means adapting the diameter depending on inter- and intraspecific changes

    float SLA=10000.0/t_LMA;
    
    t_Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(t_Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(t_Pmass*1000.0)+0.25*log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
    t_Jmaxm=pow(10.0, minf((-1.50+0.41*log10(t_Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(t_Pmass*1000.0)+0.32*log10(SLA)))); // added as a Species member variable 14-04-2015; this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
    t_Vcmax=t_Vcmaxm*t_LMA;
    t_Jmax=t_Jmaxm*t_LMA;
    t_Rdark=t_LMA*(8.5341-130.6*t_Nmass-567.0*t_Pmass-0.0137*t_LMA+11.1*t_Vcmaxm+187600.0*t_Nmass*t_Pmass)*0.001; //t_Rdark corresponds to leaf maintenance respiration. From Table 6 in Atkin et al 2015 New phytologist v.2.0
    
    t_Gamma = 38.0*iCair; // s_Gamma at 25°C computed according to von Caemmerer 2000 formula: gamma=Kc*O*0.25/(2*Ko), with Kc=260 microbar, Ko=179mbar and O=210 mbar (the last value is from Farquhar et al 1980, the first two one are from von Caemmerer 2000 table 2.3 page 45). gamma is set to 36.9 on Atkin et al 2015. Could be a global variable. v.2.0
    
    t_leaflifespan = pow(10,(2.040816*(2.579713-log10(SLA)))); //this is the expression from Reich et al. 1997 PNAS (provides probably more realistic estimates for species with high LMA).
    //t_leaflifespan=1.5+pow(10,(7.18+3.03*log10(t_LMA*0.0001)));           //this is the expression from Reich et al 1991 Oecologia (San Carlos Rio Negro).
    //t_leaflifespan=0.5+pow(10,(-2.509+1.71*log10(t_LMA)));    //this is the expression from Wright et al 2004 Nature (leaf economics spectrum).

#ifdef LL_limit
    t_leaflifespan = maxf(t_leaflifespan,3.0);
#endif
    
    t_time_young=minf(t_leaflifespan/3.0,1.0);
    t_time_mature=t_leaflifespan/3.0;
    t_time_old=t_leaflifespan-t_time_mature-t_time_young;
    t_litter = 0.0;
    t_fci = 0.0;

    t_Crown_Radius = exp(1.9472 + 0.5925*log(t_dbh));   /* this is crown allometry derived from data set provided by Jucker et al. 2016 (Global Change Biology) */

    t_Crown_Radius  *= t_intraspecific_multiplier_CR;
    t_Crown_Depth = de0 * t_intraspecific_multiplier_CD;
    
    t_ddbh=0.0;
    t_dbh_thresh = t_dmax;  // contrary to version without intraspecific variation, dmax already includes variation, so t_dbh_thresh = t_dmax
    t_dbhmature=t_dbh_thresh*0.5; // this correponds to the mean thresholds of tree size to maturity, according to Visser et al. 2016 Functional Ecology (suited to both understory short-statured species, and top canopy large-statured species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230
    //float u=genrand2();
    //t_dbhmature=maxf(0, -(t_dmax)*0.25*log((1-u)/u)+t_dmax*0.5);    // IM test 02-2017, try to introduce intra-species inter-individual variability in dbhmature, following a sigmoidal repartition function, as in Visser et al. 2016 and Wright et al. 2005
    float hrealmax = t_intraspecific_multiplier_height * 1.5*t_hmax*t_dbh_thresh/(1.5*t_dbh_thresh+t_ah);
    
#else
    //t_dbh=(t_s->s_dbh0);
#ifdef ALLOM_relative
    t_ah *= t_s->s_dmax; /* if the Michaelis Menten parameters are inferred from a relative diameter (t_dbh_rel = t_dbh/t_dmax) instead of the unscaled parameter t_dbh, then the equation h_rel = t_hmax * t_dbh_rel/ (t_dbh_rel + t_ah) is equivalent to t_hmax * t_dbh / (t_dbh + t_ah_rescaled), where ah_rescaled = t_ah * t_dmax. This means that, apart from rescaling of t_ah, everything else can be calculated as usually */
#endif
    t_dbh=t_ah*H0/(t_hmax-H0);    // we make sure that all trees start at the same height, this means adapting the diameter depending on inter- and intraspecific changes
    t_Crown_Radius  = exp(1.9472 + 0.5925*log(t_dbh));   /* this is crown allometry derived from data set provided by Jucker et al. 2016 (Global Change Biology) */
    t_Crown_Depth = de0;
    t_ddbh=0.0;
    t_dbh_thresh = ((t_s->s_dmax)-t_dbh)*flor(1.0+log(genrand2())*0.01)+t_dbh;
    t_dbhmature=t_dbh_thresh*0.5; // this correponds to the mean thresholds of tree size to maturity, according to Visser et al. 2016 Functional Ecology (suited to both understory short-statured species, and top canopy large-statured species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230
    //float u=genrand2();
    //t_dbhmature=maxf(0, -(t_s->s_dmax)*0.25*log((1-u)/u)+t_s->s_dmax*0.5);    // IM test 02-2017, try to introduce intra-species inter-individual variability in dbhmature, following a sigmoidal repartition function, as in Visser et al. 2016 and Wright et al. 2005
    float hrealmax=1.5*t_hmax*t_dbh_thresh/(1.5*t_dbh_thresh+t_ah);
#endif
    
    if(_BASICTREEFALL || _TREEFALL) t_Ct = hrealmax*flor(1.0-vC*sqrt(-log(genrand2())));
    
    t_dens=dens;

#ifdef CROWN_SHAPE
    t_Crown_Volume_layer = 0.0;
    
    int crown_base = int(t_Tree_Height - t_Crown_Depth);
    int crown_top = int(t_Tree_Height);
    int crown_center = int(t_Tree_Height-0.5*t_Crown_Depth);
    float crown_extent_top = t_Tree_Height - crown_center;
    float crown_extent_base = crown_center - (t_Tree_Height - t_Crown_Depth);
    
    if(t_Crown_Depth > 2.0){        // only consider cases, where more than two layers are reached by the crown
        /* slopes are calculated as horizontal/vertical slopes (inverse of more intuitive definition vertical/horizontal), reasoning: when shape_crown = 1.0, this slope is still defined */
        /* the slope at the base part of the crown is steeper, meaning that the trapezoid narrows down faster (limit: triangle/cone) */
        t_Crown_Slope_Top = t_Crown_Radius * (1.0 - shape_crown) / crown_extent_top;
        t_Crown_Slope_Bottom = min(2 * t_Crown_Slope_Top, t_Crown_Radius/crown_extent_base);
    } else {
        t_Crown_Slope_Top = 0.0;
        t_Crown_Slope_Bottom = 0.0;
    }
    t_Crown_Volume = 0.0;
    for(int h=crown_base;h<crown_top+1;h++){
        /* calculating the height of the current layer */
        float height_layer = 1.0;
        if(crown_top == crown_base) height_layer = t_Crown_Depth;
        else if(h == crown_top) height_layer = (t_Tree_Height-crown_top);
        else if(h == crown_base) height_layer = (crown_base+1-(t_Tree_Height-t_Crown_Depth));
        
        /* calculating the radius of the current layer depending on the respective slopes */
        float radius_layer = t_Crown_Radius;
        if(crown_top == crown_base){}
        else if(h == crown_base) radius_layer -= t_Crown_Slope_Bottom * crown_extent_base;
        else if(h < crown_center ) radius_layer -= t_Crown_Slope_Bottom * (crown_center-h);
        else radius_layer -= t_Crown_Slope_Top * (h - crown_center);            // for h = crown_center full radius
        
        float crown_area_fl = maxf(0.0,PI * radius_layer * radius_layer);
        t_Crown_Volume += crown_area_fl * height_layer;
    }

    t_youngLA=t_dens*t_Crown_Volume;
#else
    t_youngLA=t_dens*PI*t_Crown_Radius*LH*t_Crown_Radius*LH*t_Crown_Depth*LV;
#endif
    /* initially, all stems have only young leaves -- LA stands for leaf area */
    t_matureLA=0.0;           /* this is the amount of leaf area at maturity */
    t_oldLA=0.0;              /* leaf area of senescing leaves */
    t_leafarea=t_youngLA;   /* should be sum of LA young+mature+old, but the equation is correct initially */
    tempRday=0.0;

    (t_s->s_nbind)++;
    nblivetrees++;
}


/*##############################################
 ####	            Tree birth from Data    ####
 ####  called by BirthInit and UpdateTree   ####
 ##############################################*/


void Tree::BirthFromData(Species *S, int nume, int site0, float dbh_measured) {
    t_site = site0;
    t_sp_lab = nume;            /* t_sp_lab is the species label of a site. Can be defined even if the site is empty (cf. persistence function defined in Chave, Am Nat. 2001) */
    t_NPPneg=0.0;
    t_s = S+t_sp_lab;
    t_age = 1;
    t_from_Data = 1;
    t_hurt = 0;
    //t_Tree_Height = H0; // change
    t_hmax = (t_s->s_hmax);
    t_ah = (t_s->s_ah);

    if((t_s->s_dmax)*1.5 > dbh_measured) t_dbh = dbh_measured;  // UNIT ? NINO   // force dbh to be within limits of TROLL specifications
    else {
        t_dbh = (t_s->s_dmax); // Yes but, is intraspecific variability to be taken into account in this part ?  Cause it impacts the individual dbh threshold.
        cout << "Warning: DBH_measured > 1.5*DBH_max for species. DBH set to DBH_max for species \n";
    }

#ifdef INTRASPECIFIC
    int dev_rand = genrand2i()%100000;
    t_intraspecific_multiplier_height = d_intraspecific_height[dev_rand];
    t_intraspecific_multiplier_CR = d_intraspecific_CR[dev_rand];
    t_intraspecific_multiplier_N = d_intraspecific_N[dev_rand];
    t_intraspecific_multiplier_P = d_intraspecific_P[dev_rand];
    t_intraspecific_multiplier_LMA = d_intraspecific_LMA[dev_rand];
    t_intraspecific_multiplier_CD = d_intraspecific_CD[dev_rand];
    t_intraspecific_deviation_wsg = d_intraspecific_wsg[dev_rand];
    t_intraspecific_multiplier_dmax = d_intraspecific_dmax[dev_rand];

    t_Pmass = (t_s->s_Pmass) * t_intraspecific_multiplier_P;
    t_Nmass = (t_s->s_Nmass) * t_intraspecific_multiplier_N;
    t_LMA = (t_s->s_LMA) * t_intraspecific_multiplier_LMA;
    t_wsg = maxf((t_s->s_wsg) + t_intraspecific_deviation_wsg,0.1); // cutoff for normal distribution (could also be set at 0.0)
    t_dmax = (t_s->s_dmax) * t_intraspecific_multiplier_dmax;

#ifdef ALLOM_relative
    t_ah *= t_dmax; /* if the Michaelis Menten parameters are inferred from a relative diameter (t_dbh_rel = t_dbh/t_dmax) instead of the unscaled parameter t_dbh, then the equation h_rel = t_hmax * t_dbh_rel/ (t_dbh_rel + t_ah) is equivalent to t_hmax * t_dbh / (t_dbh + t_ah_rescaled), where ah_rescaled = t_ah * t_dmax. This means that, apart from rescaling of t_ah, everything else can be calculated as usually */
#endif
   // t_dbh=t_ah*H0/(t_hmax*t_intraspecific_multiplier_height-H0);    // change  we make sure that all trees start at the same height, this means adapting the diameter depending on inter- and intraspecific changes

    float SLA=10000.0/t_LMA;
    
    t_Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(t_Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(t_Pmass*1000.0)+0.25*log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
    t_Jmaxm=pow(10.0, minf((-1.50+0.41*log10(t_Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(t_Pmass*1000.0)+0.32*log10(SLA)))); // added as a Species member variable 14-04-2015; this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
    t_Vcmax=t_Vcmaxm*t_LMA;
    t_Jmax=t_Jmaxm*t_LMA;
    t_Rdark=t_LMA*(8.5341-130.6*t_Nmass-567.0*t_Pmass-0.0137*t_LMA+11.1*t_Vcmaxm+187600.0*t_Nmass*t_Pmass)*0.001; //t_Rdark corresponds to leaf maintenance respiration. From Table 6 in Atkin et al 2015 New phytologist v.2.0
    
    t_Gamma = 38.0*iCair; // s_Gamma at 25°C computed according to von Caemmerer 2000 formula: gamma=Kc*O*0.25/(2*Ko), with Kc=260 microbar, Ko=179mbar and O=210 mbar (the last value is from Farquhar et al 1980, the first two one are from von Caemmerer 2000 table 2.3 page 45). gamma is set to 36.9 on Atkin et al 2015. Could be a global variable. v.2.0
    
    t_leaflifespan = pow(10,(2.040816*(2.579713-log10(SLA)))); //this is the expression from Reich et al. 1997 PNAS (provides probably more realistic estimates for species with high LMA).
    //t_leaflifespan=1.5+pow(10,(7.18+3.03*log10(t_LMA*0.0001)));           //this is the expression from Reich et al 1991 Oecologia (San Carlos Rio Negro).
    //t_leaflifespan=0.5+pow(10,(-2.509+1.71*log10(t_LMA)));    //this is the expression from Wright et al 2004 Nature (leaf economics spectrum).

#ifdef LL_limit
    t_leaflifespan = maxf(t_leaflifespan,3.0);
#endif
    
    t_time_young=minf(t_leaflifespan/3.0,1.0);
    t_time_mature=t_leaflifespan/3.0;
    t_time_old=t_leaflifespan-t_time_mature-t_time_young;
    t_litter = 0.0;
    t_fci = 0.0;

//    t_Crown_Radius = exp(1.9472 + 0.5925*log(t_dbh));   /* this is crown allometry derived from data set provided by Jucker et al. 2016 (Global Change Biology) */

//    t_Crown_Radius  *= t_intraspecific_multiplier_CR;
//   t_Crown_Depth = de0 * t_intraspecific_multiplier_CD;
    
    t_dbh_thresh = t_dmax;  // contrary to version without intraspecific variation, dmax already includes variation, so t_dbh_thresh = t_dmax
 //   t_dbhmature=t_dbh_thresh*0.5; // this correponds to the mean thresholds of tree size to maturity, according to Visser et al. 2016 Functional Ecology (suited to both understory short-statured species, and top canopy large-statured species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230
    //float u=genrand2();
    //t_dbhmature=maxf(0, -(t_dmax)*0.25*log((1-u)/u)+t_dmax*0.5);    // IM test 02-2017, try to introduce intra-species inter-individual variability in dbhmature, following a sigmoidal repartition function, as in Visser et al. 2016 and Wright et al. 2005
    float hrealmax = t_intraspecific_multiplier_height * 1.5*t_hmax*t_dbh_thresh/(1.5*t_dbh_thresh+t_ah);
    
#else
    //t_dbh=(t_s->s_dbh0);
#ifdef ALLOM_relative
    t_ah *= t_s->s_dmax; /* if the Michaelis Menten parameters are inferred from a relative diameter (t_dbh_rel = t_dbh/t_dmax) instead of the unscaled parameter t_dbh, then the equation h_rel = t_hmax * t_dbh_rel/ (t_dbh_rel + t_ah) is equivalent to t_hmax * t_dbh / (t_dbh + t_ah_rescaled), where ah_rescaled = t_ah * t_dmax. This means that, apart from rescaling of t_ah, everything else can be calculated as usually */
#endif
    //t_dbh=t_ah*H0/(t_hmax-H0);    // we make sure that all trees start at the same height, this means adapting the diameter depending on inter- and intraspecific changes

    t_dbh_thresh = ((t_s->s_dmax)-t_dbh)*flor(1.0+log(genrand2())*0.01)+t_dbh;
 //   t_dbhmature=t_dbh_thresh*0.5; // this correponds to the mean thresholds of tree size to maturity, according to Visser et al. 2016 Functional Ecology (suited to both understory short-statured species, and top canopy large-statured species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230
    //float u=genrand2();
    //t_dbhmature=maxf(0, -(t_s->s_dmax)*0.25*log((1-u)/u)+t_s->s_dmax*0.5);    // IM test 02-2017, try to introduce intra-species inter-individual variability in dbhmature, following a sigmoidal repartition function, as in Visser et al. 2016 and Wright et al. 2005
    float hrealmax=1.5*t_hmax*t_dbh_thresh/(1.5*t_dbh_thresh+t_ah);
#endif
    
    if(_BASICTREEFALL || _TREEFALL) t_Ct = hrealmax*flor(1.0-vC*sqrt(-log(genrand2())));
    
    
    t_ddbh=0.0;
    t_dbhmature=t_dbh_thresh*0.5; // this correponds to the mean thresholds of tree size to maturity, according to Visser et al. 2016 Functional Ecology (suited to both understory short-statured species, and top canopy large-statured species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230


#ifdef INTRASPECIFIC// NINO
    t_Tree_Height = minf(t_intraspecific_multiplier_height * t_hmax * t_dbh/(t_dbh + t_ah),HEIGHT-1); // todef NINO-OK
#else
    t_Tree_Height = t_hmax * t_dbh/(t_dbh + t_ah);
#endif

t_Crown_Radius = exp(1.9472 + 0.5925*log(t_dbh)); /* this is crown allometry derived from data set provided by Jucker et al. 2016 (Global Change Biology) */
if(t_Tree_Height < 5.0) {t_Crown_Depth = de0 + 0.17 * (t_Tree_Height-H0);}
else {t_Crown_Depth = de0+0.26*(t_Tree_Height-H0) - 0.09 * (5.0-H0);} /* allometry deduced from Piste Saint-Elie dataset (unpublished) */
    
#ifdef INTRASPECIFIC
    t_Crown_Radius *= t_intraspecific_multiplier_CR;
    t_Crown_Depth *= t_intraspecific_multiplier_CD;
#endif

    t_dens=dens;
    t_leafarea=t_dens*PI*t_Crown_Radius*LH*t_Crown_Radius*LH*t_Crown_Depth;
    t_youngLA=0.25*t_leafarea;
    t_matureLA=0.5*t_leafarea;
    t_oldLA=0.25*t_leafarea;
    Fluxh(int(t_Tree_Height)+1);
    tempRday=0.0;

// useless in my case, so unchecked
#ifdef CROWN_SHAPE
    t_Crown_Volume_layer = 0.0;
    
    int crown_base = int(t_Tree_Height - t_Crown_Depth);
    int crown_top = int(t_Tree_Height);
    int crown_center = int(t_Tree_Height-0.5*t_Crown_Depth);
    float crown_extent_top = t_Tree_Height - crown_center;
    float crown_extent_base = crown_center - (t_Tree_Height - t_Crown_Depth);
    
    if(t_Crown_Depth > 2.0){        // only consider cases, where more than two layers are reached by the crown
        /* slopes are calculated as horizontal/vertical slopes (inverse of more intuitive definition vertical/horizontal), reasoning: when shape_crown = 1.0, this slope is still defined */
        /* the slope at the base part of the crown is steeper, meaning that the trapezoid narrows down faster (limit: triangle/cone) */
        t_Crown_Slope_Top = t_Crown_Radius * (1.0 - shape_crown) / crown_extent_top;
        t_Crown_Slope_Bottom = min(2 * t_Crown_Slope_Top, t_Crown_Radius/crown_extent_base);
    } else {
        t_Crown_Slope_Top = 0.0;
        t_Crown_Slope_Bottom = 0.0;
    }
    t_Crown_Volume = 0.0;
    for(int h=crown_base;h<crown_top+1;h++){
        /* calculating the height of the current layer */
        float height_layer = 1.0;
        if(crown_top == crown_base) height_layer = t_Crown_Depth;
        else if(h == crown_top) height_layer = (t_Tree_Height-crown_top);
        else if(h == crown_base) height_layer = (crown_base+1-(t_Tree_Height-t_Crown_Depth));
        
        /* calculating the radius of the current layer depending on the respective slopes */
        float radius_layer = t_Crown_Radius;
        if(crown_top == crown_base){}
        else if(h == crown_base) radius_layer -= t_Crown_Slope_Bottom * crown_extent_base;
        else if(h < crown_center ) radius_layer -= t_Crown_Slope_Bottom * (crown_center-h);
        else radius_layer -= t_Crown_Slope_Top * (h - crown_center);            // for h = crown_center full radius
        
        float crown_area_fl = maxf(0.0,PI * radius_layer * radius_layer);
        t_Crown_Volume += crown_area_fl * height_layer;
    }

#endif



    (t_s->s_nbind)++;
    nblivetrees++;
}

/*####################################################
#### Tree Initialization from Previous Simulation ####
####    New function added by Nino Page - 2017    ####
####################################################*/

void Tree::BirthFromSimulation(Species *S, int nume, int site0, float dbh_measured, float parameters[]) {
    t_site = site0;
    t_sp_lab = nume;            /* t_sp_lab is the species label of a site. Can be defined even if the site is empty (cf. persistence function defined in Chave, Am Nat. 2001) */
    t_s = S+t_sp_lab;
    t_dbh = dbh_measured;

    t_NPPneg= int(parameters[1]);
    t_age = parameters[2];
    t_hurt = parameters[3];

    t_dbh_thresh = parameters[4];
    t_hmax = parameters[5];
    t_ah = parameters[6]; // IFDEF alom_rel ???
    

    t_youngLA=parameters[7];			/* LA stands for leaf area */
    t_matureLA= parameters[8];           /* this is the amount of leaf area at maturity */
    t_oldLA=parameters[9];              /* leaf area of senescing leaves */
    t_leafarea= parameters[10];   /* should be sum of LA young+mature+old, but the equation is correct initially */
    t_dens= parameters[11]; // LAI parameter
    t_litter = parameters[12]; // 


//    if((t_s->s_dmax)*1.5 > dbh_measured) t_dbh = dbh_measured;  // UNIT ? NINO   // force dbh to be within limits of TROLL specifications
//    else {
//        t_dbh = (t_s->s_dmax); // Yes but, is intraspecific variability to be taken into account in this part ?  Cause it impacts the individual dbh threshold.
//        cout << "Warning: DBH_measured > 1.5*DBH_max for species. DBH set to DBH_max for species \n";
//    }

    t_Tree_Height = parameters[13];
    t_Crown_Radius = parameters[14];
    t_Crown_Depth = parameters[15];
    t_ddbh= parameters[16];

    if(_BASICTREEFALL || _TREEFALL) t_Ct = parameters[17];

    #ifdef INTRASPECIFIC
    t_intraspecific_multiplier_height = parameters[18];
    t_intraspecific_multiplier_CR = parameters[19];
    t_intraspecific_multiplier_CD = parameters[20];
    t_intraspecific_multiplier_N = parameters[21];
    t_intraspecific_multiplier_P = parameters[22];
    t_intraspecific_multiplier_LMA = parameters[23];
    t_intraspecific_deviation_wsg = parameters[24];
    t_intraspecific_multiplier_dmax = parameters[25];

    t_Pmass = (t_s->s_Pmass) * t_intraspecific_multiplier_P;
    t_Nmass = (t_s->s_Nmass) * t_intraspecific_multiplier_N;
    t_LMA = (t_s->s_LMA) * t_intraspecific_multiplier_LMA;

    t_wsg = maxf((t_s->s_wsg) + t_intraspecific_deviation_wsg,0.1); // cutoff for normal distribution (could also be set at 0.0)
    t_dmax = (t_s->s_dmax) * t_intraspecific_multiplier_dmax;
    // t_dmax = parameters[];


    float SLA=10000.0/t_LMA;
    
    t_Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(t_Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(t_Pmass*1000.0)+0.25*log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
    t_Jmaxm=pow(10.0, minf((-1.50+0.41*log10(t_Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(t_Pmass*1000.0)+0.32*log10(SLA)))); // added as a Species member variable 14-04-2015; this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
    t_Vcmax=t_Vcmaxm*t_LMA;
    t_Jmax=t_Jmaxm*t_LMA;
    t_Rdark=t_LMA*(8.5341-130.6*t_Nmass-567.0*t_Pmass-0.0137*t_LMA+11.1*t_Vcmaxm+187600.0*t_Nmass*t_Pmass)*0.001; //t_Rdark corresponds to leaf maintenance respiration. From Table 6 in Atkin et al 2015 New phytologist v.2.0
    t_Gamma = 38.0*iCair; // s_Gamma at 25°C computed according to von Caemmerer 2000 formula: gamma=Kc*O*0.25/(2*Ko), with Kc=260 microbar, Ko=179mbar and O=210 mbar (the last value is from Farquhar et al 1980, the first two one are from von Caemmerer 2000 table 2.3 page 45). gamma is set to 36.9 on Atkin et al 2015. Could be a global variable. v.2.0
    t_leaflifespan = pow(10,(2.040816*(2.579713-log10(SLA)))); //this is the expression from Reich et al. 1997 PNAS (provides probably more realistic estimates for species with high LMA).

        // Easily computed
    t_time_young=minf(t_leaflifespan/3.0,1.0);
    t_time_mature=t_leaflifespan/3.0;
    t_time_old=t_leaflifespan-t_time_mature-t_time_young;
    t_fci = 0.0; // parameters[]?

    #ifdef LL_limit
    t_leaflifespan = maxf(t_leaflifespan,3.0);
    #endif
    
    // Easily computed
    t_time_young=minf(t_leaflifespan/3.0,1.0);
    t_time_mature=t_leaflifespan/3.0;
    t_time_old=t_leaflifespan-t_time_mature-t_time_young;
    t_fci = 0.0;
	#endif

	//#ifdef ALLOM_relative
   	//	t_ah *= t_dmax; /* if the Michaelis Menten parameters are inferred from a relative diameter (t_dbh_rel = t_dbh/t_dmax) instead of the unscaled parameter t_dbh, then the equation h_rel = t_hmax * t_dbh_rel/ (t_dbh_rel + t_ah) is equivalent to t_hmax * t_dbh / (t_dbh + t_ah_rescaled), where ah_rescaled = t_ah * t_dmax. This means that, apart from rescaling of t_ah, everything else can be calculated as usually */
	//#endif

    t_dbhmature=t_dbh_thresh*0.5;

	#ifdef CROWN_SHAPE
    	t_Crown_Slope_Top = parameters[26];
    	t_Crown_Slope_Bottom = parameters[27];
    	t_Crown_Volume = parameters[28];
    	t_Crown_Volume_layer = parameters[29];
    #endif

    tempRday=0.0;
    (t_s->s_nbind)++;
    nblivetrees++;
}

/*################################################
 #### Contribution of each tree in LAI field  ####
 ####          called by UpdateField          ####
 #################################################*/

/* modified in v.2.3: additional contribution to voxels that are not fully occupied by the tree crown. !!!: this does not calculate LAI3D directly, this only calculates the density in each voxel belonging to a tree. The final LAI field is calculated outside of the class Tree */
/* modified in v.2.4.0: version with CROWN_EXPANSION, making growth more continuous by gradual allocation due to distance from center of crown */

#ifdef CROWN_EXPANSION
#ifdef CROWN_SHAPE
void Tree::CalcLAI() {
    if(t_age > 0){
        int crown_base = int(t_Tree_Height-t_Crown_Depth),
        crown_top = int(t_Tree_Height),
        crown_center = int(t_Tree_Height-0.5*t_Crown_Depth),
        row_trunk = t_site/cols,
        col_trunk = t_site%cols,
        col,row,site;

        float crown_extent_base = crown_center - (t_Tree_Height - t_Crown_Depth);

        for(int h=crown_base;h < crown_top+1;h++){
            /* calculating the height of the current layer */
            float height_layer = 1.0;
            if(crown_top == crown_base) height_layer = t_Crown_Depth;
            else if(h == crown_top) height_layer = (t_Tree_Height-crown_top);
            else if(h == crown_base) height_layer = (crown_base+1-(t_Tree_Height-t_Crown_Depth));
            
             /* calculating the radius of the current layer depending on the respective slopes */
            float radius_layer = t_Crown_Radius;
            if(crown_top == crown_base){}
            else if(h == crown_base) radius_layer -= t_Crown_Slope_Bottom * crown_extent_base;
            else if(h < crown_center ) radius_layer -= t_Crown_Slope_Bottom * (crown_center-h);
            else radius_layer -= t_Crown_Slope_Top * (h - crown_center);            // for h = crown_center full radius
            
            int crown_area = int(PI * radius_layer * radius_layer);     // floor of crown_area to bound area accumulation
            crown_area = max(crown_area,1);                                 // minimum area of crown (1), important for radii around 0.5
            crown_area = min(crown_area,1963);                              // maximum area of crown (radius 25)

            //if(t_Crown_Depth > 2.0) cout << crown_area << " radius_layer: " << radius_layer << endl;
            
            /* now fill up all the sites until crown_area equals theoretically calculated crown_area (floor of it) */
            for(int i=0; i < crown_area; i++){
                int site_relative = LookUp_Crown_site[i];
                row = row_trunk + site_relative/51 - 25;
                col = col_trunk + site_relative%51 - 25;
                
                if(row >= 0 && row < rows && col >= 0 && col < cols){
                    site=col+cols*row+SBORD;
                    LAI3D[h][site] += t_dens * height_layer;
                }
            }
        }
    }
}
#else

/* the idea is to loop over a LookUp table for the respective crown voxels (in order of increasing distance) until area is filled up */

void Tree::CalcLAI() {
    if(t_age > 0){
        int crown_base = int(t_Tree_Height-t_Crown_Depth),
        crown_top = int(t_Tree_Height),
        row_trunk = t_site/cols,
        col_trunk = t_site%cols,
        col,row,site;
        int crown_area = int(PI * t_Crown_Radius * t_Crown_Radius);     // floor of crown_area to bound area accumulation
        
        crown_area = max(crown_area,1);                                 // minimum area of crown (1), important for radii around 0.5
        crown_area = min(crown_area,1963);                              // maximum area of crown (radius 25)
        
        /* now fill up all the sites until crown_area equals theoretically calculated crown_area (floor of it) */
        for(int i=0; i < crown_area; i++){
            int site_relative = LookUp_Crown_site[i];
            row = row_trunk + site_relative/51 - 25;
            col = col_trunk + site_relative%51 - 25;
            
            if(row >= 0 && row < rows && col >= 0 && col < cols){
                site=col+cols*row+SBORD;
                if(crown_top-crown_base == 0) {
                    LAI3D[crown_top][site] += t_dens*t_Crown_Depth;
                } else {
                    LAI3D[crown_top][site] += t_dens*(t_Tree_Height-crown_top);
                    LAI3D[crown_base][site] += t_dens*(crown_base+1-(t_Tree_Height-t_Crown_Depth));
                    if(crown_top-crown_base>=2){
                        for(int h=crown_base+1;h <= crown_top-1;h++){
                            LAI3D[h][site] += t_dens;    // loop over the crown depth
                        }
                    }
                }
            }
        }
    }
}
#endif
#else

void Tree::CalcLAI() {
    if(t_age>0) {
        int crown_base = int(t_Tree_Height-t_Crown_Depth),
        crown_top = int(t_Tree_Height),
        crown_r =int(t_Crown_Radius+0.5),
        row_trunk = t_site/cols,
        col_trunk = t_site%cols,
        xx,yy,col,row,site;
        
        // loop over the tree crown
        for(col=max(0,col_trunk-crown_r);col<=min(cols-1,col_trunk+crown_r);col++) {
            for(row=max(0,row_trunk-crown_r);row<=min(rows-1,row_trunk+crown_r);row++) {
                xx=col_trunk-col;
                yy=row_trunk-row;
                if(xx*xx+yy*yy<=crown_r*crown_r){  // check whether voxel is within crown
                    site=col+cols*row+SBORD;
                    if(crown_top-crown_base == 0) {
                        LAI3D[crown_top][site] += t_dens*t_Crown_Depth;
                    }
                    else{
                        LAI3D[crown_top][site] += t_dens*(t_Tree_Height-crown_top);
                        LAI3D[crown_base][site] += t_dens*(crown_base+1-(t_Tree_Height-t_Crown_Depth));

                        if(crown_top-crown_base>=2){
                            for(int h=crown_base+1;h <= crown_top-1;h++)
                                LAI3D[h][site] += t_dens;    // loop over the crown depth
                        }
                    }
                }
            }
        }
    }
}

#endif


/*###################################################
 ####  Computation of PPFD right above the tree  ####
 ####    called by Tree::Birth and Growth   in     ####
 ####################################################*/

/* v.2.3.: Tree::Fluxh() computes the average light flux received by a tree crown layer at height h , and also the average VPD and T it is thriving in (modified 1/02/2016)*/
/* v.2.4.0: added FLUX_AVG that depends both on density before and in current layer */

#ifdef CROWN_EXPANSION
void Tree::Fluxh(int h) {
    t_PPFD = 0.0;
    t_VPD  = 0.0;
    t_T    = 0.0;
    float absorb_prev=0.0;
#ifdef FLUX_AVG
    float absorb_delta=0.0;
#endif

    int row_trunk = t_site/cols,
    col_trunk = t_site%cols,
    col,row,site;
  
#ifdef CROWN_SHAPE
    int crown_base = int(t_Tree_Height - t_Crown_Depth);
    int crown_top = int(t_Tree_Height);
    int crown_center = int(t_Tree_Height - 0.5 * t_Crown_Depth);
    float crown_extent_base = crown_center - (t_Tree_Height - t_Crown_Depth);

    /* flux is always calculated above the layer, so h - 1 necessary to be in accordance with CalcLAI() */
    int h_crown = h - 1;
    
    /* calculating the height of the current layer */
    float height_layer = 1.0;
    if(crown_top == crown_base) height_layer = t_Crown_Depth;
    else if(h_crown == crown_top) height_layer = (t_Tree_Height-crown_top);
    else if(h_crown == crown_base) height_layer = (crown_base+1-(t_Tree_Height-t_Crown_Depth));

    /* calculating the radius of the current layer depending on the respective slopes */
    float radius_layer = t_Crown_Radius;
    if(crown_top == crown_base){}
    else if(h_crown == crown_base) radius_layer -= t_Crown_Slope_Bottom * crown_extent_base;
    else if(h_crown < crown_center ) radius_layer -= t_Crown_Slope_Bottom * (crown_center-h_crown);
    else radius_layer -= t_Crown_Slope_Top * (h_crown - crown_center);            // for h = crown_center full radius
    
    float crown_area_fl = maxf(PI * radius_layer * radius_layer,0.0);         // crown area for weighting
    int crown_area = int(crown_area_fl);                                         // floor of crown_area to bound area accumulation
    crown_area = max(crown_area,1);                                 // minimum area of crown (1), important for radii below ~0.5
    crown_area = min(crown_area,1963);                              // maximum area of crown (radius 25)
    
    t_Crown_Volume_layer = crown_area_fl * height_layer;     // gpp weighting
#else
    int crown_area = int(PI * t_Crown_Radius * t_Crown_Radius);     // floor of crown_area to bound area accumulation
    crown_area = max(crown_area,1);                                 // minimum area of crown (1), important for radii around 0.5
    crown_area = min(crown_area,1963);                              // maximum area of crown (radius 25)
#endif
    
    /* loop over LookUp table, until crown_area is reached */
    for(int i=0; i < crown_area; i++){
        int site_relative = LookUp_Crown_site[i];
        row = row_trunk + site_relative/51 - 25;
        col = col_trunk + site_relative%51 - 25;
        
        
        if(row >= 0 && row < rows && col >= 0 && col < cols){
            site=col+cols*row+SBORD;
            absorb_prev = LAI3D[h][site];
            absorb_prev = minf(absorb_prev,19.95);
#ifdef FLUX_AVG
            absorb_delta = LAI3D[h-1][site]-LAI3D[h][site];
            absorb_delta = minf(absorb_delta,9.95);
            int intabsorb = int(absorb_prev*20.0) + 400*int(absorb_delta*20.0);
#else
            int intabsorb = int(absorb_prev*20.0);
#endif
        
            t_PPFD += Wmax*LookUp_flux[intabsorb];
            t_VPD += VPDmax*LookUp_VPD[intabsorb];
            t_T += tmax - LookUp_T[intabsorb];
        }
    }
            
    float inv_crown_area = 1.0/crown_area;
    
    t_PPFD *= inv_crown_area;
    t_VPD  *= inv_crown_area;
    t_T    *= inv_crown_area;
}

#else

void Tree::Fluxh(int h) {
    t_PPFD = 0.0;
    t_VPD  = 0.0;
    t_T    = 0.0;
    float absorb_prev=0.0;
#ifdef FLUX_AVG
    float absorb_delta=0.0;
#endif

    int count=0,
    xx,yy,radius_int;

    radius_int=int(t_Crown_Radius);

    if(radius_int == 0) {
        count=1;
        if (h < HEIGHT) absorb_prev = LAI3D[h][t_site+SBORD];
        absorb_prev = minf(absorb_prev,19.5);
        // absorb = 0.0 by default
        int intabsorb=int(absorb_prev*20.0);
        t_PPFD = Wmax*LookUp_flux[intabsorb];
        t_VPD  = VPDmax*LookUp_VPD[intabsorb];
        t_T    = tmax - LookUp_T[intabsorb];
    }
    else {
        int row0,col0;
        row0=t_site/cols;
        col0=t_site%cols;
        for(int col=max(0,col0-radius_int);col<min(cols,col0+radius_int+1);col++) {
            for(int row=max(0,row0-radius_int);row<min(rows,row0+radius_int+1);row++) {
                //loop over the tree crown
                xx=col0-col;
                yy=row0-row;
                if(xx*xx+yy*yy <= radius_int*radius_int) {
                    //is the voxel within crown?
                    count++;
                    int site=col+cols*row+SBORD;
                    absorb_prev = LAI3D[h][site];
                    absorb_prev = minf(absorb_prev,19.95);
#ifdef FLUX_AVG
                    absorb_delta = LAI3D[h-1][site]-LAI3D[h][site];
                    absorb_delta = minf(absorb_delta,9.95);
                    int intabsorb = int(absorb_prev*20.0) + 400*int(absorb_delta*20.0);
#else
                    int intabsorb = int(absorb_prev*20.0);
#endif
                    
                    t_PPFD += Wmax*LookUp_flux[intabsorb];
                    t_VPD  += VPDmax*LookUp_VPD[intabsorb];
                    t_T    += tmax - LookUp_T[intabsorb];
                }
            }
        }
    }
    t_PPFD *= 1.0/float(count);
    t_VPD  *= 1.0/float(count);
    t_T    *= 1.0/float(count);
        
}

#endif

#ifdef INTRASPECIFIC

/*############################################
 ####            GPP calculation          ####
 ####    v.2.4.0: now at tree level       ####
 #############################################*/

float Tree::DeathRateNDD(float PPFD, float dbh, int nppneg, float ndd) {
    float dr=0;
    float basal=m*(1-t_wsg);
    float dd=deltaD*ndd*(1-2*dbh/t_dmax);
    
    dr=basal;
    if (nppneg > t_leaflifespan) {
        dr+=1.0/timestep;
    }
    if (dd > 0) {
        dr+=dd;
    }
    
    return dr*timestep;
}

    
float Tree::DeathRate(float PPFD, float dbh, int nppneg) {
    float dr=0.0;
    float basal=maxf(m-m1*t_wsg,0.0);

    dr=basal;
    
    if (nppneg > t_leaflifespan) dr+=1.0/timestep;
    
    if (iter == int(nbiter-1))
        output[26]<< t_wsg << "\t" << basal << "\t"  << dbh << "\t"  << dr   <<  "\n";
    
    return dr*timestep;
}

float Tree::GPPleaf(float PPFD, float VPD, float T) {
    /* v.2.3.0: theta defined as a global variable */
    //theta=0.7;   // this is the fixed value of theta used by von Caemmerer 2000
    
    //float theta=0.76+0.018*T-0.00037*T*T;         // theta, but temperature dependent cf. Bernacchi et al 2003 PCE
    
    /* Parameters for Farquhar model, with temperature dependencies */
    int convT= int(iTaccuracy*T); // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---
    
    //if(convT>500 || isnan(convT) || convT <0) cout << convT << " | T: " << T << " | PPFD: " << PPFD << " | VPD: " << VPD << endl;
    float KmT    =LookUp_KmT[convT];
    float GammaT =LookUp_GammaT[convT];
    tempRday    +=t_Rdark*LookUp_tempRday[convT];
    float VcmaxT =t_Vcmax*LookUp_VcmaxT[convT];
    float JmaxT  =t_Jmax*LookUp_JmaxT[convT];
    
    t_fci=g1/(g1+sqrt(VPD));
    
    /* FvCB model */
    
    float I=alpha*PPFD;
    float J = (I+JmaxT-sqrt((JmaxT+I)*(JmaxT+I)-4.0*theta*JmaxT*I))*0.5/theta;
    float A = minf(VcmaxT/(t_fci+KmT),0.25*J/(t_fci+2.0*GammaT))*(t_fci-GammaT);
    
    return A;
}

/* dailyGPPleaf returns the assimilation rate (computed from Tree::GPPleaf) averaged across the daily fluctuations in climatic conditions (light, VPD and T), in micromoles C/m^2/s */

float Tree::dailyGPPleaf(float PPFD, float VPD, float T) {
    float ppfde,dailyA=0.0;
    
    for(int i=0; i<24; i++) {
        ppfde=PPFD*daily_light[i];
        if(ppfde > 0.1)
            // new v.2.3.0: compute GPP only if enough light is available threshold is arbitrary, but set to be low: in full sunlight ppfd is aroung 700 W/m2, and even at dawn, it is ca 3% of the max value, or 20 W/m2. The minimum threshold is set to 0.1 W/m2
            // Future update: compute slightly more efficiently, using 3-hourly values? This will have to be aligned with climate forcing layers (e.g. NCAR)
            dailyA+=Tree::GPPleaf(ppfde, VPD*daily_vpd[i], T*daily_T[i]);
        //the 6 lines in comment below corresponds to a finer version in which the multiplier is computed and used every 48 half hour, ie. with the corresponding environment instead of assuming a constant multiplier correponding the one at maximum incoming irradiance
        //float hhA=0;
        //hhA=GPPleaf(PPFD*daily_light[i], VPD*daily_vpd[i], T*daily_T[i]);
        //float alpha=phi*PPFD*daily_light[i]/hhA;
        //float D=klight*dens*CD;
        //hhA*=alpha/(D*(alpha-1))*log(alpha/(1+(alpha-1)*exp(-D)));
        //dailyA+=hhA;
    }
    //daily_light is the averaged (across one year, meteo station Nouragues DZ) and normalized (from 0 to 1) daily fluctuation of light, with half-hour time step, during the day time (from 7am to 7pm, ie 12 hours in total), same for daily_vpd and daily_T. Taking into account these daily variation is necessary considering the non-linearity of FvCB model
    
    dailyA*=0.0417;                                 // 0.0417=1/24 (24=12*2 = number of half hours in the 12 hours of daily light)
    tempRday*=0.0417;
    return dailyA;
}

/* faster, whole crown GPP calculation */
float Tree::dailyGPPcrown(float PPFD, float VPD, float T, float dens, float CD) {
    float ppfde,dailyA=0.0;
    
    for(int i=0; i<24; i++) {
        ppfde=PPFD*daily_light[i];
        if(ppfde > 0.1)
            // new v.2.3.0: compute GPP only if enough light is available threshold is arbitrary, but set to be low: in full sunlight ppfd is aroung 700 W/m2, and even at dawn, it is ca 3% of the max value, or 20 W/m2. The minimum threshold is set to 0.1 W/m2
            // Future update: compute slightly more efficiently, using 3-hourly values? This will have to be aligned with climate forcing layers (e.g. NCAR)
            dailyA+=Tree::GPPleaf(ppfde, VPD*daily_vpd[i], T*daily_T[i]);
        //the 6 lines in comment below corresponds to a finer version in which the multiplier is computed and used every 48 half hour, ie. with the corresponding environment instead of assuming a constant multiplier correponding the one at maximum incoming irradiance
        //float hhA=0;
        //hhA=GPPleaf(PPFD*daily_light[i], VPD*daily_vpd[i], T*daily_T[i]);
        //float alpha=phi*PPFD*daily_light[i]/hhA;
        //float D=klight*dens*CD;
        //hhA*=alpha/(D*(alpha-1))*log(alpha/(1+(alpha-1)*exp(-D)));
        //dailyA+=hhA;
    }
    //daily_light is the averaged (across one year, meteo station Nouragues DZ) and normalized (from 0 to 1) daily fluctuation of light, with half-hour time step, during the day time (from 7am to 7pm, ie 12 hours in total), same for daily_vpd and daily_T. Taking into account these daily variation is necessary considering the non-linearity of FvCB model
    
    float alpha=phi*PPFD/GPPleaf(PPFD, VPD, T);             //alpha is a non-dimensional figure used to compute the multiplier below
    float D=klight*dens*CD;                                 //D is a non-dimensional figure used to compute the multiplier below
    dailyA*=alpha/(D*(alpha-1))*log(alpha/(1+(alpha-1)*exp(-D)));  // the FvCB assimilation rate computed at the top of the tree crown is multiplied by a multiplier<1, to account for the lower rate at lower light level within the crown depth. This multiplier is computed assuming that change in photosynthetic assimilation rate within a tree crown is mainly due to light decrease due to self-shading following a Michealis-menten relationship (ie. we assume that 1/ the change is not due to changes in VPD or temperature, which are supposed homogeneous at the intra-crown scale, and 2/ that other tree contributions to light decrease is neglected).
    
    dailyA*=0.0417;                                 // 0.0417=1/24 (24=12*2 = number of half hours in the 12 hours of daily light)
    tempRday*=0.0417;
    return dailyA;
}
#endif


/*############################################
 ####            Tree growth              ####
 ####         called by UpdateTree        ####
 #############################################*/

void Tree::Growth() {
    
    /* Flux Tree variables */
    t_GPP=0.0;
    t_NPP=0.0;
    t_Rday=0.0;
    t_Rnight=0.0;
    t_Rstem=0.0;
    tempRday=0.0;
    
    /* Local environmental Tree variables */
    t_PPFD=0.0;
    t_VPD=0.0;
    t_T=0.0;
    
    /* variables for flux computations */
    
    t_age+= timestep;                               /* new v.2.2: increments are not 1 yr, but the duration of the timestep (usually 1 or <1, i.e. 1/12 if monthly, 1/365 if daily */
    
    /* computation of average t_GPP (per area) from the sum of GPP of each tree crown layer: */
    /* Farquhar is applied once per tree crown (at the top of the crown) (48 times per timestepifdef DAILIGHT, once per timestep if not), a pultiplier is net used to account for the decrease in photosynthetic rate with light decreae within the tree crown. */
    /* v.2.3.1 -- fast GPP calculation option. In addition, the daily computation of the Farquhar model is now the default option (_DAILYLIGHT is deprecated) */
    
    if(_FASTGPP){
        Fluxh(int(t_Tree_Height)+1);
#ifdef INTRASPECIFIC
        t_GPP = Tree::dailyGPPcrown(t_PPFD, t_VPD, t_T, t_dens, t_Crown_Depth);
#else
        t_GPP = (t_s->dailyGPPcrown(t_PPFD, t_VPD, t_T, t_dens, t_Crown_Depth));
#endif
        t_Rday += tempRday;
        tempRday=0.0;
    }
    else CalcRespGPP();  /* new v.2.4: t_GPP and t_Rday are updated in separate function (thus adaptable for different modules) */
        
    /* Computation of GPP. New v.2.2: assumes an efficiency of 0.5 for young and old leaves vs. 1 for mature leaves */
    
    float effLA=0.5*(t_leafarea+t_matureLA)*189.3*timestep;     // effLA is the scaling factor used for all fluxes new in v.2.3.0
    t_GPP*=effLA;
    
    /* new v.2.2. sapwood thickness (useful to compute stem respiration) */
    float sapthick=0.04;
    if (t_dbh < 0.08) sapthick=0.5*t_dbh;
    
    /* new v.2.2 stem respiration -- update lookup v230 */
    int convT= int(iTaccuracy*temp); // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---
    int convTnight= int(iTaccuracy*tnight); // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---
    t_Rstem=sapthick*(t_dbh-sapthick)*(t_Tree_Height-t_Crown_Depth)*LookUp_Rstem[convT];
    t_Rday *= effLA*0.40;
    
#ifdef INTRASPECIFIC
    t_Rnight=t_Rdark*effLA*LookUp_Rnight[convTnight];
#else
    t_Rnight=(t_s->s_Rdark)*effLA*LookUp_Rnight[convTnight];
#endif
    t_NPP = 0.75*(t_GPP - 1.5*(t_Rday+t_Rnight+t_Rstem));
    /* Rleaf=Rday+Rnight is multiplied by 1.5 to also account for fine root respiration (cf as in Fyllas et al 2014 and Malhi 2012); Rstem is multiplied by 1.5 to account for coarse root respiration (according to the shoot root biomass ratio of 0.2 - Jérôme's paper in prep- and also to branch respiration (Meir & Grace 2002, Cavaleri 2006, Asao 2005). */

    if(t_NPP<0.0){
        t_NPPneg++;
        t_NPP=0.0;  /* in case of eventlog activation, t_NPP should be set back to 0.0 after the UpdateLeafDynamics() function to allow for recording of negative values */
        /* v.2.3.0 -- Line of code below was odd. If NPP <0.0, then to ensure C balance it should be simply reset to NPP=0 at this stage */
        //t_NPP=t_GPP - 1.5*(t_Rday+t_Rnight+t_Rstem); REMOVED AS OF v.2.3.0.a4
    }
    else {
        t_NPPneg=0;
        /**** NPP allocation to wood and tree size increment *****/
        UpdateTreeBiometry();
    }
    
    /**** NPP allocation to leaves *****/
    UpdateLeafDynamics();

    /* Output for control purposes */
    /* v.2.4.0: legacy from version 2.2, may not work for all parameterisations, hence deactivated, could be cleaned up in next version */
    
//    if(!_OUTPUT_reduced){
//        if (iter == 2) OutputTreeStandard(output[28]);
//        if (iter == int(nbiter/2)) OutputTreeStandard(output[29]);
//        if (iter == int(nbiter-1)) OutputTreeStandard(output[30]);
//        
//        if (t_site==2500) OutputTreeStandard(output[11]);
//        if (t_site==10380) OutputTreeStandard(output[12]);
//        if (t_site==100950) OutputTreeStandard(output[13]);
//        if (t_site==12090) OutputTreeStandard(output[14]);
//        if (t_site==120090) OutputTreeStandard(output[15]);
//        if (t_site==150667) OutputTreeStandard(output[16]);
//    }
    
}

    /*####################################################
     ####      GPP and Respiration Calculation        ####
     ####         called by Tree::Growth              ####
     #####################################################*/

/* new in v.2.4.0: extra functions */

#ifdef CROWN_SHAPE
/* the main change for crown_shape is the weighting of GPP by crown layers, not only due to variation in vertical, but also in lateral extent */
void Tree::CalcRespGPP(){
    int crown_above_base=int(t_Tree_Height-t_Crown_Depth)+1; // for flux above crown base
    int crown_above_top=int(t_Tree_Height)+1;                // for flux above crown top
    
    for(int h=crown_above_base; h<=crown_above_top; h++) {
        Fluxh(h);
#ifdef INTRASPECIFIC
        t_GPP+=t_Crown_Volume_layer*Tree::dailyGPPleaf(t_PPFD, t_VPD, t_T);
#else
        t_GPP+=t_Crown_Volume_layer*t_s->dailyGPPleaf(t_PPFD, t_VPD, t_T);
#endif
        t_Rday+=t_Crown_Volume_layer*tempRday;
        tempRday=0.0;
    }
    float i_Crown_Volume=1.0/float(t_Crown_Volume);
    t_GPP   *=i_Crown_Volume;
    t_Rday  *=i_Crown_Volume;
}

#elif defined(GPP_PROPORTIONAL)
/* introduces weighting according to how far tree reaches into the respective flux layer */
/* ToDo: top and baselayer contribution could be integrated into loop */
void Tree::CalcRespGPP(){
    int crown_base = int(t_Tree_Height-t_Crown_Depth);
    int crown_top = int(t_Tree_Height);
    
    if(crown_base == crown_top){
        int crown_above_top = crown_top+1;
        Fluxh(crown_above_top);
#ifdef INTRASPECIFIC
        t_GPP=Tree::dailyGPPleaf(t_PPFD, t_VPD, t_T);
#else
        t_GPP=t_s->dailyGPPleaf(t_PPFD, t_VPD, t_T);
#endif
        t_Rday=tempRday;
        tempRday=0.0;
        //if(t_site == 8771) cout << t_site << " | Height! " << t_Tree_Height << " ddbh: " << t_ddbh << " t_GPP: " << t_GPP << " Depth: " << t_Crown_Depth << " t_GPP: " << t_GPP << " t_Rday: " << t_Rday << " t_dens: "<< t_dens << endl;
    } else {
        int crown_above_base = crown_base+1;
        int crown_above_top = crown_top+1;
        
        float toplayer_proportion = 0.0, baselayer_proportion = 0.0;
        toplayer_proportion = t_Tree_Height - float(crown_top);
        baselayer_proportion = float(crown_above_base) + t_Crown_Depth - t_Tree_Height;
        
        Fluxh(crown_above_base);
#ifdef INTRASPECIFIC
        t_GPP+=baselayer_proportion*Tree::dailyGPPleaf(t_PPFD, t_VPD, t_T);
#else
        t_GPP+=baselayer_proportion*t_s->dailyGPPleaf(t_PPFD, t_VPD, t_T;
#endif
        t_Rday+=baselayer_proportion*tempRday;
        tempRday=0.0;

        Fluxh(crown_above_top);
#ifdef INTRASPECIFIC
        t_GPP+=toplayer_proportion*Tree::dailyGPPleaf(t_PPFD, t_VPD, t_T);
#else
        t_GPP+=toplayer_proportion*t_s->dailyGPPleaf(t_PPFD, t_VPD, t_T);
#endif
        t_Rday+=toplayer_proportion*tempRday;
        tempRday=0.0;

        for(int h=crown_above_base+1; h<crown_above_top; h++) {
          //if(t_site == 8771) cout << t_site << " | Height! " << t_Tree_Height << " ddbh: " << t_ddbh << " t_GPP: " << t_GPP << " Depth: " << t_Crown_Depth << " base: " << 1.0 << " top: " << 1.0 << " t_GPP: " << t_GPP << " t_Rday: " << " t_dens: "<< t_dens << endl;
          Fluxh(h);
#ifdef INTRASPECIFIC
          t_GPP+=Tree::dailyGPPleaf(t_PPFD, t_VPD, t_T);
#else
          t_GPP+=t_s->dailyGPPleaf(t_PPFD, t_VPD, t_T);
#endif
          t_Rday+=tempRday;
          tempRday=0.0;
        }
        float inv_Crown_Depth=1.0/t_Crown_Depth;
        t_GPP   *=inv_Crown_Depth;
        t_Rday  *=inv_Crown_Depth;

    }
}

#else

void Tree::CalcRespGPP(){
    int crown_above_base=int(t_Tree_Height-t_Crown_Depth)+1; // for flux above crown base
    int crown_above_top=int(t_Tree_Height)+1;                // for flux above crown top
    for(int h=crown_above_base; h<=crown_above_top; h++) {
        Fluxh(h);
#ifdef INTRASPECIFIC
        t_GPP+=Tree::dailyGPPleaf(t_PPFD, t_VPD, t_T);
#else
        t_GPP+=t_s->dailyGPPleaf(t_PPFD, t_VPD, t_T);
#endif
        t_Rday+=tempRday;
        tempRday=0.0;
    }
    float inb_layer=1.0/float(crown_above_top-crown_above_base+1);
    t_GPP   *=inb_layer;
    t_Rday  *=inb_layer;
}
    
#endif
                                                      
    
/*####################################################
 ####       Leaf dynamics and C allocation        ####
 ####         called by Tree::Growth              ####
 #####################################################*/

void Tree::UpdateLeafDynamics() {
    
    // make a standalone function for leaf dynamics & litter?
    
    /**** NPP allocation to leaves *****/                                       /* rk: in this current scheme of leaf demography and phenology in three leaf age classes: only the old leaves generate treefall, and the dynamic of leaves cycle is generated by the dynamic of NPP, with a total leaf biomass varying - as opposed to De Weirdt et al 2012 in ORCHIDEE, but as in Wu et al 2016 but importantly without prescribing litterfall- */
    
    /* The following line is to convert the NPP allocated to leaves (falloccanopy is the fraction of biomass assumed to be alloacted to canopy (leaves+reproductive organs+twigs) at each timestep - Malhi et al 2011-, 68% of which is allocated to leaves - chave et al 2008, Chave et al 2010-), in new leaf area (2 is to convert carbon mass in biomass and LMA to convert leaf biomass into leaf area).*/
#ifdef INTRASPECIFIC
    float flush=2.0*maxf(t_NPP,0.0)*falloccanopy*0.68/t_LMA;
    
    /* litter module */
    float flush_fine = flush/float(leafdem_resolution);
    float lambda_young = 1.0/(t_time_young * leafdem_resolution) ;
    float lambda_mature = 1.0/(t_time_mature * leafdem_resolution) ;
    float lambda_old = 1.0/(t_time_old * leafdem_resolution) ;
#else
    float flush=2.0*maxf(t_NPP,0.0)*falloccanopy*0.68/(t_s->s_LMA);
    
    /* litter module */
    float flush_fine = flush/float(leafdem_resolution);
    float lambda_young = 1.0/(t_s->s_time_young * leafdem_resolution) ;
    float lambda_mature = 1.0/(t_s->s_time_mature * leafdem_resolution) ;
    float lambda_old = 1.0/(t_s->s_time_old * leafdem_resolution) ;
#endif
    
    t_litter = 0.0;
    
    for (int i = 0; i < leafdem_resolution; i++){
        float new_litter = lambda_old * t_oldLA ;
        float new_young = flush_fine;
        float new_mature = t_youngLA * lambda_young ;
        float new_old    = t_matureLA * lambda_mature ;

        t_youngLA += new_young - new_mature ;
        t_matureLA += new_mature - new_old ;
        t_oldLA += new_old - new_litter ;
        t_litter += new_litter;
    }
    t_leafarea = t_youngLA + t_matureLA + t_oldLA ;
    
    /* update t_dens */
    
#ifdef INTRASPECIFIC
    t_litter*=t_LMA;
#else
    t_litter*=t_s->s_LMA;
#endif
    
#ifdef CROWN_SHAPE
    t_dens=t_leafarea/t_Crown_Volume;
#else
    float crownvolume=PI*t_Crown_Radius*LH*t_Crown_Radius*LH*t_Crown_Depth*LV;
    t_dens=t_leafarea/crownvolume;
#endif
}

void Tree::UpdateTreeBiometry(){
    /* New standalone function in v.2.3.0 */
    
    /* Tree dbh increment */
    t_ddbh=0.0;

#ifdef INTRASPECIFIC
    /* volume in m^3: the first factor of 2 is to convert C into biomass. the 1/s_ wsg to convert biomass into volume. the 1e-6 term converts cm^3 into m^3 (the sole metric unit in the model). fallocwood is the fraction of biomass allocated to aboveground wood (stem + branches) growth. For the time being, we shall assume that a fixed proportion of NPP is allocated into AGB production. Currently, 0.20=%biomasse allocated to stem increment could be a global variable, even though this % allocation could in fact vary with resouce variation/co-limitation*/
    
    /* taking into account wood elements recycling (ex. fallen branches etc...) */
    /* t_ddbh = flor( volume* 4.0/( 3.0*PI*t_dbh*LH*t_Tree_Height*LV ) )* NH; */
    
    float volume=2.0*t_NPP/(t_wsg) * fallocwood * 1.0e-6;
    if (t_dbh>t_dbh_thresh) volume*=flor(3.0-2.0*t_dbh/t_dbh_thresh);
    t_ddbh = flor(volume/(0.559*t_dbh*LH*t_Tree_Height*LV*(3.0-t_dbh/(t_dbh+t_ah))) )* NH;
    
    /* With V=pi*r^2*h, increment of volume = dV = 2*pi*r*h*dr + pi*r^2*dh */
    /* With isometric growth assumption (ddbh/dbh=dh/h)and dbh=2*r: dV=3/4*pi*dbh*h*ddbh, ddbh in m, it follows: ddbh = 4/3 * V = 4/3 * 1/(pi*dbh*h)   */
    
    t_dbh += t_ddbh;
    
    /* update of tree height */
    /* alternative calculation in concordance with isometric growth assumption: dh = h*ddbh/dbh */
    /* t_Tree_Height += t_Tree_Height*t_ddbh/t_dbh; */
    
    t_Tree_Height = minf(t_intraspecific_multiplier_height * t_hmax * t_dbh/(t_dbh + t_ah),HEIGHT-1);
#else
    float volume=2.0*t_NPP/(t_s->s_wsg) * fallocwood * 1.0e-6;
    if (t_dbh>t_dbh_thresh) volume*=flor(3.0-2.0*t_dbh/t_dbh_thresh);
    t_ddbh = flor( volume/(0.559*t_dbh*LH*t_Tree_Height*LV*(3.0-t_dbh/(t_dbh+t_ah))) )* NH;
    t_dbh += t_ddbh;
    t_Tree_Height = t_hmax * t_dbh/(t_dbh + t_ah);
#endif
  
    /* update of tree crown depth -- allometry deduced from Piste Saint-Elie dataset */

    if(t_Tree_Height < 5.0) {t_Crown_Depth = de0 + 0.17 * (t_Tree_Height-H0);}
    else {t_Crown_Depth = de0+0.26*(t_Tree_Height-H0) - 0.09 * (5.0-H0);} /* allometry deduced from Piste Saint-Elie dataset (unpublished) */
    
    t_Crown_Radius = exp(1.9472 + 0.5925*log(t_dbh)); /* this is crown allometry derived from data set provided by Jucker et al. 2016 (Global Change Biology) */
    // t_Crown_Radius  = 0.80+10.47*t_dbh-3.33*t_dbh*t_dbh; /* allometry deduced from Piste Saint-Elie dataset */

#ifdef INTRASPECIFIC
    t_Crown_Radius *= t_intraspecific_multiplier_CR;
    t_Crown_Depth *= t_intraspecific_multiplier_CD;
#endif
 
#ifdef CROWN_SHAPE
    /* update the crown slopes */
    int crown_base = int(t_Tree_Height - t_Crown_Depth);
    int crown_top = int(t_Tree_Height);
    int crown_center = int(t_Tree_Height-0.5*t_Crown_Depth);
    float crown_extent_top = t_Tree_Height - crown_center;
    float crown_extent_base = crown_center - (t_Tree_Height - t_Crown_Depth);
    
    if(t_Crown_Depth > 2.0){        // only consider cases, where more than two layers are reached by the crown
        /* slopes are calculated as horizontal/vertical slopes (inverse of more intuitive definition vertical/horizontal), reasoning: when shape_crown = 1.0, this slope is still defined */
        /* the slope at the base part of the crown is steeper, meaning that the trapezoid narrows down faster (limit: triangle/cone) */
        t_Crown_Slope_Top = t_Crown_Radius * (1.0 - shape_crown) / crown_extent_top;
        t_Crown_Slope_Bottom = min(2 * t_Crown_Slope_Top, t_Crown_Radius/crown_extent_base);
    } else {
        t_Crown_Slope_Top = 0.0;
        t_Crown_Slope_Bottom = 0.0;
    }
    
    t_Crown_Volume = 0.0;

    for(int h=crown_base;h<crown_top+1;h++){
        /* calculating the height of the current layer */
        float height_layer = 1.0;
        if(crown_top == crown_base) height_layer = t_Crown_Depth;
        else if(h == crown_top) height_layer = (t_Tree_Height-crown_top);
        else if(h == crown_base) height_layer = (crown_base+1-(t_Tree_Height-t_Crown_Depth));

        /* calculating the radius of the current layer depending on the respective slopes */
        float radius_layer = t_Crown_Radius;
        if(crown_top == crown_base){}
        else if(h == crown_base) radius_layer -= t_Crown_Slope_Bottom * crown_extent_base;
        else if(h < crown_center ) radius_layer -= t_Crown_Slope_Bottom * (crown_center-h);
        else radius_layer -= t_Crown_Slope_Top * (h - crown_center);            // for h = crown_center full radius

        float crown_area_fl = maxf(PI * radius_layer * radius_layer,0.0);
        
        t_Crown_Volume += crown_area_fl* height_layer;
    }
#endif
}


/*####################################################
 ####           Death of the tree                ####
 ####         called by Tree::Update             ####
 ####################################################*/

void Tree::Death() {
    
    /* This function records basic properties of a tree that has died and sempties the variables at the tree site */
    /* It does not, however, apply the destructor to the tree object */
    
    /* new v.2.4: statistics are now calculated inside the Death() function */
    /* tree death statistics */
    nbdead_n1++;
    if(t_dbh*LH>0.1) nbdead_n10++;
    if(t_dbh*LH>0.3) nbdead_n30++;
    
    /* New v.2.2. new outputs */
    if(!_OUTPUT_reduced) {
        if(iter == 2) output[23] << "N\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
        if(iter == int(nbiter/2)) output[24]<< "N\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
        if(iter == int(nbiter-1)) output[25]<< "N\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
    }

    t_sp_lab = 0;
    t_age = 0;
    t_hurt = 0;
    t_NPP=t_GPP=t_Rday=t_Rnight=t_Rstem=t_PPFD=t_VPD=t_T=0.0; /* new v.2.3 */
    t_dbh = t_Tree_Height = t_Crown_Radius = t_Crown_Depth= 0.0;
    
    if(_BASICTREEFALL) t_Ct = 0.0;
    
    if ((t_s->s_nbind)>0) (t_s->s_nbind)--;
    nblivetrees--;
    t_s = NULL;
    
}


/*#################################
 ####      Seed dispersal      ####
 ####  called by UpdateField   ####
 #################################*/

void Tree::DisperseSeed(){
    /* New v.2.0 reproduction can only occur for trees that receive enough
     light (twice the LCP) */
    /* New v.2.1 threshold of maturity is defined as a size threshold
     (and not age as before), following Wright et al 2005 JTE */
    if((t_dbh>=t_dbhmature)&&(t_PPFD>2.0*(t_s->s_LCP))) {
        float rho,theta_angle;
        int nbs=0;
        if(_SEEDTRADEOFF){
            nbs=int(t_NPP*2.0*falloccanopy*0.08*0.5*(t_s->s_iseedmass)*0.05);    /* nbs is the number of seeds produced at this time step; it is computed from the NPP (in g) allocated to reproductive organs -fruits and seeds-, *2 is to convert in biomass,  * 0.40 is to obtain the NPP allocated to canopy (often measured as litterfall), drawn from Malhi et al 2011 Phil. trans roy. Soc. and 0.08 is the part of litterfall corresponding the fruits+seeds, drawn from Chave et al 2008 JTE; assumed to be twice the biomass dedicated to seeds only (ie without fruits), and then divided by the mass of a seed to obtain the number of seeds */
            //nbs=(int)nbs;
        }
        else nbs=nbs0;
        //else nbs=int(t_NPP*2*falloccanopy*0.08*0.5); /* test 17/01/2017: use a factor to translate NPP into seeds produced, but not species specific, not linked to mass of grains */
        
        for(int ii=0;ii<nbs;ii++){                                                 /* Loop over number of produced seeds */
            
            rho = 2.0*((t_s->s_ds)+t_Crown_Radius)*float(sqrt(fabs(log(genrand2()*iPi))));    /* Dispersal distance rho: P(rho) = rho*exp(-rho^2) */
            theta_angle = float(twoPi*genrand2());                                                /* Dispersal angle theta */
            t_s->FillSeed(flor(int(rho*cos(theta_angle))+t_site%cols), /* column */               /* Update of field s_Seed */
                          flor(int(rho*sin(theta_angle))+t_site/cols));      /* line */
            
        }
    }
    
}
//#endif


/*##################################
 ####   Tree death and growth   ####
 ####   called by UpdateTree    ####
 ##################################*/

void Tree::Update() {

    int death;

    if(t_age) {
        if(t_dbh > 0.1) nbtrees_n10++;
        
        /* v.2.4.0: outputs have been moved to Death() function */
        
#ifdef INTRASPECIFIC
        if(_NDD)
            death = int(genrand2()+DeathRateNDD(t_PPFD, t_dbh, t_NPPneg, t_NDDfield[t_sp_lab]));
        else
            death = int(genrand2()+DeathRate(t_PPFD, t_dbh, t_NPPneg));
#else
        if(_NDD)
            death = int(genrand2()+t_s->DeathRateNDD(t_PPFD, t_dbh, t_NPPneg, t_NDDfield[t_sp_lab]));
        else
            death = int(genrand2()+t_s->DeathRate(t_PPFD, t_dbh, t_NPPneg));
#endif
        
        if(death){
            Death();
        } else {
            Growth();   // v.2.4: t_hurt is now updated in the TriggerTreefallSecondary() function
        }
    }
}
    
    
/*##################################
 ####           Treefall         ####
 #### called by UpdateTreefall   ####
 ####################################*/
    

/* NEW in TROLL v.2.4: FallTree() function has become Treefall() function, calculation of angle and treefall outside of function, and damages are now added up from several treefalls */
    
void Tree::Treefall(float angle) {
    
    /* treefall statistics */
    nbTreefall1++;
    if(t_dbh*LH>0.1) nbTreefall10++;
    if(t_dbh*LH>0.3) nbTreefall30++;
    
    int xx,yy;
    int row0,col0,h_int, r_int;
    float h_true = t_Tree_Height*LV;
    h_int = int(h_true*NH);
    row0=t_site/cols;
    col0=t_site%cols;
    
    /* update of Thurt field at the site of the tree, for consistency */
    Thurt[0][t_site+sites] = int(t_Tree_Height);
    
    /* fallen stem destructs other trees */
    for(int h=1;h<h_int;h++) {                      // loop on the fallen stem (horizontally)
        xx=int(flor(col0+h*cos(angle)));          // get projection in col (= xx) direction, where xx is absolute location
        if(xx<cols){
            yy=   int(row0+h*sin(angle));         // get projection in row (= yy) direction, where yy is absolute location
            Thurt[0][xx+(yy+rows)*cols] += int(t_Tree_Height);
            // Thurt[0] where the stem fell, calculation: xx+(yy+rows)*cols= xx + yy*cols + rows*cols = xx + yy*cols + sites / NEW in v.2.4: addition of damage instead of setting equal in order to account for cumulative damage (several treefalls hitting the same site)
        }
    }
    
    /* fallen crown destructs other trees, less damaging than stem */
    xx=col0+int((h_true*NH-t_Crown_Radius)*cos(angle));
    yy=row0+int((h_true*NH-t_Crown_Radius)*sin(angle));
    r_int = int(t_Crown_Radius);
    for(int col=max(0,xx-r_int);col<min(cols,xx+r_int+1);col++) { // loop on the fallen crown (horizontally)
        for(int row=yy-r_int;row<yy+r_int+1;row++) {
            if((col-xx)*(col-xx)+(row-yy)*(row-yy)<r_int*r_int) Thurt[0][col+(row+rows)*cols] += int((t_Tree_Height-t_Crown_Radius*NV*LH)*0.5); // less severe damage than stem / NEW in v.2.4: addition of damage instead of setting equal in order to account for cumulative damage (several treefalls hitting the same site)
        }
    }
    /* v.2.4.0: outputs have been moved to Death() function */
    Death();
}
    
void Tree::Couple(float &c_forceflex, float &angle) {
    int site2,quadist,h0,xx,yy, radius_int,h_int;
    float fx, fy,temp,lai;
    radius_int = int(t_Crown_Radius);
    h_int = int(t_Tree_Height);
    h0 = int(t_Tree_Height-t_Crown_Depth);
    if(radius_int){
        int row0,col0;
        row0=t_site/cols;
        col0=t_site%cols;
        fx = fy = 0.0;
        for(int col=max(0,col0-radius_int);col<min(cols,col0+radius_int+1);col++) {
            for(int row=row0-radius_int;row<=row0+radius_int;row++) {
                xx=col0-col;
                yy=row0-row;
                quadist = xx*xx+yy*yy;
                if((quadist<=radius_int*radius_int)&&quadist) {
                    //site2 = col+cols*(row+RMAX); //modif 23/03/2011
                    site2 = col+cols*row+SBORD;
                    for(int h=h0;h<=h_int;h++) {
                        if(h_int<HEIGHT) lai = LAI3D[h_int][site2]-LAI3D[h_int+1][site2];
                        else  lai = LAI3D[h_int][site2];
                        if(lai>dens) { // needs to be changed when TREEFALL is revised
                            temp = 1.0/sqrt(float(quadist));
                            if(temp>0.0) {
                                fx += float(xx)*temp;
                                fy += float(yy)*temp;
                            }
                        }
                    }
                }
            }
        }
        c_forceflex = int(sqrt(fx*fx+fy*fy)*t_Tree_Height);
        if(fx!=0.0) angle=atan2(fy,fx);
        else angle = Pis2*sgn(fy);
    }
    else{c_forceflex = 0; angle = 0.0; }
}
    
/*##################################
 ####   NINO   Tree felling     ####
 #### called by SelectiveLogging ####
 ####################################*/


void Tree::FellTree() {

	float t_angle = float(twoPi*genrand2()); // random angle
	float h_true = t_Tree_Height*LV;
    int xx, yy, row0, col0, r_int, h_int=int(h_true*NH);    

    if(t_dbh*LH>0.1) nbTreefall10++;
    Thurt[0][t_site+sites] = int(t_Tree_Height); // Thurt[0] saves the integer tree height, here exactly at the place where the tree fell...
    Tlogging[2][t_site] = 1;
    row0=t_site/cols;       /* fallen stem destructs other trees */
    col0=t_site%cols;
    for(int h=1;h<h_int;h++) { // loop on the fallen stem (horizontally)
        xx=int(flor(col0+h*cos(t_angle))); // get projection in col (= xx) direction, where xx is absolute location
        if(xx<cols){
            yy=int(row0+h*sin(t_angle)); // get projection in row (= yy) direction, where yy is absolute location
            Thurt[0][xx+(yy+rows)*cols] = int(t_Tree_Height); // Thurt[0] where the stem fell, calculation: xx+(yy+rows)*cols= xx + yy*cols + rows*cols = xx + yy*cols + sites
            if(xx<cols && yy<rows) Tlogging[2][xx+yy*cols] = 1; // Tfell[0] is used to represent gaps for gaps damages modelling
        }
    }
    xx=col0+int((h_true*NH-t_Crown_Radius)*cos(t_angle)); // where crown ends/starts fallen crown destructs other trees */
    yy=row0+int((h_true*NH-t_Crown_Radius)*sin(t_angle));
    r_int = int(t_Crown_Radius);
    for(int col=max(0,xx-r_int);col<min(cols,xx+r_int+1);col++) { /* loop on the fallen crown (horizontally) */
        for(int row=yy-r_int;row<yy+r_int+1;row++) {
            if((col-xx)*(col-xx)+(row-yy)*(row-yy)<r_int*r_int) {
            	Thurt[0][col+(row+rows)*cols] = int((t_Tree_Height-t_Crown_Radius*NV*LH)*0.5);
            	if(xx<cols && yy<rows) Tlogging[2][xx+yy*cols] = 1;
            }
        }
    }
    Death();
}    
  /* DOES IT HAVE TO BE UPDATED  */
/*#####################################################
 ####      For Average and OutputField             ####
 ######################################################*/

void Tree::Average() {
    if(t_age>0) {
        if(t_dbh*LH >= 0.1) {
            (t_s->s_output_field[1])++;
            t_s->s_output_field[6] += t_dbh*LH*t_dbh*LH;
        }
        if(t_dbh*LH >= 0.3) (t_s->s_output_field[2])++;
        t_s->s_output_field[3] += t_dbh*LH*t_dbh*LH;
        t_s->s_output_field[4] += t_NPP*1.0e-6;
        t_s->s_output_field[5] += t_GPP*1.0e-6;
#ifdef INTRASPECIFIC
        t_s->s_output_field[7] += 0.0673*pow(t_wsg*t_Tree_Height*LV*t_dbh*t_dbh*LH*LH*10000, 0.976);  // this is the allometric equ 4 in Chave et al. 2014 Global Change Biology to compute above ground biomass
#else
        t_s->s_output_field[7] += 0.0673*pow(t_s->s_wsg*t_Tree_Height*LV*t_dbh*t_dbh*LH*LH*10000, 0.976);  // this is the allometric equ 4 in Chave et al. 2014 Global Change Biology to compute above ground biomass
#endif
        t_s->s_output_field[8] += t_Rday*1.0e-6;
        t_s->s_output_field[9] += t_Rnight*1.0e-6;
        t_s->s_output_field[10] += t_Rstem*1.0e-6;
        t_s->s_output_field[11] += t_litter*1.0e-6;
    }
}

void Tree::histdbh() {
    if(t_age) nbdbh[int(100.*t_dbh*LH)]++;
    // compute the diameter distribution density
    // where dbh is in cm (it is in number of horizontal cells throughout the code)
    // values are always rounded down (so nbdbh[30] gives you trees with more than 30 cm dbh, and less than 31))
}


/*#####################################################
 ####      User Output for Control Purposes        ####
 ######################################################*/

void Tree::OutputTreeStandard(fstream& output){
    output << iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight << "\t" << t_site << "\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)+1][t_site+SBORD];
}

void Tree::OutputTreeStandard(){
    cout << iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight << "\t" << t_site << "\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)+1][t_site+SBORD] << endl;
}

/* Class objects */

Species *S=NULL;
Tree *T=NULL;


/*############################################
 ############################################
 ############     MAIN PROGRAM    ###########
 ############################################
 ############################################*/

int main(int argc,char *argv[]) {
    
    if(_DAILYLIGHT == 1) cout << "Activated Module: Dailylight" << endl;
    if(_FASTGPP == 1) cout << "Activated Module: FastGPP" << endl;
    if(_BASICTREEFALL == 1) cout << "Activated Module: BASICTREEFALL" << endl;
    if(_TREEFALL == 1) cout << "Activated Module: TREEFALL" << endl;
    if(_NDD == 1) cout << "Activated Module: NDD" << endl;
    if(_SEEDTRADEOFF == 1) cout << "Activated Module: SEEDTRADEOFF" << endl;
    if(_FromData == 1) cout << "Activated Module: FromData" << endl;
    if(_INPUT_fullFinal == 1) cout << "Activated Module: INPUT_fullFinal" << endl;
    if(_OUTPUT_reduced == 1) cout << "Activated Module: OUTPUT_reduced" << endl;
    if(_OUTPUT_last100 == 1) cout << "Activated Module: OUTPUT_last100" << endl;
    if(_OUTPUT_fullFinal == 1) cout << "Activated Module: OUTPUT_fullFinal" << endl;
    if(_DISTURBANCE == 1) cout << "Activated Module: Disturbance" << endl;
    if(_LOGGING == 1) cout << "Activated Module: Logging" << endl;
    
    /***********************/
    /*** Initializations ***/
    /***********************/
#ifdef MPI   /* Lookup processor number / total number of processors */
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&p_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
#else
    mpi_rank = 0;
    mpi_size = 1;
#endif
    easympi_rank = 0;
    
    for(int argn=1;argn<argc;argn++){ /* Arguments of the input and output files */
        if(*argv[argn] == '-'){
            switch(*(argv[argn]+1)){
                case 'i':
                    bufi = argv[argn]+2;
                    break;
                case 'm':                       /* new v.2.4; initialisation of climate parameters from separate file */
                    bufi_climate = argv[argn]+2;
                    break;
                case 'o':
                    buf = argv[argn]+2;
                    break;
                case 'f':                      /* new v.2.3: initialisation from field, 'f' for "forest", "field data" */
                    bufi_data = argv[argn]+2;
                    break;
                 case 's':                      /*  NINO new sylviculture module: initialisation from field, 's' for "sylviculture" parameters */
                    bufi_sylviculture = argv[argn]+2;
                    break;  
                case 'n':
                    easympi_rank=atoi(argv[argn]+2); /* new v.2.2 */
            }
        }
    }
    cout << "Easy MPI rank: " << easympi_rank << endl;

    int t = (int) time(NULL);
    int seed = 3*t+2*(easympi_rank+1)+1;

    if(_NONRANDOM == 1) seed = 1;   /* new v.2.4 */
    
    sgenrand2i(seed);
    sgenrand2(seed);
    
#if defined (INTRASPECIFIC) || defined(DCELL)
    // Stuff for constant number generator
    const gsl_rng_type *Trandgsl;
    gsl_rng_env_setup();
    Trandgsl = gsl_rng_default;
    gslrand = gsl_rng_alloc (Trandgsl);
    
    unsigned long int t2 = (unsigned long int) time(NULL);
    unsigned long int seed2 = 3*t2 + 2*(easympi_rank+1)+1;
    
    if(_NONRANDOM == 1) seed2 = 1;
    
    gsl_rng_set(gslrand, seed2);
    
#endif
    
    cout<< "On proc #" << easympi_rank << " seed: " << seed <<  " rng: "<< genrand2() << endl;
    cout << genrand2()<< endl;
    cout << genrand2()<< endl;
    
    // input files
    
    sprintf(inputfile,"%s",bufi);
    sprintf(inputfile_climate,"%s",bufi_climate);
    if(_FromData){
        sprintf(inputfile_data,"%s",bufi_data);
    }
        if(_DISTURBANCE || _LOGGING){
        sprintf(inputfile_sylviculture,"%s",bufi_sylviculture);
    } /* NINO*/ 
    
    sprintf(outputinfo,"%s_%i_par.txt",buf, easympi_rank);                     /* Files with general output info */
    cout<< "On proc #" << easympi_rank << " seed: " << seed <<  " rng: "<< genrand2() <<  endl;
    out.open(outputinfo, ios::out);
    if(!out) cerr<< "ERROR with par file"<< endl;
    sprintf(outputinfo,"%s_%i_info.txt",buf, easympi_rank);
    out2.open(outputinfo, ios::out);
    if(!out2) cerr<< "ERROR with info file"<< endl;
    sprintf(outputinfo,"%s_%i_paramsylviculture.txt",buf, easympi_rank);		/* sylviculture module parameters NINO */
    out3.open(outputinfo, ios::out);
    if(!out3) cerr<< "ERROR with sylviculture par file"<< endl;
    Initialise();           /* Read global parameters */
    
    
#ifdef INTRASPECIFIC
    float max_intraspecific_height=0.0, min_intraspecific_height=1000.0,
    max_intraspecific_CR=0.0, min_intraspecific_CR = 1000.0,
    max_intraspecific_CD=0.0, min_intraspecific_CD = 1000.0,
    max_intraspecific_P=0.0, min_intraspecific_P = 1000.0,
    max_intraspecific_N=0.0, min_intraspecific_N = 1000.0,
    max_intraspecific_LMA=0.0, min_intraspecific_LMA = 1000.0,
    max_intraspecific_wsg=0.0, min_intraspecific_wsg = 1000.0,
    max_intraspecific_dmax=0.0, min_intraspecific_dmax = 1000.0;
    float variation_height=0.0, variation_CR=0.0, variation_CD=0.0, variation_P=0.0, variation_N=0.0, variation_LMA=0.0,variation_wsg=0.0,variation_dmax=0.0;
    for(int i=0;i<100000;i++){
#ifdef INTRASPECIFIC_covariance
        double var_height_dbl, var_CR_dbl;
        gsl_ran_bivariate_gaussian(gslrand, sigma_height, sigma_CR, corr_CR_height, &var_height_dbl, &var_CR_dbl);
        gsl_ran_multivariate_gaussian(gslrand, mu_N_P_LMA, mcov_N_P_LMA, variation_N_P_LMA);
        variation_height = float(var_height_dbl);
        variation_CR = float(var_CR_dbl);
        variation_N = gsl_vector_get(variation_N_P_LMA, 0);
        variation_P = gsl_vector_get(variation_N_P_LMA, 1);
        variation_LMA = gsl_vector_get(variation_N_P_LMA, 2);
        variation_CD = float(gsl_ran_gaussian(gslrand, sigma_CD));
        variation_wsg = float(gsl_ran_gaussian(gslrand, sigma_wsg));
        variation_dmax = float(gsl_ran_gaussian(gslrand, sigma_dmax));
#else
        variation_height = float(gsl_ran_gaussian(gslrand, sigma_height));
        variation_CR = float(gsl_ran_gaussian(gslrand, sigma_CR));
        variation_N = float(gsl_ran_gaussian(gslrand, sigma_N));
        variation_P = float(gsl_ran_gaussian(gslrand, sigma_P));
        variation_LMA = float(gsl_ran_gaussian(gslrand, sigma_LMA));
        variation_CD = float(gsl_ran_gaussian(gslrand, sigma_CD));
        variation_wsg = float(gsl_ran_gaussian(gslrand, sigma_wsg));
        variation_dmax = float(gsl_ran_gaussian(gslrand, sigma_dmax));
#endif
        d_intraspecific_height[i] = exp(variation_height);
        d_intraspecific_CR[i] = exp(variation_CR);
        d_intraspecific_N[i] = exp(variation_N);
        d_intraspecific_P[i] = exp(variation_P);
        d_intraspecific_LMA[i] = exp(variation_LMA);
        d_intraspecific_CD[i] = exp(variation_CD);
        d_intraspecific_wsg[i] = variation_wsg;             // normal, not log-normal
        d_intraspecific_dmax[i] = exp(variation_dmax);
        //        d_intraspecific_height[i] = exp(float(gsl_ran_gaussian(gslrand, sigma_height)));
        max_intraspecific_height = maxf(max_intraspecific_height,d_intraspecific_height[i]);
        min_intraspecific_height = minf(min_intraspecific_height,d_intraspecific_height[i]);
        max_intraspecific_CR = maxf(max_intraspecific_CR,d_intraspecific_CR[i]);
        min_intraspecific_CR = minf(min_intraspecific_CR,d_intraspecific_CR[i]);
        max_intraspecific_N = maxf(max_intraspecific_N,d_intraspecific_N[i]);
        min_intraspecific_N = minf(min_intraspecific_N,d_intraspecific_N[i]);
        max_intraspecific_P = maxf(max_intraspecific_P,d_intraspecific_P[i]);
        min_intraspecific_P = minf(min_intraspecific_P,d_intraspecific_P[i]);
        max_intraspecific_LMA = maxf(max_intraspecific_LMA,d_intraspecific_LMA[i]);
        min_intraspecific_LMA = minf(min_intraspecific_LMA,d_intraspecific_LMA[i]);
        max_intraspecific_CD = maxf(max_intraspecific_CD,d_intraspecific_CD[i]);
        min_intraspecific_CD = minf(min_intraspecific_CD,d_intraspecific_CD[i]);
        max_intraspecific_wsg = maxf(max_intraspecific_wsg,d_intraspecific_wsg[i]);
        min_intraspecific_wsg = minf(min_intraspecific_wsg,d_intraspecific_wsg[i]);
        max_intraspecific_dmax = maxf(max_intraspecific_dmax,d_intraspecific_dmax[i]);
        min_intraspecific_dmax = minf(min_intraspecific_dmax,d_intraspecific_dmax[i]);
    }
    cout << "Max and min allometry deviation, lognormal (height): " << max_intraspecific_height << " | " << min_intraspecific_height << endl;
    cout << "Max and min allometry deviation, lognormal (crown radius): " << max_intraspecific_CR << " | " << min_intraspecific_CR << endl;
    cout << "Max and min trait deviation, lognormal (N): " << max_intraspecific_N << " | " << min_intraspecific_N << endl;
    cout << "Max and min trait deviation, lognormal (P): " << max_intraspecific_P << " | " << min_intraspecific_P << endl;
    cout << "Max and min trait deviation, lognormal (LMA): " << max_intraspecific_LMA << " | " << min_intraspecific_LMA << endl;
    cout << "Max and min allometry deviation, normal (crown depth): " << max_intraspecific_CD << " | " << min_intraspecific_CD << endl;
    cout << "Max and min trait deviation, normal (wsg): " << max_intraspecific_wsg << " | " << min_intraspecific_wsg << endl;
    cout << "Max and min trait deviation, lognormal (dmax): " << max_intraspecific_dmax << " | " << min_intraspecific_dmax << endl;
#endif
    
    AllocMem();             /* Memory allocation */
    
    if(_FromData){
        InitialiseFromData();   /* Initial configuration of the forest, read from data */
    }
    
    BirthInit();            /* Initial configuration of the forest */
    out.close();
    
    cout << "klight is: " << klight << endl;
    cout << "CO2 concentration is: " << Cair << endl << endl;
    
    /***********************/
    /*** Evolution loop  ***/
    /***********************/
    
    double start_time,stop_time, duration=0.0;           /* for simulation duration */
    stop_time = clock();
    for(iter=0;iter<nbiter;iter++) {
        start_time = stop_time;
        
        Evolution();
        stop_time = clock();
        duration +=flor(stop_time-start_time);
        
        /****** Output for 100 year analysis ******/
        
        if(_OUTPUT_last100 && nbiter>100 && iter > (nbiter-101)) OutputSnapshotDetail(output[33]);  // 100 years development
        
        /****** Final Outputs ******/
        
        if(iter == nbiter-2){
            OutputSnapshot(output[10]);                         // Final Pattern
            if(_OUTPUT_reduced) OutputSnapshot10cm(output[18]);
            if(!_OUTPUT_reduced)  OutputSpeciesParameters(output[18]);
            if(_OUTPUT_fullFinal)	OutputSnapshotFullFinal(output[37]);				// Full Final Pattern (TROLL save) NINO : check if output 37 available
        }
        
    }
    
    /*************************/
    /*** End of simulation ***/
    /*************************/
    
    float durf = duration/double(CLOCKS_PER_SEC);        /* output of the effective CPU time */
#ifdef MPI
    MPI_Reduce(&durf,&durf,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
#endif
    if(!mpi_rank) {
        cout << "\n";
#ifdef MPI
        out2 << "Number of processors : "<< mpi_size << "\n";
#endif
        out2 << "Average computation time : "<< durf/float(mpi_size) << " seconds.\n";
        out2 << "End of simulation.\n";
        cout << "\nNumber of processors : "<< mpi_size << "\n";
        cout << "Average computation time : "<< durf/float(mpi_size) << " seconds.\n";
        cout << "End of simulation.\n";
    }
    out2.close();
    
    /*Free dynamic memory */ /* added in oct2013 */
    FreeMem();
    
#ifdef easyMPI
    MPI::Finalize();
#endif
    
    exit(0);
}


/*###########################################
 ############################################
 #######  Initialization routines    ########
 ############################################
 ############################################*/

void Initialise() {
    
    /*** Initialization of the simulation parameters ***/
    /***************************************************/
    
    fstream In(inputfile, ios::in);
    
    if (In) {
        for(int line=0;line<4;line++) In.getline(buffer,128,'\n');
        
        /*General parameters */
        In >> cols; In.getline(buffer,128,'\n');
        In >> rows; In.getline(buffer,128,'\n');
        sites = rows*cols;
#ifdef DCELL
        In >> length_dcell; In.getline(buffer,128,'\n');
        sites_per_dcell = length_dcell*length_dcell;
        nbdcells = int(sites/sites_per_dcell);
        linear_nb_dcells = int(cols/length_dcell);
        cout << "Number of dcells" << nbdcells << endl;
        cout << "Lin number of dcells" << linear_nb_dcells << endl;
#endif
        In >> nbiter; In.getline(buffer,128,'\n');
        In >> iterperyear; In.getline(buffer,128,'\n');
        timestep=1.0/float(iterperyear);
        cout << "iterperyear " << iterperyear << endl;
        In >> NV; In.getline(buffer,128,'\n');
        In >> NH; In.getline(buffer,128,'\n');
        LV = 1.0/NV;
        LH = 1.0/NH;
        In >> nbout; In.getline(buffer,128,'\n');
        if(nbout) freqout = nbiter/nbout;
        In >> numesp; In.getline(buffer,128,'\n');
        In >> p_nonvert; In.getline(buffer,128,'\n');
        In.getline(buffer,128,'\n');
        In.getline(buffer,128,'\n');
        
        /*Characters shared by species */
        In >> klight; In.getline(buffer,128,'\n');
        theta = 0.70; /* v.2.3.0 This could be added in the input file, just like klight */
        In >> phi; In.getline(buffer,128,'\n');
        In >> g1; In.getline(buffer,128,'\n');
        In >> vC; In.getline(buffer,128,'\n');
        In >> DBH0; In.getline(buffer,128,'\n');
        In >> H0; In.getline(buffer,128,'\n');
        In >> ra0; In.getline(buffer,128,'\n');
        In >> ra1; In.getline(buffer,128,'\n');
        In >> de0; In.getline(buffer,128,'\n');
        In >> de1; In.getline(buffer,128,'\n');
        In >> dens; In.getline(buffer,128,'\n');
        In >> fallocwood; In.getline(buffer,128,'\n');
        In >> falloccanopy; In.getline(buffer,128,'\n');
        In >> Cseedrain; In.getline(buffer,128,'\n');
        In >> nbs0; In.getline(buffer,128,'\n');
        
#ifdef INTRASPECIFIC
        In >> sigma_height; In.getline(buffer,128,'\n');
        In >> sigma_CR; In.getline(buffer,128,'\n');
        In >> sigma_CD; In.getline(buffer,128,'\n');
        In >> sigma_P; In.getline(buffer,128,'\n');
        In >> sigma_N; In.getline(buffer,128,'\n');
        In >> sigma_LMA; In.getline(buffer,128,'\n');
        In >> sigma_wsg; In.getline(buffer,128,'\n');
        In >> sigma_dmax; In.getline(buffer,128,'\n');
#ifdef INTRASPECIFIC_covariance
        In >> corr_CR_height; In.getline(buffer,128,'\n');
        In >> corr_N_P; In.getline(buffer,128,'\n');
        In >> corr_N_LMA; In.getline(buffer,128,'\n');
        In >> corr_P_LMA; In.getline(buffer,128,'\n');
        
        // convert correlations to covariances
        cov_N_P = corr_N_P * sigma_N * sigma_P;
        cov_N_LMA = corr_N_LMA * sigma_N * sigma_LMA;
        cov_P_LMA = corr_P_LMA * sigma_P * sigma_LMA;
        
        // Initialise covariance matrix for N, P, LMA
        
        mcov_N_P_LMA = gsl_matrix_alloc(3,3);
        gsl_matrix_set(mcov_N_P_LMA,0,0,sigma_N*sigma_N);
        gsl_matrix_set(mcov_N_P_LMA,0,1,cov_N_P);
        gsl_matrix_set(mcov_N_P_LMA,0,2,cov_N_LMA);
        gsl_matrix_set(mcov_N_P_LMA,1,0,cov_N_P);
        gsl_matrix_set(mcov_N_P_LMA,1,1,sigma_P*sigma_P);
        gsl_matrix_set(mcov_N_P_LMA,1,2,cov_P_LMA);
        gsl_matrix_set(mcov_N_P_LMA,2,0,cov_N_LMA);
        gsl_matrix_set(mcov_N_P_LMA,2,1,cov_P_LMA);
        gsl_matrix_set(mcov_N_P_LMA,2,2,sigma_LMA*sigma_LMA);
        
        cout << "\nCovariance matrix N,P,LMA: " << endl;
        for(int mrow=0; mrow<3;mrow++){
            cout << gsl_matrix_get(mcov_N_P_LMA, mrow, 0) << "\t" << gsl_matrix_get(mcov_N_P_LMA, mrow, 1) << "\t" << gsl_matrix_get(mcov_N_P_LMA, mrow, 2) << endl;
        }
        
        // Cholesky decomposition for multivariate draw
        gsl_linalg_cholesky_decomp1(mcov_N_P_LMA);

        cout << "\nCovariance matrix N,P,LMA (after Cholesky decomposition) " << endl;
        for(int mrow=0; mrow<3;mrow++){
            cout << gsl_matrix_get(mcov_N_P_LMA, mrow, 0) << "\t" << gsl_matrix_get(mcov_N_P_LMA, mrow, 1) << "\t" << gsl_matrix_get(mcov_N_P_LMA, mrow, 2) << endl;
        }
        
        // allocating mean and result vectors for multivariate draw
        mu_N_P_LMA = gsl_vector_alloc(3);
        for(int j=0;j<3;j++){
            gsl_vector_set(mu_N_P_LMA,j,0.0);
        }
        variation_N_P_LMA = gsl_vector_alloc(3);
        
#endif
#endif
        
        /* new in v.2.4.0 */
        In >> leafdem_resolution; In.getline(buffer,128,'\n');
        In >> p_tfsecondary; In.getline(buffer,128,'\n');
        In >> hurt_decay; In.getline(buffer,128,'\n');
        
        In >> m; In.getline(buffer,128,'\n');
        In >> m1; In.getline(buffer,128,'\n');
        In >> Cair; In.getline(buffer,128,'\n');
        iCair = 1.0/Cair;
        if (_NDD) {
            In >> R; In.getline(buffer,128,'\n');
            In >> deltaR; In.getline(buffer,128,'\n');
            In >> deltaD; In.getline(buffer,128,'\n');
        }
        
        DBH0 *= NH;
        H0   *= NV;
        ra0  *= NH;
        de0  *= NV;
        alpha = 4.0*phi;
        /* apparent quantum yield to electron transport in mol e-/mol photons */
        // see Mercado et al 2009 , the conversion of the apparent quantum yield in micromolCO2/micromol quantum into micromol e-/micxromol quantum is done by multipliyng by 4, since four electrons are needed to regenerate RuBP.
        /* alpha is fixed at 0.3 mol e-/mol photons in Medlyn et al 2002
         but see equ8 and Appendix 1 in Farquahr et al 1980: it seems that alpha should vary with leaf thickness: there is a fraction of incident light which is lost by absorption by other leaf parts than the chloroplast lamellae, and this fraction f may increase with leaf thickness.
         With the values of the paper: alpha= 0.5*(1-f)=0.5*(1-0.23)=0.385, but this is a theoretical value and observations often report lower values (see ex discussion in medlyn et al 2005 Tree phsyiology, Lore Veeryckt values, Mercado et al 2009 Table 10, Domingues et al. 2014)*/
#ifdef CROWN_SHAPE
        shape_crown = 1.0;  // crown shape parameter (1.0 corresponds to cylindric shape, 0 to cones in top and bottom part)
#endif
        
    }
    
    else {
        cout<< "ERROR with the input file of parameters" << endl;
    }
    
    /*** Information in file info ***/
    /***********************************/
    
    if(!mpi_rank){
        out2 << "\nTROLL simulator\n\n";
        out2 << "\n   2D discrete network: horizontal step = " << LH
        << " m, one tree per "<< LH*LH << " m^2 \n\n";
        out2 << "\n   Tree : (t_dbh,t_Tree_Height,t_Crown_Radius,t_Crown_Depth) \n\n";
        out2 << "\n            + one species label \n\n";
        out2 << " Number of sites      : "<<rows<<"x"<<cols<<"\n";
        out2 << " Number of iterations : "<<nbiter<<"\n";
        out2 << " Duration of timestep : "<<timestep<<" years\n";
        out2 << " Number of Species    : "<<numesp << "\n\n";
        out2.flush();
    }
    
    
    /*** Initialization of trees ***/
    /*******************************/
    
    if(NULL==(T=new Tree[sites])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Tree" << endl;
    }
    
    if(_NDD){
        for(int site=0;site<sites;site++) {
            if (NULL==(T[site].t_NDDfield = new float[numesp+1])) cerr<<"!!! Mem_Alloc\n";
            for(int ii=0;ii<(numesp+1);ii++) T[site].t_NDDfield[ii]=0;
        }
    }
    
    
    /*** Initialization of species ***/
    /*********************************/
    
    int sp;
    if(NULL==(S=new Species[numesp+1])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Species" << endl;
    }
    
    for(int line=0;line<3;line++) In.getline(buffer,128,'\n');                           /* Read species parameters (ifstream In) */
    for(sp=1;sp<=numesp;sp++) {
        S[sp].Init(sp,In);
    }
    
    /* Close ifstream In */
    In.close();
    
    /*** Initialization of environmental variables ***/
    /*************************************************/
    
    fstream InClim(inputfile_climate, ios::in);
    
    if (InClim) {
        for(int line=0;line<4;line++) InClim.getline(buffer,128,'\n');
        
        /* normalized daily variation in PPFD, VPD, T */
        for (int i=0; i<=23; i++) InClim >> daily_light[i];
        InClim.getline(buffer,128,'\n');
        for (int i=0; i<=23; i++) InClim >> daily_vpd[i];
        InClim.getline(buffer,128,'\n');
        for (int i=0; i<=23; i++) InClim >> daily_T[i];

        InClim.getline(buffer,128,'\n');
        InClim.getline(buffer,128,'\n');
        InClim.getline(buffer,128,'\n');
        
        /* monthly averages */
    
        if(NULL==(Temperature=new float[iterperyear])) {                                // rk, the current structure of code suppose that environment is periodic (a period = a year), if one want to make climate vary, with interannual variation and climate change along the simulation, one just need to provide the full climate input of the whole simulation (ie number of columns=iter and not iterperyear) and change iterperyear by nbiter here.
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc Temperature" << endl;
        }
        
        for (int i=0; i<iterperyear; i++) InClim >> Temperature[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "Temperature" << i << "\t"  << Temperature[i] <<  "\t";
        //cout<<endl;
        
        if(NULL==(DailyMaxTemperature=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc DailyMaxTemperature" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> DailyMaxTemperature[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "DailyMaxTemperature" << i << "\t"  << DailyMaxTemperature[i] << "\t";
        //cout<<endl;
        
        if(NULL==(NightTemperature=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc NightTemperature" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> NightTemperature[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "NightTemperature" << i << "\t"  << NightTemperature[i] << "\t";
        //cout<<endl;
        
        if(NULL==(Rainfall=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc Rainfall" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> Rainfall[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "Rainfall" << i << "\t"  << Rainfall[i] << "\t";
        //cout<<endl;
        
        if(NULL==(WindSpeed=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc WindSpeed" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> WindSpeed[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "WindSpeed" << i << "\t"  << WindSpeed[i] << "\t";
        //cout<<endl;
        
        if(NULL==(MaxIrradiance=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc Irradiance" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> MaxIrradiance[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "MaxIrradiance" << i << "\t"  << MaxIrradiance[i] << "\t";
        //cout<<endl;
        
        if(NULL==(MeanIrradiance=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc Irradiance" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> MeanIrradiance[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "MeanIrradiance" << i << "\t"  << MeanIrradiance[i] << "\t";
        //cout<<endl;
        
        if(NULL==(SaturatedVapourPressure=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc SaturatedVapourPressure" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> SaturatedVapourPressure[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "SaturatedVapourPressure" << i << "\t"  << SaturatedVapourPressure[i] << "\t";
        //cout<<endl;
        
        if(NULL==(VapourPressure=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc VapourPressure" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> VapourPressure[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "VapourPressure" << i << "\t"  << VapourPressure[i] << "\t";
        //cout<<endl;
        
        if(NULL==(VapourPressureDeficit=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc VapourPressureDeficit" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> VapourPressureDeficit[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "VapourPressureDeficit" << i << "\t"  << VapourPressureDeficit[i] << "\t";
        //cout<<endl;
        
        if(NULL==(DailyVapourPressureDeficit=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc DailyVapourPressureDeficit" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> DailyVapourPressureDeficit[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "DailyVapourPressureDeficit" << i << "\t"  << DailyVapourPressureDeficit[i] << "\t";
        //cout<<endl;
        
        if(NULL==(DailyMaxVapourPressureDeficit=new float[iterperyear])) {
            cerr<<"!!! Mem_Alloc\n";
            cout<<"!!! Mem_Alloc DailyMaxVapourPressureDeficit" << endl;
        }
        for (int i=0; i<iterperyear; i++) InClim >> DailyMaxVapourPressureDeficit[i];
        InClim.getline(buffer,128,'\n');
        
        //for (int i=0; i<iterperyear; i++) cout<< "DailyMaxVapourPressureDeficit" << i << "\t"  << DailyMaxVapourPressureDeficit[i] << "\t";
        //cout<<endl;
        
        temp=Temperature[iter%iterperyear];
        tmax=DailyMaxTemperature[iter%iterperyear];
        tnight=NightTemperature[iter%iterperyear];
        precip=Rainfall[iter%iterperyear];
        WS=WindSpeed[iter%iterperyear];
        Wmax=MaxIrradiance[iter%iterperyear]*1.678;       // 1.678 is to convert irradiance from W/m2 to micromol of PAR /s /m2, the unit used in the FvCB model of photosynthesis
        Wmean=MeanIrradiance[iter%iterperyear];            // still in W/m2
        e_s=SaturatedVapourPressure[iter%iterperyear];
        e_a=VapourPressure[iter%iterperyear];
        VPDbasic=VapourPressureDeficit[iter%iterperyear];
        VPDday=DailyVapourPressureDeficit[iter%iterperyear];
        VPDmax=DailyMaxVapourPressureDeficit[iter%iterperyear];
    } else {
        cout<< "ERROR with the input file of climate parameters" << endl;
    }

    /* Close ifstream InClim */
    InClim.close();
    
    nbTbins=500;
    float Taccuracy=0.1;
    iTaccuracy=1.0/Taccuracy;
    cout << "Built-in maximal temperature: " << float(nbTbins)*Taccuracy <<endl;
    if(NULL==(LookUp_KmT=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_KmT" << endl;
    if(NULL==(LookUp_GammaT=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_GammaT" << endl;
    if(NULL==(LookUp_tempRday=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_tempRday" << endl;
    if(NULL==(LookUp_VcmaxT=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_VcmaxT" << endl;
    if(NULL==(LookUp_JmaxT=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_JmaxT" << endl;
    if(NULL==(LookUp_Rstem=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_Rstem" << endl;
    if(NULL==(LookUp_Rnight=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_Rnight" << endl;
    for(int i=0;i<nbTbins;i++) { // loop over "T" in GPPleaf()
        float temper=float(i)*Taccuracy;
        LookUp_KmT[i] = 404.0*exp(((temper-25.0)/(298*0.00831*(273+temper)))*59.36)*
        (1+210*1.0/248.0*exp(-(temper-25.0)/(298*0.00831*(273+temper))*35.94))*iCair;
        LookUp_GammaT[i]=37.0*exp(((temper-25.0)/(298*0.00831*(273+temper)))*23.4)*iCair;
        LookUp_tempRday[i]=exp((temper-25.0)*0.1*log(3.09-0.0215*(25.0+temper)));
        LookUp_VcmaxT[i]=exp(26.35-65.33/(0.00831*(temper+273.15)));
        LookUp_JmaxT[i]=exp(17.57-43.54/(0.00831*(temper+273.15)));
        LookUp_Rstem[i]=39.6*378.7*PI*timestep*exp((temper-25.0)/10.0*log(2.0));
        LookUp_Rnight[i]=exp((temper-25.0)*0.1*log(3.09-0.0215*(25.0+temper)));
        /* exp((temp-25)/10*log(2)) is the temperature dependency of Rstem, supposing a constant Q10=2, according to Ryan et al 1994 and Meir & Grace 2002
         exp((tnight-25)*0.1*log(3.09-0.0215*(25+tnight))) is the temperature dependencies used by Atkin 2015 (equ1) */
        
        
    }

    
#ifdef FLUX_AVG
    /* look up table for flux averaging/integration */
    /* division into absorption prior to current voxel (absorb_prev) and absorption in current voxel (absorb_delta) */
    /* prior absorption has a maximum of 20 m2/m3, while absorption within one voxel can maximally reach 10 m2/m3*/
    
    if(NULL==(LookUp_flux=new float[80000])) cerr<<"!!! Mem_Alloc LookUp_flux" << endl;
    if(NULL==(LookUp_VPD=new float[80000])) cerr<<"!!! Mem_Alloc LookUp_VPD" << endl;
    if(NULL==(LookUp_T=new float[80000])) cerr<<"!!! Mem_Alloc LookUp_VPD" << endl;
    for(int i=0;i<400;i++) { // loop over "absorb" in Fluxh()
        for(int j=0;j<200;j++){
            float absorb_prev=float(i)/20.0;
            float absorb_delta=float(j)/20.0;
            if(absorb_delta==0){
                // if the voxel does not contain any plant matter, values are constant across the voxel, e.g. just the top value calculated from absorb_prev
                LookUp_flux[i+400*j] = exp(-klight*absorb_prev);
                LookUp_VPD[i+400*j] = 0.25 + sqrt(maxf(0.0 , 0.08035714*(7.0-absorb_prev)));
                LookUp_T[i+400*j] = 0.4285714 * (minf(7.0,absorb_prev));
            } else {
                // if the voxel does contain plant matter, then average values will be computed
                // for voxels of 1 unit length depth, this corresponds just to the integral over LAI, which can be decomposed into a constant absorb_prev and a linearly increasing absorb_delta
                // once LAI reaches the critical value of 7.0, VPD And T do not decrease anymore, hence the distinction between two cases
                LookUp_flux[i+400*j] = exp(-klight * absorb_prev) * (1.0 - exp(-klight * absorb_delta)) / (klight * absorb_delta);
                if(absorb_prev+absorb_delta >= 7){
                    // this could be calculated more exactly, but difference is negligible
                    LookUp_VPD[i+400*j] = 0.25;
                    LookUp_T[i+400*j] = 3.0;  // 0.4285714 * 7.0
                } else {
                    LookUp_VPD[i+400*j] = 0.25 + (0.188982/absorb_delta) * (pow((7.0 - absorb_prev),1.5) - pow((7.0 - absorb_prev - absorb_delta),1.5));
                    LookUp_T[i+400*j] = 0.4285714 * (absorb_prev + 0.5 * absorb_delta);
                }
            }
        }
    }

#else
    if(NULL==(LookUp_flux=new float[400])) cerr<<"!!! Mem_Alloc LookUp_flux" << endl;
    if(NULL==(LookUp_VPD=new float[400])) cerr<<"!!! Mem_Alloc LookUp_VPD" << endl;
    if(NULL==(LookUp_T=new float[400])) cerr<<"!!! Mem_Alloc LookUp_T" << endl;
    for(int i=0;i<400;i++) { // loop over "absorb" in Fluxh()
        float absorbance=float(i)/20.0;
        LookUp_flux[i]=exp(-klight*absorbance);
        LookUp_VPD[i]= 0.25+sqrt(maxf(0.0 , 0.08035714*(7.0-absorbance)));
        // this expressions results from fit of observations of relationships between VPD and height within dense canopy (HOBO data on COPAS tower, Shuttleworth et al 1985; Camargo & Kapos 1995 journal of Tropical Ecology)
        LookUp_T[i] =  0.4285714*(minf(7.0,absorbance));
        // 0.4285714=3/7, assuming deltaT between the top canopy and dense understorey is constant = 3°C, could be refined.
        
    }
#endif
    
#ifdef CROWN_EXPANSION
    /* new in v.2.4: LookUp table that gives the sites within a crown in order of distance from the center */
    /* crowns can thus be assembled from inside out, in a radial fashion */
    
    int Crown_dist[2601];                                           // this saves the distances from the center of the crown
    int extent = 25;                                                // maximum extent of crowns (after test: either enlarge or allocate dynamically)
    int extent_full = 2*25 + 1;
    int index_crown = 0, xx, yy, site_rel, dist;
    Crown_dist[index_crown] = 0;                                    // this is the distance of the center of the crown from the center (0)
    LookUp_Crown_site[index_crown] = extent + extent_full * extent;      // this is the label of the site at the center of the crown ( x = extent, y = extent)
    
    /* loop over crown */
    for(int col = 0; col < extent_full; col++){
        for(int row = 0; row < extent_full; row++){
            xx = col - extent;                                      // distance from center (x = extent) in x direction
            yy = row - extent;                                      // distance from center (y = extent) in y direction
            if(!(xx == 0 & yy == 0)){
                site_rel = col + extent_full * row;
                dist = xx*xx + yy*yy;
                /* now order the arrays according to distance from center */
                /* index_crown saves last filled position in array */
                /* for every voxel we run through the array from position zero to last filled position, and check where to put the new value */
                for(int i=0; i < index_crown + 1; i++){
                    int temp = Crown_dist[i];
                    int site_rel_temp = LookUp_Crown_site[i];
                    if(dist <= temp){                               // if distance is smaller than at current array position, move everything up
                        Crown_dist[i] = dist;
                        LookUp_Crown_site[i] = site_rel;
                        dist = temp;
                        site_rel = site_rel_temp;
                    }
                }
                Crown_dist[index_crown+1] = dist;                   // the last value that has been pushed into the dist variable fills a new position
                index_crown = index_crown + 1;
            }
        }
    }
    
#endif
    
    In.open(inputfile, ios::in);
    if(!mpi_rank) {
        do{
            In.getline(buffer,256,'\n');
            out << buffer <<endl;
        }while (!In.eof()) ;
    }
    
    In.close();
    /* Close ifstream In */
    
     if(_DISTURBANCE || _LOGGING){							// new sylviculture module parameters
    	fstream In(inputfile_sylviculture, ios::in);
    
    	/*** Initialization of the sylviculture parametres ***/
    	/***************************************************/
    
    	if (In) {
        	for(int ligne=0;ligne<4;ligne++) In.getline(buffer,128,'\n');
        	/* general parameters */
        	In >> disturb_iter; In.getline(buffer,128,'\n');
       		In.getline(buffer,128,'\n');
       		/* disturbance parameters */
        	In >> disturb_intensity; In.getline(buffer,128,'\n');
        	In.getline(buffer,128,'\n');
        	/* logging parameters */
        	In >> designated_volume; In.getline(buffer,128,'\n');
        	designated_volume *= cols*NV*rows*NH/10000;
        	In >> harvested_volume; In.getline(buffer,128,'\n');
        	harvested_volume *= cols*NV*rows*NH/10000;
        	In >> numespharvestable; In.getline(buffer,128,'\n');
        	/* species parameters */
        	for(int jump=0;jump<3;jump++) In.getline(buffer,128,'\n');                           /* Read species parameters (ifstream In) */
    		for(sp=1;sp<=numespharvestable;sp++) {
        		InitialiseLogging(In);
    		}
    	}
    	else {
        	cout<< "ERROR with the input file of sylviculture parameters" << endl;
    	}
    	In.close();

    	/*** Sylviculture parameters in par file ***/
    	/******************************************/
    	In.open(inputfile_sylviculture, ios::in);
        	do{
            	In.getline(buffer,256,'\n');
            	out3 << buffer <<endl;
        	}while (!In.eof()) ;

    	In.close();
	}
    
    
       
    /*** Initialization of output streams ***/
    /****************************************/
    
    char nnn[200];
    if(!mpi_rank) {
        if(_OUTPUT_reduced){
            sprintf(nnn,"%s_%i_outputs.txt",buf, easympi_rank);
            output[0].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_final_pattern.txt",buf, easympi_rank);
            output[10].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_final_pattern_10cm.txt",buf, easympi_rank);
            output[18].open(nnn, ios::out);

            if(_OUTPUT_fullFinal){
        		sprintf(nnn, "%s_%i_fullfinal.txt", buf, easympi_rank);
        		output[37].open(nnn, ios::out);
            } // NINO
        }
        else{
            sprintf(nnn,"%s_%i_abund.txt",buf, easympi_rank);
            output[0].open(nnn, ios::out);
            if (!output[0]) {
                cout<< "ERROR with abund file"<< endl;
            }
            sprintf(nnn,"%s_%i_abu10.txt",buf, easympi_rank);
            output[1].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_abu30.txt",buf, easympi_rank);
            output[2].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_ba.txt",buf, easympi_rank);
            output[3].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_npp.txt",buf, easympi_rank);
            output[4].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_gpp.txt",buf, easympi_rank);
            output[5].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_ba10.txt",buf, easympi_rank);
            output[6].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_ppfd0.txt",buf, easympi_rank);
            output[7].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_death.txt",buf, easympi_rank);
            output[8].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_state.txt",buf, easympi_rank);
            output[9].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_final_pattern.txt",buf, easympi_rank);
            output[10].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_site1.txt",buf, easympi_rank);
            
            // used to be added to explore potential link with Belassen curve, could be suppressed, but maybe useful to have an idea of the magnitude and distribution of increment of dbh
            
            output[11].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_site2.txt",buf, easympi_rank);
            output[12].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_site3.txt",buf, easympi_rank);
            output[13].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_site4.txt",buf, easympi_rank);
            output[14].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_site5.txt",buf, easympi_rank);
            output[15].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_site6.txt",buf, easympi_rank);
            output[16].open(nnn, ios::out);
            
            //output[17] for parameter space is valid both for reduced and standard output, defined below
            
            sprintf(nnn,"%s_%i_sp_par.txt",buf, easympi_rank);
            output[18].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_agb.txt",buf, easympi_rank);
            output[19].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_Rday.txt",buf, easympi_rank);
            output[20].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_Rnight.txt",buf, easympi_rank);
            output[21].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_Rstem.txt",buf, easympi_rank);
            output[22].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_death1.txt",buf, easympi_rank);
            output[23].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_death2.txt",buf, easympi_rank);
            output[24].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_death3.txt",buf, easympi_rank);
            output[25].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_deathrate.txt",buf, easympi_rank);
            output[26].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_litterfall.txt",buf, easympi_rank);
            output[27].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_leafdens1.txt",buf, easympi_rank);
            output[28].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_leafdens2.txt",buf, easympi_rank);
            output[29].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_leafdens3.txt",buf, easympi_rank);
            output[30].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_dbh.txt",buf, easympi_rank);
            output[31].open(nnn,ios::out);
            sprintf(nnn,"%s_%i_vertd.txt",buf, easympi_rank);
            output[32].open(nnn,ios::out);
            sprintf(nnn,"%s_%i_100yearsofsolitude.txt",buf, easympi_rank);
            output[33].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_NDDfield.txt",buf, easympi_rank);
            output[34].open(nnn, ios::out);
            sprintf(nnn,"%s_%i_cica.txt",buf, easympi_rank);
            output[35].open(nnn, ios::out);
            // output[36] to register killed trees during a disturbance event
            if(_DISTURBANCE || _LOGGING){
				sprintf(nnn,"%s_%i_disturbance.txt",buf, easympi_rank);
            	output[36].open(nnn, ios::out);

            }  // NINO
            if(_OUTPUT_fullFinal){
        	sprintf(nnn, "%s_%i_fullfinal.txt", buf, easympi_rank);
        	output[37].open(nnn, ios::out);
             } // NINO
        }
    
//            sprintf(nnn,"%s_%i_paramspace.txt",buf, easympi_rank);
//            output[17].open(nnn, ios::out);
//            
//            output[17] << "proc"  << "\t" << easympi_rank << endl;
//            output[17] << "phi"  << "\t" << phi << endl;
//            output[17] << "k"  << "\t" <<klight << endl;
//            output[17] << "fallocwood"  << "\t" << fallocwood << endl;
//            output[17] << "falloccanopy"  << "\t" << falloccanopy << endl;
//            output[17] << "m"  << "\t" << m << endl;
//            output[17] << "m1"  << "\t" << m1 << endl;
        
    }
}

/*****************************************
 *** Initialisation logging parameters *** NINO
 ****************************************/

void InitialiseLogging(fstream& is){

	char name[256];
	int sp, interest;
	float dbhmin, dbhmax;
	
    /*** Read parameters ***/
    is  >> name >> sp >> dbhmin >> dbhmax >> interest;

    /* Adding parameters to the species */
    S[sp].s_harvestable=1;
    S[sp].s_dbhmin=dbhmin;
    S[sp].s_dbhmax=dbhmax;
    S[sp].s_interest=interest;        
}

/***************************************
 *** Initialisation from inventories ***
 ***************************************/


void InitialiseFromData(){


    // prepare data set beforehand: col/x_dim, row/y_dim, dbh_measured, species_label from inventories
    // or use fullfinal.txt output to load TROLL simulation in its full state
    
    fstream In(inputfile_data, ios::in);  
                                                  // input stream
    float col_data, row_data, dbh_measured, sp_lab_data;                            // values needed for tree initialisation
    
// New full-final input : the number of columns in the input files depends on the activated macros 
// We assume that the user wil always run simulations from fullfinal with the same activated macros as those of the simulation that resulted in the input data. 
// you need to use it in another way, you have to adapt the output/input routines. 
// Defining the number of parameters
	
	int nbparameters;

    if(_INPUT_fullFinal){
    	

		if(_BASICTREEFALL || _TREEFALL){
			#if defined(INTRASPECIFIC) && defined(CROWN_SHAPE)
				nbparameters = 29;
			#elif defined(INTRASPECIFIC) && !defined(CROWN_SHAPE)
				nbparameters = 25;
			#else // corresponds to !defined(INTRASPECIFIC) && !defined(CROWN_SHAPE) because I don't use crown_shape without intraspecific
				nbparameters = 17;
			#endif

		}
		else{
			#if defined(INTRASPECIFIC) && defined(CROWN_SHAPE)
				nbparameters = 28;
			#elif defined(INTRASPECIFIC) && !defined(CROWN_SHAPE)
				nbparameters = 24;
			#else // !defined(INTRASPECIFIC) && !defined(CROWN_SHAPE)
				nbparameters = 16;
			#endif
		}
    }
    /* Defining an array than will contain every parameter for a given tree, except species_lab, coordinates, and dbh. */
    /* This is the main argument of the new Tre::BirthFromSimulation() function*/

    float fullparameters[nbparameters];

    int col_int, row_int, data_read=0, data_initialised=0; // fullfinal=0;                      // diagnostics
    float height_max=0;                                                             // diagnostics
    
    nblivetrees=0;
     
    cout << "Reading from file " << inputfile_data << "\n";
    
    In.getline(buffer,256,'\n');    
    string header;
    getline(In, header);
    cout << "Header line skipped \n";

    if(_INPUT_fullFinal){ // Improved the previous Dirty Hack (NINO)
    	cout << "Full final will be used to load entirely the previous TROLL modelled forest." << endl;
    	cout << "Header size" << header.size() << endl;
    }

    //if(header.size() < 100){ // ! DIRTY HACK NEED TO BE IMPROVED !.    	
    //	cerr << "Your data input is apparently not a Full final pattern from a previous simulation" << endl;
    //}
                                                // skip header lines
    

    // NINO : new version of the data reading procedure

    string line; // The currently read line
    int colA;
    float x; // The parameter index, varying betwen 0 and nbparameters-1 ; the variable-receiver.
    
    while (!In.eof() && getline(In,line) && data_read < sites){       // restricting to data sets with a maximum number of values corresponding to no. of sites
    
    	std::istringstream streamA(line);
    	streamA >> col_data >> row_data >> dbh_measured  >> sp_lab_data;
        	
			if(_INPUT_fullFinal){ 
				colA = 0;
				while(streamA >> x){ // && colA > nbparameters){
					fullparameters[colA] = x;
					colA++;
				}
			}
       // In.getline(buffer, 256, '\n'); // reads additional information into buffer
        if((sp_lab_data > 0) && (sp_lab_data <= numesp) && (col_data >= 0) && (col_data < cols) && (row_data >= 0) && (row_data < rows)){
            
            // read only species that are rendered in TROLL (species label to be added via R / comparison input file and data)
            // cout << "col: " << round(col_data) << " row: " << round(row_data) << " species: " << sp_lab_data << " dbh: " << dbh_measured << " data measured \n";
            dbh_measured = 0.001*dbh_measured;          //here given in mm, converting to m
            col_int = (int) (col_data+0.5f);            //rounding, works since negatives have been eliminated before
            row_int = (int) (row_data+0.5f);
            
            // immediate tree birth
            
            if(T[col_int+row_int*cols].t_age==0){

            	 if(_INPUT_fullFinal){
            	 	T[col_int+row_int*cols].BirthFromSimulation(S,sp_lab_data,col_int+row_int*cols,dbh_measured, fullparameters); 
            	 }
            	 else{
            		T[col_int+row_int*cols].BirthFromData(S,sp_lab_data,col_int+row_int*cols,dbh_measured);            	 	
            	 }
            	
            }
            else{
            	cout << "You asshole have duplicates in your tree grid at coordinates " << "X;" << col_int <<" ; Y:" << row_int << endl;	
            } 
            
            if(height_max<T[col_int+row_int*cols].t_Tree_Height){
            	height_max = T[col_int+row_int*cols].t_Tree_Height;
            } 
            
            // first attempt: simple, only trees with coordinates, only known species
            // other possibilities: not spatially explicit and/or assign species randomnly to trees whose species are not known
            
            data_initialised++;
        }
        data_read++;
    }
    
    cout << "\n" << data_read << " rows read from file. " << data_initialised << " rows usable for initialisation from data. \n";
    cout << "Maximum height of trees included is: " << height_max << "\n";
    cout << "NBtrees from Data:\t" << nblivetrees << "\n";
    cout << "Initialisation from data finished \n";
    
    In.close();
}



/***************************************
 *** Field dynamic memory allocation ***
 ***************************************/

void AllocMem() {
    
    int spp,h,i; // haut became h
    
    // HEIGHT = 80;
    // instead of static definition of height, calculate according to maximum height that can be reached by trees
#ifdef INTRASPECIFIC
    HEIGHT = 80;
#else
    HEIGHT = 0;
    
    for(spp=1;spp<=numesp;spp++) {
        HEIGHT = max(HEIGHT,int(S[spp].s_hmax*S[spp].s_dmax*1.5/(S[spp].s_dmax*1.5+S[spp].s_ah)));  //   in number of vertical cells
    }
#endif
    
    cout << "HEIGHT " << HEIGHT << endl;
    
    float d,r;
    d = 0.0;
    r = 0.0;
    for(spp=1;spp<=numesp;spp++){
        d = maxf(d,S[spp].s_dmax*1.5);
        /* in number of horizontal cells */
        r = maxf(r,ra0+S[spp].s_dmax*1.5*ra1);
        /* in number of horizontal cells */
    }
    
    RMAX = int(r+p_nonvert*NH*LV*HEIGHT);
    //  RMAX = int(r);
    SBORD = cols*RMAX;
    dbhmaxincm = int(100.*d);
    
    if(!mpi_rank) {
        /*cout << "HEIGHT : " << HEIGHT
         << " RMAX : " << RMAX << " DBH : " << DBH <<"\n";
         cout.flush();*/
        
        if(RMAX>rows){
            /* Consistency tests */
            cerr << "Error : RMAX > rows \n";
            exit(-1);
        }
        if(HEIGHT > rows){
            cerr << "Error : HEIGHT > rows \n";
            exit(-1);
        }
    }
    
    
    /*** Initialization of dynamic Fields ***/
    /****************************************/
    
    if (NULL==(nbdbh=new int[dbhmaxincm])) cerr<<"!!! Mem_Alloc\n";                         /* Field for DBH histogram */
    if (NULL==(layer=new float[HEIGHT+1])) cerr<<"!!! Mem_Alloc\n";                         /* Field for variables averaged by vertical layer */
    if (NULL==(SPECIES_GERM=new int[numesp+1])) cerr<<"!!! Mem_Alloc\n";                           /* Field for democratic seed germination */
    if(_SEEDTRADEOFF){
        if (NULL==(PROB_S=new float[numesp+1])) cerr<<"!!! Mem_Alloc\n";
    }
    if(_NDD){
        if (NULL==(PROB_S=new float[numesp+1])) cerr<<"!!! Mem_Alloc\n";
    }
    //  if (NULL==(persist=new long int[nbiter])) cerr<<"!!! Mem_Alloc\n";                  /* Field for persistence */
    //  if (NULL==(distr=new int[cols])) cerr<<"!!! Mem_Alloc\n";
    
    if (NULL==(LAI3D=new float*[HEIGHT+1]))                                                   /* Field 3D */
        cerr<<"!!! Mem_Alloc\n";                                                            /* Trees at the border of the simulated forest need to know the canopy occupancy by trees in the neighboring processor.*/
    for(h=0;h<(HEIGHT+1);h++)                                                          /* For each processor, we define a stripe above (labelled 0) and a stripe below (1). Each stripe is SBORD in width.*/
        if (NULL==(LAI3D[h]=new float[sites+2*SBORD]))                                   /* ALL the sites need to be updated.*/
            cerr<<"!!! Mem_Alloc\n";
    for(h=0;h<(HEIGHT+1);h++)
        for(int site=0;site<sites+2*SBORD;site++)
            LAI3D[h][site] = 0.0;
    
    if (NULL==(Thurt[0]=new unsigned short[3*sites]))                                       /* Field for treefall impacts */
        cerr<<"!!! Mem_Alloc\n";
    for(i=1;i<3;i++)
        if (NULL==(Thurt[i]=new unsigned short[sites]))
            cerr<<"!!! Mem_Alloc\n";
    if(_LOGGING)													/* global variable for logging module erased after use NINO :Maybe we should ad some ? */
    	if (NULL==(Tlogging[0]=new unsigned short[sites]))			/* Field for tree felling */
            cerr<<"!!! Mem_Alloc\n";
        if (NULL==(Tlogging[1]=new unsigned short[sites]))			/* Field for loging tracks */
            cerr<<"!!! Mem_Alloc\n";
		if (NULL==(Tlogging[2]=new unsigned short[sites]))			/* Field for gap damages */
            cerr<<"!!! Mem_Alloc\n";
    
    
#ifdef DCELL
    /* MAP_DCELL is a list of vectors -- there is one vector per dcell, and the vector has as many entries as the number of sites per dcell, and returns the absolute site value of these sites -- this is useful to distribute seeds across the dcell */
    int x0,y0,x,y,site0,site,dcol,drow;
    if (NULL==(MAP_DCELL=new int*[nbdcells])) cerr<<"!!! Mem_Alloc\n";
    for(int dcell=0;dcell<nbdcells;dcell++){
        if (NULL==(MAP_DCELL[dcell]=new int[sites_per_dcell])) cerr<<"!!! Mem_Alloc\n";
        for(x0=0;x0<length_dcell;x0++)
            for(y0=0;y0<length_dcell;y0++){
                site0=x0+y0*length_dcell;
                dcol=dcell%linear_nb_dcells;
                drow=dcell/linear_nb_dcells;
                x=dcol*length_dcell+x0;
                y=drow*length_dcell+y0;
                //cerr << "x_MAP " << x << "\ty_MAP " << y << "\t";
                site=x+y*cols;
                MAP_DCELL[dcell][site0] = site;
            }
    }
    if (NULL==(prior_DCELL=new double[sites_per_dcell])) cerr<<"!!! Mem_Alloc\n";
    for(int i=0;i<sites_per_dcell;i++) prior_DCELL[i]=1.0;
    if (NULL==(post_DCELL=new unsigned int[sites_per_dcell])) cerr<<"!!! Mem_Alloc\n";
    for(int i=0;i<sites_per_dcell;i++) post_DCELL[i]=0;
    
    if (NULL==(prior_GERM=new double[numesp+1])) cerr<<"!!! Mem_Alloc\n";
    for(int i=0;i<=numesp;i++) prior_GERM[i]=0.0;
    if (NULL==(post_GERM=new unsigned int[numesp+1])) cerr<<"!!! Mem_Alloc\n";
    for(int i=0;i<=numesp;i++) post_GERM[i]=0;
    
#endif
    
#ifdef MPI                                                                                      /* Fields for MPI operations */
    for(i=0;i<2;i++){                                                                       /*  Two fields: one for the CL north (0) one for the CL south (1) */
        if (NULL==(LAIc[i]=new unsigned short*[HEIGHT+1]))                                    /* These fields contain the light info in the neighboring procs (2*SBORD in width, not SBORD !). They are used to update local fields */
            cerr<<"!!! Mem_Alloc\n";
        for(h=0;h<(HEIGHT+1);h++)
            if (NULL==(LAIc[i][h]=new unsigned short[2*SBORD]))
                cerr<<"!!! Mem_Alloc\n";
    }
#endif
}


/***************************************
 **** Initial non-local germination ****
 ***************************************/

void BirthInit() {
    
    if(!_FromData){
        nblivetrees=0;
    }
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    cout<<endl;
}


/*###############################################
 ###############################################
 #######  Evolution at each timestep    ########
 ###############################################
 ###############################################*/

void Evolution() {
    
    if(_DISTURBANCE) Disturbance();		                    /* Create a disturbance for a given iteration - NINO */
    if(_LOGGING) SelectiveLogging();		                /* Simulate selective logging for a given iteration - NINO*/
    UpdateField();                                              /* Update light fields and seed banks */
    if(_BASICTREEFALL || _TREEFALL) TriggerTreefallSecondary(); /* Compute and distribute Treefall events, caused by previous treefalls */
    UpdateTree();                                               /* Update trees */
    if(_BASICTREEFALL || _TREEFALL) TriggerTreefall();          /* Compute and distribute Treefall events */
    Average();                                                  /* Compute averages for outputs */
    if(!_OUTPUT_reduced) OutputField();                         /* Output the statistics */
}


//// BLOCK SYLVAIN NINO
/*##################################
 ####    Compute a disturbance   ###
 ##################################*/

void Disturbance() {

    if(iter == disturb_iter) {
        cout << "Disturbance of " << disturb_intensity * 100 << "% of BA." << endl;

        int site, row, col;
    	float dbh=0.0, disturb_dbh=0.0;

        for(site=0;site<sites;site++)
        	dbh += T[site].t_dbh;  

        while (disturb_dbh/dbh < disturb_intensity) {
        	int site=floor(genrand2()*sites);

        	if(T[site].t_age != 0) {
        		disturb_dbh += T[site].t_dbh;
        		// saving killed tree in disturbance.txt file
        		row = floor(site/cols);
        		col = site - (row*cols);
        		output[36] << "L" << "\t" << col << "\t" << row << "\t" << T[site].t_age << "\t" << T[site].t_dbh << "\t" << T[site].t_Tree_Height << "\t" << T[site].t_Crown_Radius << "\t" << T[site].t_Crown_Depth << "\t" << T[site].t_sp_lab << endl;
            	T[site].Death();
        	} 
        }
    }
}

/*#######################################
 ####    Simulate selective logging   ###
 #######################################*/

void SelectiveLogging() {

    if(iter == disturb_iter) {

    	int site;
    	for(site=0;site<sites;site++){
    		Tlogging[0][site]=0;		// tree felling
    		Tlogging[1][site]=0;		// tracks
    		Tlogging[2][site] = 0;			// gaps
    	}

    	cout << "###   Selective Logging   ###" << endl;
    	Designate();
    	Select();
    	Rot();
    	Fell();
    	MainTracks();
    	SecondaryTracks();
        cout << "### Selective Logging done ###" << endl;
    }

    if(iter == (disturb_iter+iterperyear)) {
    	GapDamages();  
    	int i;
    	for (i=0; i<3; i++) delete [] Tlogging[i];	 // free memory
    }
}

/********************************/
/* SelectiveLogging submodules */
/******************************/


void Designate() {

	//int site, col, row, sp, sph=0, designated;
	int site, sp, sph=0, designated;
	float volume, dbh_min[numespharvestable], min_dbh_min, max_dbh_max=0.0;

	/* getting species vector of minimum harvestable diameter */
	for(sp=1;sp<=numesp;sp++)
		if(S[sp].s_harvestable){
			dbh_min[sph]=S[sp].s_dbhmin; 
			sph++;
			if(S[sp].s_dbhmax > max_dbh_max)
				max_dbh_max = S[sp].s_dbhmax;
		}

	/* getting minimum value of minimum harvestable diameter among species*/
	min_dbh_min = dbh_min[0];
	for(sph=1;sph<numespharvestable;sph++)
		if(dbh_min[sph] < min_dbh_min)
			min_dbh_min = dbh_min[sph];

	/* designating tree, increasing minimum harvestable dbh if needed to br under the objective */
	//for(min_dbh_min; min_dbh_min < max_dbh_max; min_dbh_min += 0.1){
	//for(mindiam = min_dbh_min; mindiam < max_dbh_max; mindiam += 0.1){
	for(float increment = 0; increment < max_dbh_max - min_dbh_min; increment += 0.1){
		volume=0.0;
		designated=0;
		for(site=0;site<sites;site++){
        	if(T[site].t_age > 0										/*alive tree*/
        		&& S[T[site].t_sp_lab].s_harvestable 					/*harvestable species*/
        		&& T[site].t_dbh >= S[T[site].t_sp_lab].s_dbhmin+increment		/*reached minimum dbh*/
			  //&& T[site].t_dbh >= S[T[site].t_sp_lab].s_dbhmin
        		&& T[site].t_dbh <= S[T[site].t_sp_lab].s_dbhmax){		/*under maximum dbh*/
        		Tlogging[0][site] = 1;
        		volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        		designated++;
        	}
        }
        //if(volume < designated_volume)
        if(volume < designated_volume + 5 && volume > designated_volume - 5) // Nino
        	break;														/*if the volume is under the objective we can stop */
        //else
        	//for(sp=1;sp<=numesp;sp++)
        		//if(S[sp].s_harvestable)
        			//S[sp].s_dbhmin += 0.01;								/*if the volume is greater than the objective we need to derease minimum harvestable diameter for all species */
	}

	cout << designated << " trees have been designated, representing " << volume << " m3." << endl;
    cout << "dbh min is now " << min_dbh_min << endl;
}        
        

void Select() {
	
	//int site, sp, i, rank, rankmax=0, unselected=0;
	int site, sp, rank, rankmax=0, unselected=0;
	float volume=0.0;

	/* Calculating designated volume */
	for(site=0;site<sites;site++) 
		if(Tlogging[0][site] == 1)
			volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
    
	if(volume <= harvested_volume)
        cout << "All designated trees will be harvested." << endl;

	if(volume > harvested_volume){

		/* determining maximal interest rank (leat valuable species) */
		for(sp=1;sp<=numesp;sp++)
        	if(S[sp].s_interest > rankmax)
        		rankmax = S[sp].s_interest;	

        /* determining determining headcount for each rank */
        int rank_nb[rankmax];
        for(rank=0;rank<rankmax;rank++) rank_nb[rank]=0;
        for(site=0;site<sites;site++)
        	if(S[sp].s_harvestable)
        		rank_nb[S[T[site].t_sp_lab].s_interest]++;

        /* removing tree untill wanted volume is reached starting by highest rank */
        for(rank=rankmax-1;rank>=0;rank--){
        	while(rank_nb[rank]>0){
        		site=floor(genrand2()*sites);
        		if(Tlogging[0][site]==1 && S[T[site].t_sp_lab].s_interest==rank){
        			Tlogging[0][site]=0;
        			rank_nb[rank]--;
        			unselected++;
        			volume -= -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        			if(volume <= harvested_volume) break;
        		}
        	}
        	if(volume <= harvested_volume) break;
        }

        cout << unselected << " trees have been unselected, volume is now of " << volume << " m3." << endl;
	}
}


void Rot() {

	int site, rotten=0;
	float protten, volume=0.0;

	/* Calculating selected volume */
	for(site=0;site<sites;site++) 
		if(Tlogging[0][site] == 1)
			volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/

	/* evaluates each tree probability to be rotten, and remove it if randomly in the risk to be rotten*/
    for(site=0;site<sites;site++){
       	if(Tlogging[0][site]==1){
       		protten = 1 / (1 + exp(-(-5.151 + 0.042*T[site].t_dbh*100))); /*Probability to be rotten*/
       		if(genrand2() < protten){
       			Tlogging[0][site]=0;
       			rotten++;
           		volume -= -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
       		}
       	}
       } 
    cout << rotten << " trees are rotten, volume is now of " << volume << " m3." << endl;
}


void Fell() {

	int site, row, col, felled=0;
	float volume=0.0;

	/* fell the selected tree not rotten */
    for(site=0;site<sites;site++){
        if(Tlogging[0][site]==1){
        	row = floor(site/cols);
        	col = site-(row*cols);
        	volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        	output[36] << "L" << "\t" << col << "\t" << row << "\t" << T[site].t_age << "\t" << T[site].t_dbh << "\t" << T[site].t_Tree_Height << "\t" << T[site].t_Crown_Radius << "\t" << T[site].t_Crown_Depth << "\t" << T[site].t_sp_lab << endl;
           	T[site].FellTree();
           	felled ++;
        }
    } 
    cout << felled << " trees have been felled representing " << volume << " m3." << endl;
}


void MainTracks() {

    int site, row, col, individuals=0;
    float volume=0.0;
           
    for(row=0;row<(rows/2);row++){
        for(col=((cols/2)-3);col<((cols/2)+3);col++){
        	site = col+row*cols;
        	Tlogging[1][site] = 1;
        	if(T[site].t_age != 0) {
        		volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        		output[36] << "MT" << "\t" << col << "\t" << row << "\t" << T[site].t_age << "\t" << T[site].t_dbh << "\t" << T[site].t_Tree_Height << "\t" << T[site].t_Crown_Radius << "\t" << T[site].t_Crown_Depth << "\t" << T[site].t_sp_lab << endl;
            	T[site].Death();
            	individuals++;
        	}
       	}
    }    
    cout << individuals << " trees have been killed for the main track representing " << volume << " m3." << endl;
}


void SecondaryTracks() {

	    int load[sites], tracks[sites], individuals=0, felt=0;
        int site, col, row, site0, row0, col0, siteT, rowT, colT;
        float d, d0, volume=0.0;

        /* Counting number of felt trees to skid */
        for(site=0;site<sites;site++) felt += Tlogging[0][site];
        
        while(felt > 0){

        	/*Computing loadings and tracks distance for each tree*/
        	for(site0=0;site0<sites;site0++){ 
        		load[site0]=0;
        		tracks[site0]=rows*rows + cols*cols;
        		row0 = floor(site0/cols);
        		col0 = site0-(row0*cols);
        		for(site=0;site<sites;site++){
        			if(Tlogging[0][site]==1 || Tlogging[1][site]==1){ // compute distance if the site is a felt tree or a track
        				row = floor(site/cols);
        				col = site-(row*cols);
        				d = (row - row0)*(row - row0) + (col - col0)*(col - col0);
        				if(Tlogging[0][site]==1 && d <= (30*30)) // site can evacuate the tree if is at a distance smaller than 30 meters
        					load[site0]++;
        				if(Tlogging[1][site]==1 && d < tracks[site0]) // save the track distance if it's closest than the previously saved one
        					tracks[site0]=d;
        			}
        		}
        	}

        	/*Seeking the best place to start the secondary track*/
        	site0=0;
        	for(site=0;site<sites;site++){
        		if(load[site]>load[site0]) // best candidate is the one which can evacuate maximum number of trees
        			site0=site;
        		if(load[site]==load[site0] && tracks[site]<tracks[site0]) // for equal loadings, best candidate is the one with a minimum distance to join an existing track
        			site0=site;
        	}	

        	/*Seeking for the closest track*/
        	row0 = floor(site0/cols);
        	col0 = site0-(row0*cols);
        	d0 = rows*rows+cols*cols;
        	for(site=0;site<sites;site++){
        		if(Tlogging[1][site]==1){ // if it's a track compute distance to the track
        			rowT = floor(site/cols);
        			colT = site-(rowT*cols);
        			d = (row0 - rowT)*(row0 - rowT) + (col0 - colT)*(col0 - colT);
        			if(d<d0){ // if the track is closer than the previously saved one, keep the location
        				siteT=site;
        				d0=d;
        			}
        		}
        	}

        	/*Trace the secondary track*/
        	rowT = floor(siteT/cols);
        	colT = siteT-(rowT*cols);
        	do{
        		do {
            		for(int i=-2;i<=2;i++){ 
        				for(int j=-2;j<=2;j++){
        					site = (col0+i)+(row0+j)*cols;
        					if(site>=0 && site<sites) Tlogging[1][site]=1; //flag the track with a size of 4 meters
        				}
        			}
        			for(site=0;site<sites;site++){ 
        				if(Tlogging[0][site]==1){
        					row = floor(site/cols);
        					col = site-(row*cols);
        					d = (row - row0)*(row - row0) + (col - col0)*(col - col0);
        					if(d <= (33*33)){ //unflag served trees in a radius of 30 meters
        						Tlogging[0][site]=0;
        						felt--;
        					}
        				}
        			}
        			if(col0 > colT) col0--; //move in direction of the closest existing track
        			if(col0 < colT) col0++;
        			if(row0 > rowT) row0--;
        			if(row0 < rowT) row0++;
        		} while(row0 != rowT); //stop when we reach the closest existing track
        	} while(col0 != colT);
        	cout << "A secondary track have been traced, " << felt << " trees still need to be evacuated." << endl; //!LONG! computation, console output to follow advancement
		}

		/* Removing trees on secondary tracks */
        for(site=0;site<sites;site++){ 
        	if(Tlogging[1][site] == 1 && T[site].t_age != 0){
        		row = (site/cols);
        		col = site-(row*cols);
        		volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        		output[36] << "ST" << "\t" << col << "\t" << row << "\t" << T[site].t_age << "\t" << T[site].t_dbh << "\t" << T[site].t_Tree_Height << "\t" << T[site].t_Crown_Radius << "\t" << T[site].t_Crown_Depth << "\t" << T[site].t_sp_lab << endl;
            	T[site].Death();
            	individuals ++;
        	}
        }
    cout << individuals << " trees have been killed for secondary tracks representing " << volume << " m3." << endl;
}


void GapDamages() {

	cout << "###   Selective Logging Long Term Damages   ###" << endl;
	int site, row, col, siteG, rowG, colG;
    float deathrate, gaps_deathrate, gaps_hurt, d, dgaps[sites];

    /*Initialise dgaps to maximum distance to a gaps (null effect)*/
    for(site=0;site<sites;site++) dgaps[site] = rows*rows + cols*cols;

    /*Compute for each tree the distance to the closest gap*/
    for(siteG=0;siteG<sites;siteG++){
        if(Tlogging[2][siteG] == 1){
        	if(T[siteG].t_age > 1){ // New trees could have been recruited over a year
        		rowG = floor(siteG/cols);
				colG = siteG-(rowG*cols);
				for(site=0;site<sites;site++){
					row = floor(site/cols);
        			col = site-(row*cols);
         			d = (row - rowG)*(row - rowG) + (col - colG)*(col - colG);
         			if(d < dgaps[site])	dgaps[site] = d;
         		}
			}
		}
    }

    /*Hurt trees depending on their distance to a gaps following an allometry fitted with Paracou data*/
    for(site=0;site<sites;site++){
       	if(T[site].t_age != 0 && T[site].t_dbh > 0.1){ //tree with dbh<10 have not an increased mortality closed to gaps, on the contrary they'll have a tendency 
       		gaps_deathrate = -4.441 + 0.762*exp(0.064*sqrt(dgaps[site]));
       		gaps_deathrate = exp(gaps_deathrate) / (1 + exp(gaps_deathrate)); // Allometry representing gaps damages
       		deathrate = T[site].t_s->DeathRate(T[site].t_PPFD, T[site].t_dbh, T[site].t_NPPneg);
       		if(gaps_deathrate > deathrate){
       			gaps_hurt = T[site].t_Tree_Height/(2*(gaps_deathrate - deathrate));
       			T[site].t_hurt += gaps_hurt;
       		}
       	}
    }        
}
//////////// ENDBLOCK LOGGING NINO

/*##################################
 ####    Compute field LAI 3D    ###
 ####    Compute field Seed     ####
 ##################################*/

void UpdateField() {

    int site;
    int spp=0;
    
    /* set the iteration environment -- nb: the current structure of code suppose that environment is periodic (a period = a year), if one wants to input a variable climate, with interannual variation and climate change along the simulation, a full climatic input needs to be input (ie number of columns=iter and not iterperyear) and change iterperyear by nbiter here. */
    //CURRENTLY NOT USED: precip, WS, Wmean, e_s, e_a,VPDbasic,VPDday
    temp=Temperature[iter%iterperyear];
    tmax=DailyMaxTemperature[iter%iterperyear];
    tnight=NightTemperature[iter%iterperyear];
    precip=Rainfall[iter%iterperyear];
    WS=WindSpeed[iter%iterperyear];
    Wmax=MaxIrradiance[iter%iterperyear]*1.678;       // 1.678 is to convert irradiance from W/m2 to micromol of PAR /s /m2, the unit used in the FvCB model of photosynthesis
    Wmean=MeanIrradiance[iter%iterperyear];            // still in W/m2
    e_s=SaturatedVapourPressure[iter%iterperyear];
    e_a=VapourPressure[iter%iterperyear];
    VPDbasic=VapourPressureDeficit[iter%iterperyear];
    VPDday=DailyVapourPressureDeficit[iter%iterperyear];
    VPDmax=DailyMaxVapourPressureDeficit[iter%iterperyear];
    
    /***  Compute Field LAI3D  ***/
    /*****************************/
    
    
#ifdef MPI
    /* Reinitialize field LAI3D */
    for(i=0;i<2;i++)
        for(h=0;h<(HEIGHT+1);h++)
            for(site=0;site<2*SBORD;site++)
                LAIc[i][h][site] = 0;
#endif
    
    /* v.2.4.0: no changes, but proposition: move CalcLAI() further downwards? after T.[site].Birth()? Currently freshly born trees do not allocate their plant material to LAI3Dfield (but of course, only for the first iteration) */
    int sbsite;
    for(int h=0;h<(HEIGHT+1);h++)
        for(sbsite=0;sbsite<sites+2*SBORD;sbsite++)
            LAI3D[h][sbsite] = 0.0;
    
    for(site=0;site<sites;site++)                                    /* Each tree contribues to LAI3D */
        T[site].CalcLAI();
    
    for(int h=HEIGHT;h>0;h--){                                 /* LAI is computed by summing LAI from the canopy top to the ground */
        for(site=0;site<sites;site++){
            sbsite=site+SBORD;
            LAI3D[h-1][sbsite] += LAI3D[h][sbsite];
            if (LAI3D[h-1][sbsite] < 0) T[site].OutputTreeStandard();
        }
    }
    
#ifdef MPI
    /* Communicate border of field */
    /*MPI_ShareField(LAI3D,LAIc,2*SBORD);
     This MPI command no longer exists in openMPI
     Action 20/01/2016 TODO: FIX THIS */
    MPI_ShareField(LAI3D,LAIc,2*SBORD);
    for(int h=0;h<(HEIGHT+1);h++){
        /* Add border effects in local fields */
        if(mpi_rank)
            for(site=0;site<2*SBORD;site++)
                LAI3D[h][site] += LAIc[0][h][site];
        if(mpi_rank<mpi_size-1)
            for(site=0;site<2*SBORD;site++)
                LAI3D[h][site+sites] += LAIc[1][h][site];
    }
#endif
    
    /*** Evolution of the field Seed **/
    /*********************************/
    
    /* Pass seeds across processors => two more fields to be communicated between n.n. (nearest neighbor) processors.
     NB: dispersal distance is bounded by the value of 'rows'.
     At least 99 % of the seeds should be dispersed within the stripe or on the n.n. stripe.
     Hence rows > 4.7*max(dist_moy_dissemination),for an exponential dispersal kernel.*/
    
#ifdef DCELL
    if(iter%iterperyear == 0){
        for(spp=1;spp<=numesp;spp++) {                              /* External seed rain: constant flux from the metacommunity */
            for(int dcell=0;dcell<nbdcells;dcell++){
                // loop over dcells and add exactly s_nbext seeds to s_DCELL
                S[spp].s_DCELL[dcell]=S[spp].s_nbext;
                //if(spp<10) cerr << "BeforeDis, dcell#\t" << dcell << "\ts_DCELL[dcell]\t" << S[spp].s_DCELL[dcell] << "\tS[spp].s_nbext\t" << S[spp].s_nbext <<endl;
            }
        }
        for(site=0;site<sites;site++)                                       /* disperse seeds produced by mature trees */
            if(T[site].t_age)
                T[site].DisperseSeed();
#else
        for(site=0;site<sites;site++)                                       /* disperse seeds produced by mature trees */
            if(T[site].t_age)
                T[site].DisperseSeed();
#endif
#ifdef DCELL
    }
#endif
    
#ifdef MPI                                                              /* Seeds passed across processors */
    for(spp=1;spp<=numesp;spp++) {
        MPI_ShareSeed(S[spp].s_Gc,sites);
        S[spp].AddSeed();
    }
#endif
    
#ifdef DCELL  // This entire section is commented as it is no longer needed in DCELL
#else
    if(_SEEDTRADEOFF){
        if(!mpi_rank || S[spp].s_nbind*sites > 50000){
            for(spp=1;spp<=numesp;spp++) {                              /* External seed rain: constant flux from the metacommunity */
                for(int ii=0;ii<S[spp].s_nbext;ii++){
                    site = genrand2i()%sites;
                    S[spp].s_Seed[site]++;
                }
            }
        }
    }
    else {
        if(!mpi_rank || S[spp].s_nbind*sites > 50000){
            for(spp=1;spp<=numesp;spp++) {                              /* External seed rain: constant flux from the metacommunity */
                for(int ii=0;ii<S[spp].s_nbext;ii++){
                    site = genrand2i()%sites;
                    if(S[spp].s_Seed[site]!=1)
                        S[spp].s_Seed[site] = 1; /* check for optimization */
                }
            }
        }
    }
#endif
    
    
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if(_NDD){
        /*** Evolution of the field NDDfield **/
        /**************************************/
        
        float normBA=10000.0/(0.001+PI*R*R*BAtot);
        for(site=0;site<sites;site++) {
            
            for(spp=1;spp<=numesp;spp++) {
                /*if ((iter == int(nbiter-1))&&(site>80000)&&(site<85000))  {
                 sor[142]<< T[site].t_NDDfield[spp] << "\t" ;
                 }*/
                T[site].t_NDDfield[spp]=0;
            }
            //if (iter == int(nbiter-1))  sor[142]<< "\n";
            
            int row0=T[site].t_site/cols;
            int col0=T[site].t_site%cols;
            for(int col=max(0,col0-R);col<=min(cols-1,col0+R);col++) {
                for(int row=max(0,row0-R);row<=min(rows-1,row0+R);row++) {                   /* loop over the neighbourhood */
                    int xx=col0-col;
                    int yy=row0-row;
                    float d=sqrt(xx*xx+yy*yy);
                    if((d<=R)&&(d>0))  {                                             /* is the voxel within the neighbourhood?  */
                        int j=cols*row+col;
                        if(T[j].t_age) T[site].t_NDDfield[T[j].t_sp_lab]+= PI* T[j].t_dbh * T[j].t_dbh*0.25*normBA;
                        
                    }
                    
                    
                }
                
            }
        }
    }
}


/*##############################
 ####	Birth, Growth       ####
 ####	and death of trees  ####
 ##############################*/

void UpdateTree() {

    int site,spp;
    float flux;
    
#ifdef DCELL
    
    for(spp=1;spp<=numesp;spp++) S[spp].UpdateSeed();
    
    // A much simpler version of seed recruitment is implemented again based on the idea
    // of a multinomial sampling of the s_Seed field.
    
    for(site=0;site<sites;site++){
        int nbseeds=0;
        for(spp=1;spp<=numesp;spp++) nbseeds+=S[spp].s_Seed[site];
        if(T[site].t_age == 0 && nbseeds > 0){
            //int summ=0;
            //for(spp=1;spp<=numesp;spp++) summ+=S[spp].s_Seed[site];
            //cerr << site << "\t" << summ << endl;
            prior_GERM[0]=0.0;post_GERM[0]=0;
            for(spp=1;spp<=numesp;spp++){
                prior_GERM[spp]=double(S[spp].s_Seed[site]*S[spp].s_seedmass);
                post_GERM[spp]=0;
            }
            //for(spp=1;spp<=numesp;spp++) cerr << prior_GERM[spp] << " ";
            
            gsl_ran_multinomial(gslrand,numesp+1,1,prior_GERM,post_GERM);
            float sumprior=0.0,sumpost=0;
            for(int i=1;i<=numesp;i++){
                sumprior+=prior_GERM[i];sumpost+=post_GERM[i];
            }
            //cerr << numesp << "\t"<< sumprior << "\t"<< sumpost << "\n";
            
            //for(spp=1;spp<=numesp;spp++)
            //    cerr << spp << "\tpost_GERM\t" << post_GERM[spp] << endl;
            int selected_species=0;
            for(spp=1;spp<=numesp;spp++) if(post_GERM[spp]==1){
                selected_species=spp;
                break;}
            //cerr << "selected species:\t" << spp <<endl;
            flux = Wmax*exp(-flor(LAI3D[0][site+SBORD])*klight);
            if(flux>(S[selected_species].s_LCP))
                T[site].Birth(S,selected_species,site);
            /* If enough light, germination, initialization of NPP (LCP is the species light compensation point -- here, light is the sole environmental resources tested as a limiting factor for germination, but we should think about adding nutrients (N,P) and water conditions... */
            for(spp=1;spp<=numesp;spp++) S[spp].s_Seed[site]=0;
            
            //}
        }
    }
#else
    int iii;
    if(_SEEDTRADEOFF){
        for(site=0;site<sites;site++) {                                             /***** Local germination *****/
            if(T[site].t_age == 0) {
                iii=0;
                float tot=0;                                                        /* _SEEDTRADEOFF */
                for(spp=1;spp<=numesp;spp++){                                       /* lists all the species with a seed present at given site... */
                    if(S[spp].s_Seed[site]) {
                        SPECIES_GERM[iii]=spp;
                        if(_NDD){
                            float p=S[spp].s_Seed[site]*S[spp].s_seedmass/(T[site].t_NDDfield[spp]*10000+1);
                            PROB_S[iii]=tot+ p;
                            tot+=p;
                        }
                        else{
                            PROB_S[iii]=tot+ S[spp].s_Seed[site]*S[spp].s_seedmass;
                            tot+=S[spp].s_Seed[site]*S[spp].s_seedmass;
                        }
                        iii++;
                    }
                }
                if(iii) {                                                           /* ... and then randomly select one of these species */
                    double p=genrand2();                                    /* if SEEDTRADEOFF is defined, probability of species recruit are proportional to the number of seeds time the seed mass, if NDD is also defined the probablility is also inversly proportional to the species NDDfield term at that site */
                    float itot=1.0/tot;
                    int s=0;
                    while (p>PROB_S[s]*itot) {s++;}
                    spp=SPECIES_GERM[s];
                    flux = Wmax*exp(-flor(LAI3D[0][site+SBORD])*klight);
                    if(flux>(S[spp].s_LCP)){
                        /* If enough light, germination, initialization of NPP (LCP is the species light compensation point*/
                        /* here, light is the sole environmental resources tested as a limiting factor for germination, but we should think about adding nutrients (N,P) and water conditions... */
                        T[site].Birth(S,spp,site);
                    }
                }
            }
        }
    }
    
    else {
        
        for(site=0;site<sites;site++) {                                     /***** Local germination *****/
            if(T[site].t_age == 0) {
                iii=0;
                
                float tot=0;
                
                for(spp=1;spp<=numesp;spp++){                               /* lists all the species with a seed present at given site... */
                    if(S[spp].s_Seed[site]) {
                        SPECIES_GERM[iii]=spp;
                        
                        if(_NDD){
                            float p=1.0/(1.0+deltaR*T[site].t_NDDfield[spp]);
                            PROB_S[iii]=tot+ p;
                            tot+=p;
                        }
                        
                        iii++;
                    }
                }
                if(iii) {                                                   /* ... and then randomly select one of these species */
                    
                    if(_NDD){
                        double p=genrand2();                                    /* if SEEDTRADEOFF is define, probability of species recruit are proportional to the number of seeds time the seed mass, if NDD is  defines the probablility is also inversly proportional to the species NDDfield term at that site */
                        float itot=1/tot;
                        int s=0;
                        while (p>PROB_S[s]*itot) {s++;}
                        spp=SPECIES_GERM[s];
                    }
                    else {
                        
                        spp = SPECIES_GERM[rand()%iii];
                        
                    }
                    
                    /* otherwise all species with seeds present are equiprobable */
                    flux = Wmax*exp(-flor(LAI3D[0][site+SBORD])*klight);
                    if(flux>(S[spp].s_LCP)){
                        /* If enough light, germination, initialization of NPP (LCP is the species light compensation point*/
                        /* here, light is the sole environmental resources tested as a limiting factor for germination, but we should think about adding nutrients (N,P) and water conditions... */
                        T[site].Birth(S,spp,site);
                    }
                }
            }
            else{
                for(spp=1;spp<=numesp;spp++) S[spp].s_Seed[site]=0;
            }
        }
    }
#endif
    
    nbtrees_n10=nbdead_n1=nbdead_n10=nbdead_n30=0;

    for(site=0;site<sites;site++) {
        /***** Tree evolution: Growth or death *****/
        T[site].Update();
    }
    
#ifdef DCELL
    // This loop has been moved upwards before the tree birth loop
#else
    for(spp=1;spp<=numesp;spp++) S[spp].UpdateSeed();
#endif
    
}


/******************************
 *** Treefall gap formation *** //NINO : to check
 ******************************/

/* change in v.2.4: resetting Thurt[0] field is done in TriggerSecondaryTreefall() at the beginning of each iteration */
/* further changes: separation of _BASICTREEFALL and _TREEFALL and rewriting of Tree::FallTree() which is now Tree::Treefall(angle) */
/* t_hurt can now persist longer, so new treefall events are added to older damages (that, in turn are decaying) */
    
void TriggerTreefall(){
    int site;
    
    for(site=0;site<sites;site++)
        if(T[site].t_age) {
            float angle = 0.0, c_forceflex = 0.0;
            
            /* treefall is triggered given a certain flexural force */
            /* _BASICTREEFALL: just dependent on height + random uniform distribution */
            /* _TREEFALL: dependent on crown interaction, calculated deterministically by function Couple */
            
            if(_BASICTREEFALL){
                c_forceflex = genrand2()*T[site].t_Tree_Height;     // probability of treefall = 1-t_Ct/t_Tree_Height
                angle = float(twoPi*genrand2());                    // random angle
            }
            
            if(_TREEFALL) T[site].Couple(c_forceflex, angle);       // c_forceflex and angle are passed by reference and changed through function

            /* above a given stress threshold the tree falls */
            if(c_forceflex > T[site].t_Ct){
                T[site].Treefall(angle);
            }
        }
    
#ifdef MPI
    /* Treefall field passed to the n.n. procs */
    MPI_ShareTreefall(Thurt, sites);
#endif
    
    for(site=0;site<sites;site++){
    /* Update of Field hurt */
        if(T[site].t_age) {
            T[site].t_hurt +=Thurt[0][site+sites];                 // NEW in v.2.4: addition of damages, alternative: max()
            
            
#ifdef MPI
            if(mpi_rank) T[site].t_hurt = max(T[site].t_hurt,Thurt[1][site]);               // ? v.2.4: Update needed, Thurt[1], why max?
            if(mpi_rank<mpi_size-1) T[site].t_hurt = max(T[site].t_hurt,Thurt[2][site]);
#endif
        }
    }
}


/* NEW in v.2.4: TriggerSecondaryTreefall(), called at the beginning of each iteration */
/* translates damages from previous round into tree deaths, partly treefalls, partly removing them only (e.g. splintering) */
/* in the limit of p_tfsecondary = 0.0, this is equivalent to the previous computation */


void TriggerTreefallSecondary(){
    
    nbTreefall1 = 0;
    nbTreefall10 = 0;
    nbTreefall30 = 0;
    
    int site;
    
    for(site=0;site<sites;site++){
        Thurt[0][site] = Thurt[0][site+2*sites] = 0;
        Thurt[0][site+sites] = 0;
    }
    
    for(site=0;site<sites;site++){
        if(T[site].t_age){
            if(2.0*T[site].t_hurt*genrand2() > T[site].t_Tree_Height) {        // check whether tree dies
                if(p_tfsecondary > genrand2()){                              // check whether tree falls or dies otherwise
                   float angle = float(twoPi*genrand2());                    // random angle
                   T[site].Treefall(angle); // NINO : angle could be used to do a directional treefall
                } else {
                    T[site].Death();
                }
            } else {
                T[site].t_hurt *= hurt_decay;                                // reduction of t_hurt according to hurt_decay, could be moved to Tree::Growth() function and made dependent on the tree's carbon gain
            }
        }
    }
    
#ifdef MPI
    /* Treefall field passed to the n.n. procs */
    MPI_ShareTreefall(Thurt, sites);
#endif
}

    

/*###############################################
 ###############################################
 #######        Output routines         ########
 ###############################################
 ###############################################*/

/*********************************************************
 *** Calculation of the global averages every timestep ***
 *********************************************************/

void Average(void){
    
    int site,spp,i;
    float sum1=0.0,sum10=0.0,sum30=0.0,ba=0.0,npp=0.0,gpp=0.0, ba10=0.0, agb=0.0, rday=0.0, rnight=0.0, rstem=0.0, litterfall=0.0;
    
    if(!mpi_rank) {
        for(spp=1;spp<=numesp;spp++)
            for(i=0;i<12;i++)
                S[spp].s_output_field[i]=0;
        
        float inbcells = 1.0/float(sites*mpi_size);
        float inbhectares = inbcells*NH*NH*10000.0;
        
        if(_OUTPUT_reduced){
            output[0] << iter << "\t";
            for(spp=1;spp<=numesp;spp++)
                output[0] << float(S[spp].s_nbind)*inbhectares << "\t";
        }
        else{
            for(i=0;i<7;i++) output[i] << iter << "\t";
            for(i=20;i<23;i++) output[i] << iter << "\t";
            for(spp=1;spp<=numesp;spp++) output[0] << float(S[spp].s_nbind)*inbhectares << "\t";
        }
        for(site=0;site<sites;site++)T[site].Average();
        
        
        for(spp=1;spp<=numesp;spp++) {
            S[spp].s_output_field[1] *= inbhectares;
            //species number of trees with dbh>10cm
            S[spp].s_output_field[2] *= inbhectares;
            //species number of trees with dbh>30cm
            S[spp].s_output_field[3] *= 3.1415*0.25*inbhectares;
            //species basal area
            S[spp].s_output_field[4] *= inbhectares;
            //species total NPP (sum of t_NPP) in MgC (per timestep)
            S[spp].s_output_field[5] *= inbhectares;
            //species total GPP (sum of t_GPP) in MgC (per timestep)
            S[spp].s_output_field[6] *= 3.1415*0.25*inbhectares;
            //species basal area; with only trees with dbh>10cm
            S[spp].s_output_field[7] *= inbhectares;
            //species aboveground biomass
            S[spp].s_output_field[8] *= inbhectares;
            /* species leaf Rday in MgC (per timestep) */
            S[spp].s_output_field[9] *= inbhectares;
            //species leaf Rnight in MgC (per timestep)
            S[spp].s_output_field[10] *= inbhectares;
            //species Rstem  in MgC (per timestep)
            S[spp].s_output_field[11] *= inbhectares;
            //species litterfall  in MgC (per timestep)
            sum1+= float(S[spp].s_nbind)*inbhectares;
            sum10 += S[spp].s_output_field[1];
            sum30 += S[spp].s_output_field[2];
            ba += S[spp].s_output_field[3];
            npp += S[spp].s_output_field[4];
            gpp += S[spp].s_output_field[5];
            ba10 += S[spp].s_output_field[6];
            agb += S[spp].s_output_field[7];
            rday += S[spp].s_output_field[8];
            rnight += S[spp].s_output_field[9];
            rstem += S[spp].s_output_field[10];
            litterfall += S[spp].s_output_field[11];
            
            if(!_OUTPUT_reduced){
                for(i=1;i<7;i++) output[i] << S[spp].s_output_field[i] << "\t";
                output[19] << S[spp].s_output_field[7] << "\t";
                output[20] << S[spp].s_output_field[8] << "\t";
                output[21] << S[spp].s_output_field[9] << "\t";
                output[22] << S[spp].s_output_field[10] << "\t";
                output[27] << S[spp].s_output_field[11] << "\t";
            }
        }
        
        if (_NDD) {
            BAtot=ba;
        }
        
        cout << iter << "\tNBtrees\t"<<nblivetrees << "\t(" << sum1 << " | " << sum10 << " | " << sum30 << ") *** treefalls: " << nbTreefall1 << " | " << nbTreefall10 << " | " << nbTreefall30 << endl;
        
        if(_OUTPUT_reduced){
            output[0] << sum1 << "\t";                                                     //total number of trees (dbh>1cm=DBH0)
            
            for(spp=1;spp<=numesp;spp++) output[0] << S[spp].s_output_field[1] << "\t";
            output[0] << sum10 << "\t";                                                    //total number of trees (dbh>10cm=DBH0)
            
            for(spp=1;spp<=numesp;spp++) output[0] << S[spp].s_output_field[2] << "\t";
            output[0] << sum30 << "\t";                                                    //total number of trees (dbh>30cm=DBH0)
            
            for(spp=1;spp<=numesp;spp++) output[0] << S[spp].s_output_field[6] << "\t";
            output[0] << ba10 << "\t";
            
            for(spp=1;spp<=numesp;spp++) output[0] << S[spp].s_output_field[4] << "\t";
            output[0] << npp << "\t";
            
            for(spp=1;spp<=numesp;spp++) output[0] << S[spp].s_output_field[5] << "\t";
            output[0] << gpp << "\t";
            
            for(spp=1;spp<=numesp;spp++) output[0] << S[spp].s_output_field[7] << "\t";
            output[0] << agb << endl;
        }
        else{
            
            output[0] << sum1 << endl;          //total number of trees (dbh>1cm=DBH0)
            output[1] << sum10 << endl;         //total number of trees with dbh>10cm
            output[2] << sum30 << endl;         //total number of trees with dbh>30cm
            output[3] << ba << endl;            //total basal area
            output[4] << npp << endl;           //total NPP in MgC per ha (per time step)
            output[5] << gpp << endl;           //total GPP in MgC par ha (per time step)
            output[6] << ba10 << endl;          //total basal area with only trees with dbh>10cm
            output[19] << agb << endl;          //total above ground biomass
            output[20] << rday << endl;         //total leaf day respiration
            output[21] << rnight << endl;       //total leaf night respiration
            output[22] << rstem << endl;        //total stem respiration
            output[27] << litterfall << endl;   //total litterfall
            
            float tototest=0.0, tototest2=0.0, flux;
            for(site=0;site<sites;site++) {
                flux = Wmax*exp(-flor(LAI3D[0][site+SBORD])*klight);
                tototest += flux;
                tototest2 += flux*flux;
            }
            tototest /=float(sites*LH*LH);                              // Average light flux (PPFD) on the ground
            tototest2 /=float(sites*LH*LH);
            if(iter)
                output[7] << iter<< "\tMean PPFDground\t"<< tototest << "\t" << sqrt(tototest2-tototest*tototest) << "\n";
            
            if(!_OUTPUT_reduced){
                if(_BASICTREEFALL) output[8] << iter << "\t" << nbdead_n1*inbhectares << "\t" << nbdead_n10*inbhectares<< "\t" << nbTreefall1*inbhectares << "\t" << nbTreefall10*inbhectares << "\t" << endl;
                else output[8] << iter << "\t" << nbdead_n1*inbhectares << "\t" << nbdead_n10*inbhectares<< "\t" << endl;
            }
        }
    }
    
#ifdef MPI
    MPI_Reduce(&(S[spp].s_nbind),&sind,1,
               MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(S[spp].s_output_field,S[spp].s_output_field,5,
               MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(Mortality,Mortality,4,
               MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&S[spp].s_output_field[6],&S[spp].s_output_field[6],5,
               MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
#endif
    
    if(!mpi_rank) {
        if(iter == 200){                                        // State at the 200th iteration (for all living trees: dbh, height, crown radius and depth and dbh increment)
            for(site=0;site<sites;site++) {
                if(T[site].t_dbh>0.0)
                    output[9] << T[site].t_dbh*LH*100 << "\t" << T[site].t_Tree_Height
                    << "\t" << T[site].t_Crown_Radius*LH
                    << "\t" << T[site].t_Crown_Depth*LV
                    << "\t" << T[site].t_ddbh*LH*100 << "\n";
            }
        }
    }
    
#ifdef MPI
    /* This section corresponds to the parallel version of
     the reporting of the global diagnostic variables. Since much work has been done
     on routine Average over the past years, this would need a full rewrite
     !!!!Action 20/01/2016: rework the parallel version of function Average!!!!
     
     MPI_Reduce(&(S[spp].s_nbind),&sind,1,
     MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(S[spp].s_output_field,S[spp].s_output_field,5,
     MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(Mortality,Mortality,4,
     MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&S[spp].s_output_field[6],&S[spp].s_output_field[6],5,
     MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
     */
#endif
    cout.flush();
}


/*********************************
 **** Output of Fields Shares ****
 *********************************/


void OutputField(){
    
    /* output for dbh histogram and mean LAI by height class */
    
    int site,h;
    
    if((nbout)&&((iter%freqout)==freqout-1)) {
        // output fields, nbout times during simulation (every freqout iterations)
        
        int d;
        for(d=0;d<dbhmaxincm;d++) nbdbh[d]=0;
        for(site=0;site<sites;site++) T[site].histdbh();
        
        for(h=0;h<(HEIGHT+1);h++){
            layer[h] = 0;
            for(site=0;site<sites;site++) layer[h] += LAI3D[h][site+SBORD];
        }
        
#ifdef MPI
        MPI_Status status;
        MPI_Reduce(nbdbh,nbdbh,dbhmaxincm,
                   MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
        
        MPI_Reduce(layer,layer,HEIGHT,
                   MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
#endif
        
        if(!mpi_rank) {
            
            // output of the dbh histograms (output[31])
            for(d=1;d<dbhmaxincm;d++) output[31] << d << "\t" << nbdbh[d]  << "\n";
            output[31] <<  "\n";
            
            // output of the mean LAI per height class (output[32])
            float norm = 1.0/float(sites*LH*LH*mpi_size);
            for(h=0;h<(HEIGHT+1);h++) output[32] << h*LV << "\t" << layer[h]*norm << "\n";
            output[32] <<  "\n";
        }
    }
}

/*********************************
 ******** Special Outputs ********
 *********************************/

/* This provides a snapshot of the forest whenever called. Per default, this function is used to create the final pattern */

void OutputSnapshot(fstream& output){
    for(int row=0;row<rows;row++)
        for(int col=0;col<cols;col++){
            output << col << "\t" << row << "\t" << T[col + cols*row].t_age << "\t" << T[col + cols*row].t_dbh << "\t" << T[col + cols*row].t_Tree_Height << "\t" << T[col + cols*row].t_Crown_Radius << "\t" << T[col + cols*row].t_Crown_Depth << "\t" << T[col + cols*row].t_sp_lab << endl;
        }
}

/* This provides a snapshot of all trees above 10cm in greater detail */

void OutputSnapshot10cm(fstream& output){
    for(int row=0;row<rows;row++){
        for(int col=0;col<cols;col++){
            if(T[col + cols*row].t_age > 0 && T[col + cols*row].t_dbh > 0.1){
                output << col << "\t" << row << "\t" << T[col + cols*row].t_age << "\t" << T[col + cols*row].t_dbh << "\t" << T[col + cols*row].t_Tree_Height << "\t" << T[col + cols*row].t_Crown_Radius << "\t" << T[col + cols*row].t_Crown_Depth << "\t" << T[col + cols*row].t_sp_lab << "\t" << T[col + cols*row].t_s->s_name << "\t" << "\t" <<
#ifdef INTRASPECIFIC
                    T[col + cols*row].t_LMA << "\t" << T[col + cols*row].t_Nmass << T[col + cols*row].t_Pmass  <<  "\t" << T[col + cols*row].t_wsg << "\t" << T[col + cols*row].t_dmax << "\t" << T[col + cols*row].t_intraspecific_multiplier_height << "\t" << T[col + cols*row].t_intraspecific_multiplier_CR << "\t" << T[col + cols*row].t_intraspecific_multiplier_CD << "\t" << T[col + cols*row].t_intraspecific_multiplier_LMA << "\t" << T[col + cols*row].t_intraspecific_multiplier_N << "\t" << T[col + cols*row].t_intraspecific_multiplier_P << "\t" << T[col + cols*row].t_intraspecific_deviation_wsg <<  "\t" <<  T[col + cols*row].t_intraspecific_multiplier_dmax <<
#else
                T[col + cols*row].t_s->s_LMA << "\t" << T[col + cols*row].t_s->s_Nmass << T[col + cols*row].t_s->s_Pmass  <<  "\t" << T[col + cols*row].t_s->s_wsg << "\t" << T[col + cols*row].t_s->s_dmax <<
#endif
                endl;
            }
        }
    }
}

/* This can be used to take snapshots of the forest in more detail and track its development over time. Currently used for 100 years development analysis. */

void OutputSnapshotDetail(fstream& output){
    for(int row=0;row<rows;row++)
        for(int col=0;col<cols;col++){
            output << iter << "\t" << col+cols*row << "\t" << col << "\t" << row << "\t" << T[col + cols*row].t_age << "\t" << T[col + cols*row].t_sp_lab << "\t" << T[col + cols*row].t_dbh << "\t" << T[col + cols*row].t_Tree_Height << "\t" << T[col + cols*row].t_Crown_Radius << "\t" << T[col + cols*row].t_Crown_Depth << "\t" << T[col + cols*row].t_leafarea << "\t" << T[col + cols*row].t_dens << "\t" << T[col + cols*row].t_GPP << "\t" << T[col + cols*row].t_hurt << endl;
        }
}

void OutputSnapshotFullFinal(fstream& output){ // To redefine here ? NINO

//	output << "\t" << "col" << "\t" << "row"<< "\t" << "dbh" << "\t" << "sp_lab" << "\t" << "NPPneg" << "\t" << "age" << "\t" << "hurt" << "\t" << "dbh_thresh" << "\t" << "hmax" << "\t" << "ah" << "\t" << "youngLA" << "\t" << "matureLA" << "\t" << "oldLA" << "\t" << "leafarea" << "\t" << "dens" << "\t" << "litter" << "\t" << "Tree_Height" << "\t" << "Crown_Radius" << "\t" << "Crown_Depth" << "\t" << "ddbh" << "\t" << "Ct" << "\t" << "intraspecific_multiplier_height" << "\t" << "intraspecific_multiplier_CR" << "\t" << "intraspecific_multiplier_CD" << "\t" << "intraspecific_multiplier_N" << "\t" << "intraspecific_multiplier_P" << "\t" << "intraspecific_multiplier_LMA" << "\t" << "intraspecific_deviation_wsg" << "\t" << "intraspecific_multiplier_dmax" << "\t" << "Crown_Slope_Top" << "\t" << "Crown_Slope_Bottom" << "\t" << "Crown_Volume" << "\t" << "Crown_Volume_layer" << endl;

	#if defined(INTRASPECIFIC) && defined(CROWN_SHAPE)
		output << "col" << "\t" << "row"<< "\t" << "dbh" << "\t" << "sp_lab" << "\t" << "NPPneg" << "\t" << "age" << "\t" << "hurt" << "\t" << "dbh_thresh" << "\t" << "hmax" << "\t" << "ah" << "\t" << "youngLA" << "\t" << "matureLA" << "\t" << "oldLA" << "\t" << "leafarea" << "\t" << "dens" << "\t" << "litter" << "\t" << "Tree_Height" << "\t" << "Crown_Radius" << "\t" << "Crown_Depth" << "\t" << "ddbh" << "\t" << "Ct" << "\t" << "intraspecific_multiplier_height" << "\t" << "intraspecific_multiplier_CR" << "\t" << "intraspecific_multiplier_CD" << "\t" << "intraspecific_multiplier_N" << "\t" << "intraspecific_multiplier_P" << "\t" << "intraspecific_multiplier_LMA" << "\t" << "intraspecific_deviation_wsg" << "\t" << "intraspecific_multiplier_dmax" << "\t" << "Crown_Slope_Top" << "\t" << "Crown_Slope_Bottom" << "\t" << "Crown_Volume" << "\t" << "Crown_Volume_layer" << endl;

		for(int row=0;row<rows;row++){
			for(int col=0;col<cols;col++){
				//if(!(_BASICTREEFALL || _TREEFALL)) T[col + cols*row].t_Ct = 0;

				output << col << "\t" << row << "\t" << T[col + cols*row].t_dbh << "\t" << T[col + cols*row].t_sp_lab << "\t" << T[col + cols*row].t_NPPneg << "\t" << T[col + cols*row].t_age << "\t" << T[col + cols*row].t_hurt << "\t" << T[col + cols*row].t_dbh_thresh << "\t" << T[col + cols*row].t_hmax << "\t" << T[col + cols*row].t_ah << "\t" << T[col + cols*row].t_youngLA << "\t" << T[col + cols*row].t_matureLA << "\t" << T[col + cols*row].t_oldLA << "\t" << T[col + cols*row].t_leafarea << "\t" << T[col + cols*row].t_dens << "\t" << T[col + cols*row].t_litter << "\t" << T[col + cols*row].t_Tree_Height << "\t" << T[col + cols*row].t_Crown_Radius << "\t" << T[col + cols*row].t_Crown_Depth << "\t" << T[col + cols*row].t_ddbh << "\t" << T[col + cols*row].t_Ct << "\t" << T[col + cols*row].t_intraspecific_multiplier_height << "\t" << T[col + cols*row].t_intraspecific_multiplier_CR << "\t" << T[col + cols*row].t_intraspecific_multiplier_CD << "\t" << T[col + cols*row].t_intraspecific_multiplier_N << "\t" << T[col + cols*row].t_intraspecific_multiplier_P << "\t" << T[col + cols*row].t_intraspecific_multiplier_LMA << "\t" << T[col + cols*row].t_intraspecific_deviation_wsg << "\t" << T[col + cols*row].t_intraspecific_multiplier_dmax << "\t" << T[col + cols*row].t_Crown_Slope_Top << "\t" << T[col + cols*row].t_Crown_Slope_Bottom << "\t" << T[col + cols*row].t_Crown_Volume << "\t" << T[col + cols*row].t_Crown_Volume_layer << endl;
			}
		}
        		
	#elif defined(INTRASPECIFIC) && !defined(CROWN_SHAPE)
		output << "col" << "\t" << "row"<< "\t" << "dbh" << "\t" << "sp_lab" << "\t" << "NPPneg" << "\t" << "age" << "\t" << "hurt" << "\t" << "dbh_thresh" << "\t" << "hmax" << "\t" << "ah" << "\t" << "youngLA" << "\t" << "matureLA" << "\t" << "oldLA" << "\t" << "leafarea" << "\t" << "dens" << "\t" << "litter" << "\t" << "Tree_Height" << "\t" << "Crown_Radius" << "\t" << "Crown_Depth" << "\t" << "ddbh" << "\t" << "Ct" << "\t" << "intraspecific_multiplier_height" << "\t" << "intraspecific_multiplier_CR" << "\t" << "intraspecific_multiplier_CD" << "\t" << "intraspecific_multiplier_N" << "\t" << "intraspecific_multiplier_P" << "\t" << "intraspecific_multiplier_LMA" << "\t" << "intraspecific_deviation_wsg" << "\t" << "intraspecific_multiplier_dmax" << endl;

		for(int row=0;row<rows;row++){
			for(int col=0;col<cols;col++){
				//if(!(_BASICTREEFALL || _TREEFALL)) T[col + cols*row].t_Ct = 0;
				output << col << "\t" << row << "\t" << T[col + cols*row].t_dbh << "\t" << T[col + cols*row].t_sp_lab << "\t" << T[col + cols*row].t_NPPneg << "\t" << T[col + cols*row].t_age << "\t" << T[col + cols*row].t_hurt << "\t" << T[col + cols*row].t_dbh_thresh << "\t" << T[col + cols*row].t_hmax << "\t" << T[col + cols*row].t_ah << "\t" << T[col + cols*row].t_youngLA << "\t" << T[col + cols*row].t_matureLA << "\t" << T[col + cols*row].t_oldLA << "\t" << T[col + cols*row].t_leafarea << "\t" << T[col + cols*row].t_dens << "\t" << T[col + cols*row].t_litter << "\t" << T[col + cols*row].t_Tree_Height << "\t" << T[col + cols*row].t_Crown_Radius << "\t" << T[col + cols*row].t_Crown_Depth << "\t" << T[col + cols*row].t_ddbh << "\t" << T[col + cols*row].t_Ct << "\t" << T[col + cols*row].t_intraspecific_multiplier_height << "\t" << T[col + cols*row].t_intraspecific_multiplier_CR << "\t" << T[col + cols*row].t_intraspecific_multiplier_CD << "\t" << T[col + cols*row].t_intraspecific_multiplier_N << "\t" << T[col + cols*row].t_intraspecific_multiplier_P << "\t" << T[col + cols*row].t_intraspecific_multiplier_LMA << "\t" << T[col + cols*row].t_intraspecific_deviation_wsg<< "\t" << T[col + cols*row].t_intraspecific_multiplier_dmax << endl;						
			}
		}
	#else // corresponds to !defined(INTRASPECIFIC) && !defined(CROWN_SHAPE) because I don't use crown_shape without intraspecific
		output << "col" << "\t" << "row"<< "\t" << "dbh" << "\t" << "sp_lab" << "\t" << "NPPneg" << "\t" << "age" << "\t" << "hurt" << "\t" << "dbh_thresh" << "\t" << "hmax" << "\t" << "ah" << "\t" << "youngLA" << "\t" << "matureLA" << "\t" << "oldLA" << "\t" << "leafarea" << "\t" << "dens" << "\t" << "litter" << "\t" << "Tree_Height" << "\t" << "Crown_Radius" << "\t" << "Crown_Depth" << "\t" << "ddbh" << "\t" << "Ct"  << endl;
		
		for(int row=0;row<rows;row++){
			for(int col=0;col<cols;col++){ // OK
				//if(!(_BASICTREEFALL || _TREEFALL)) T[col + cols*row].t_Ct = 0;
					output << col << "\t" << row << "\t" << T[col + cols*row].t_dbh << "\t" << T[col + cols*row].t_sp_lab << "\t" << T[col + cols*row].t_NPPneg << "\t" << T[col + cols*row].t_age << "\t" << T[col + cols*row].t_hurt << "\t" << T[col + cols*row].t_dbh_thresh << "\t" << T[col + cols*row].t_hmax << "\t" << T[col + cols*row].t_ah << "\t" << T[col + cols*row].t_youngLA << "\t" << T[col + cols*row].t_matureLA << "\t" << T[col + cols*row].t_oldLA << "\t" << T[col + cols*row].t_leafarea << "\t" << T[col + cols*row].t_dens << "\t" << T[col + cols*row].t_litter << "\t" << T[col + cols*row].t_Tree_Height << "\t" << T[col + cols*row].t_Crown_Radius << "\t" << T[col + cols*row].t_Crown_Depth << "\t" << T[col + cols*row].t_ddbh << "\t" << T[col + cols*row].t_Ct << endl;						
			}
		}
			#endif
}


// Basic BirthFromData: col, row, dbh_measured (in mm), species_label
// Tree attributes added: NPPneg, dbh_thresh, hmax, dbhmature, Tree_Height, Crown_Depth, Crown_Radius, ddbh, age, youngLA, matureLA, oldLA, leafarea, dens, litter, hurt


/* This provides relevant species parameters whenever called */

void OutputSpeciesParameters(fstream& output){
    for(int sp=1;sp<=numesp;sp++) output << S[sp].s_name << "\t" << S[sp].s_Nmass << "\t" << S[sp].s_Pmass << "\t" << S[sp].s_LMA << "\t" << S[sp].s_Vcmax << "\t" << S[sp].s_Jmax << "\t" << S[sp].s_Rdark << "\t" << S[sp].s_LCP << "\n";
}

/* This creates a Canopy Height Model and LAD profile whenever called */

void OutputLAI(fstream& output_CHM, fstream& output_LAI){
    for(int s=0;s<sites;s++){
        int height_canopy=0;
        for(int h=0;h<(HEIGHT+1);h++){
            if(LAI3D[h][s+SBORD] > 0.0) height_canopy = max(h,height_canopy);
        }
        output_CHM << s << "\t" << int(s/cols) << "\t" << int(s%cols) << "\t"  << height_canopy+1 << "\t" << LAI3D[0][s+SBORD] << endl;
    }

    float isites = 1.0/sites;
    for(int h=0;h<(HEIGHT);h++){
        float LAI3D_avg = 0.0;
        for(int s=0;s<sites;s++){
            if((LAI3D[h][s+SBORD]-LAI3D[h+1][s+SBORD]) < 0) cerr << "Be careful negative PAD!";
            LAI3D_avg+=LAI3D[h][s+SBORD]-LAI3D[h+1][s+SBORD];
        }
        output_LAI << h << "\t" << LAI3D_avg * isites << endl;
    }
}


/* This just writes the whole 3D LAI voxel field to file */

void OutputLAIFull(fstream& output_LAI3D){
    for(int s=0;s<sites;s++) {
        for(int h=0; h<HEIGHT;h++){
            output_LAI3D << s << "\t" << int(s/cols) << "\t" << int(s%cols) << "\t" << h << "\t" << LAI3D[h][s+SBORD] << endl;
        }
    }
}



/*********************************
 ********** MPI Routines *********
 *********************************/


#ifdef MPI

/* Communication of border fields in the parallel version of the code */
/* Only if the MPI option has been enabled */
void MPI_ShareSeed(unsigned char **c, int n) {
    
    MPI_Status status;
    
    if(p_rank == size-1)
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_CHAR,size-2,0,
                     c[3],n,MPI_UNSIGNED_CHAR,0,0,MPI_COMM_WORLD,&status);
    if(p_rank == 0)
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_CHAR,size-1,0,
                     c[3],n,MPI_UNSIGNED_CHAR,1,0,MPI_COMM_WORLD,&status);
    if((p_rank) && (p_rank < size-1))
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_CHAR,p_rank-1,0,
                     c[3],n,MPI_UNSIGNED_CHAR,p_rank+1,0,MPI_COMM_WORLD,&status);
    
    if(p_rank == 0)
        MPI_Sendrecv(c[1],n,MPI_UNSIGNED_CHAR,1,1,
                     c[2],n,MPI_UNSIGNED_CHAR,size-1,1,MPI_COMM_WORLD,&status);
    if(p_rank == size-1)
        MPI_Sendrecv(c[1],n,MPI_UNSIGNED_CHAR,0,1,
                     c[2],n,MPI_UNSIGNED_CHAR,size-2,1,MPI_COMM_WORLD,&status);
    if((p_rank) && (p_rank < size-1))
        MPI_Sendrecv(c[1],n,MPI_UNSIGNED_CHAR,p_rank+1,1,
                     c[2],n,MPI_UNSIGNED_CHAR,p_rank-1,1,MPI_COMM_WORLD,&status);
}

void MPI_ShareField(unsigned short **cl, unsigned short ***cp, int n) {
    
    MPI_Status status;
    for(int h=0;h<(HEIGHT+1);h++) {
        if(p_rank == 0)
            MPI_Sendrecv(cl[h],n,MPI_UNSIGNED_SHORT,size-1,h,
                         cp[1][h],n,MPI_UNSIGNED_SHORT,1,h,
                         MPI_COMM_WORLD,&status);
        if(p_rank == size-1)
            MPI_Sendrecv(cl[h],n,MPI_UNSIGNED_SHORT,size-2,h,
                         cp[1][h],n,MPI_UNSIGNED_SHORT,0,h,
                         MPI_COMM_WORLD,&status);
        if((p_rank) && (p_rank < size-1))
            MPI_Sendrecv(cl[h],n,MPI_UNSIGNED_SHORT,p_rank-1,h,
                         cp[1][h],n,MPI_UNSIGNED_SHORT,p_rank+1,h,
                         MPI_COMM_WORLD,&status);
        
        if(p_rank == 0)
            MPI_Sendrecv(cl[h]+sites,n,MPI_UNSIGNED_SHORT,1,h+HEIGHT,
                         cp[0][h],n,MPI_UNSIGNED_SHORT,size-1,h+HEIGHT,
                         MPI_COMM_WORLD,&status);
        if(p_rank == size-1)
            MPI_Sendrecv(cl[h]+sites,n,MPI_UNSIGNED_SHORT,0,h+HEIGHT,
                         cp[0][h],n,MPI_UNSIGNED_SHORT,size-2,h+HEIGHT,
                         MPI_COMM_WORLD,&status);
        if((p_rank) && (p_rank < size-1))
            MPI_Sendrecv(cl[h]+sites,n,MPI_UNSIGNED_SHORT,p_rank+1,h+HEIGHT,
                         cp[0][h],n,MPI_UNSIGNED_SHORT,p_rank-1,h+HEIGHT,
                         MPI_COMM_WORLD,&status);
    }
}

void MPI_ShareTreefall(unsigned short **c, int n) {
    
    MPI_Status status;
    if(p_rank == 0)
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_SHORT,size-1,0,
                     c[2],n,MPI_UNSIGNED_SHORT,1,0,MPI_COMM_WORLD,&status);
    if(p_rank == size-1)
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_SHORT,size-2,0,
                     c[2],n,MPI_UNSIGNED_SHORT,0,0,MPI_COMM_WORLD,&status);
    if((p_rank) && (p_rank < size-1))
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_SHORT,p_rank-1,0,
                     c[2],n,MPI_UNSIGNED_SHORT,p_rank+1,0,MPI_COMM_WORLD,&status);
    
    if(p_rank == 0)
        MPI_Sendrecv(c[0]+2*n,n,MPI_UNSIGNED_SHORT,1,1,
                     c[1],n,MPI_UNSIGNED_SHORT,size-1,1,MPI_COMM_WORLD,&status);
    if(p_rank == size-1)
        MPI_Sendrecv(c[0]+2*n,n,MPI_UNSIGNED_SHORT,0,1,
                     c[1],n,MPI_UNSIGNED_SHORT,size-2,1,MPI_COMM_WORLD,&status);
    if((p_rank) && (p_rank < size-1))
        MPI_Sendrecv(c[0]+2*n,n,MPI_UNSIGNED_SHORT,p_rank+1,1,
                     c[1],n,MPI_UNSIGNED_SHORT,p_rank-1,1,MPI_COMM_WORLD,&status);
}

#endif



/******************************************
 ******************************************
 *******  Free dynamic memory  ************
 ******************************************
 ******************************************/

void FreeMem () {
    
    delete [] T;
    delete [] S;
    delete [] nbdbh;
    delete [] layer;
    delete [] SPECIES_GERM;
    if(_SEEDTRADEOFF){
        delete [] PROB_S;
    }

    if(_NDD){
        delete [] PROB_S;
    }
    
    int h;
    for (h=0; h<(HEIGHT+1); h++) {
        delete [] LAI3D[h];
    }
    
    delete [] LAI3D;
    
    for (int i=0; i<3; i++) {
        delete [] Thurt[i];
    }
}


/***********************************
 ***********************************
 ***** RANDOM NUMBER GENERATOR *****
 ***********************************
 ***********************************/

/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

//#include<stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializing the array with a NONZERO seed */
void sgenrand2(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

/* generating reals */
/* unsigned long */ /* for integer generation */
double genrand2()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= N) { /* generate N words at one time */
        int kk;
        
        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand2(4357); /* a default initial seed is used   */
        
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];
        
        mti = 0;
    }
    
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    
    return ( (double)y / (unsigned long)0xffffffff ); /* reals */
    /* return y; */ /* for integer generation */
}


/* initializing the array with a NONZERO seed */
void sgenrand2i(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

unsigned long genrand2i()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= N) { /* generate N words at one time */
        int kk;
        
        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand2i(4357); /* a default initial seed is used   */
        
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];
        
        mti = 0;
    }
    
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    
    return y;
}
