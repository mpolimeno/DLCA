#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdbool>
#include <cmath>

#define MAXFNLEN      128                        // Maximum length of a filename
#define DIM           3                          // Dimensionality of the space
#define MAXNPARCL     100                        // Maximum number of particles in an aggreegate
#define TOTNPARCL     100                        // Total number of parcels in the system
#define BOXL          128.                       // This would be a fixed volume fraction
#define NAGG          1.                         // Number of aggregates allowed at final time
#define STOP          1400000 //set it large for now                   // Number of timesteps for the stopping criterion for N=100
#define R             1.                         // Radius of a single sphere
#define diam          (2.0*R)                    // Diameter of the sphere
#define D_T           0.5                        // D_T=0.5 for a single sphere
#define D_R           ((3./4.)*D_T)                    // Rotational diffusion coefficient for single solid sphere  (units of frequency)
#define Dt            0.01                       // Time-step
#define EPS           1.e-4                      // Small number to avoid issues with floating point arithmetic
#define Pi            (4.*atan(1.))              // Definition of Pi
#define C_u           1.1                       // this is the prefactor (2/9)*pi^2*g*\Delta{\rho}(1e-24/4e-21)
#define U_0           1.05                      // bias of a single sphere 2*pi*g*(1e-24*\Delta{\rho}/12e-21)

extern void combine_aggs(double *agg,int idx_1,int &nparcl1,int idx_2,int &nparcl2);
extern double compute_distsq(double *agg,int idx_1,int nparcl1,int idx_2,int nparcl2);
extern void init_aggs(double *aggs,int nparcl[],int nagg);
extern void print_aggs_data(FILE *fp,double *aggs,int nparcl[],int nagg);
extern void compute_com(double *agg,int idx,int nparcl,double cm[DIM]);
extern double find_max_rad(double maxrad[],int idx,double *agg,double com[DIM],int nparcl);
extern double find_gyr_rad(double Grad[],int idx,double *agg,double com[DIM],int nparcl);
extern void compute_xrel(double *agg, int idx, int nparcl,double com[DIM],double *xrel);
extern void compute_xlab(double *xrel,int idx,double com[DIM],int nparcl,double Q[][DIM],double *xlab);


#define ATTACH (diam+EPS)                       // Max distance for aggregagation to be allowed
#define SEPDISTSQ (ATTACH*ATTACH)               // Initial minimum distance squared allowed  b/w parcels
