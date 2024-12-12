#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>

enum{P_L,P_dt,P_dx,P_k,P_lambda,P_eta,P_nc,P_beta,P_nparticles,Nparameters};

double L,dx,dt;
long int Nt;
int Nx,NSample,NSkip,NParticles,Nc;
double k,Lx,eta;
double Beta;
double meanrho;

double eps = 1.e-100;

double Parameters[Nparameters];

int nproc, myid, errorcode;
ptrdiff_t local_n0,local_n0_start;

static inline void Set_Parameters()
{
    dx = 2.*L/(double)Nx;

    Parameters[P_L] = L;
    Parameters[P_dx] = dx;
    Parameters[P_dt] = dt;
    Parameters[P_k] = k;
    Parameters[P_lambda] = Lx;
    Parameters[P_eta] = eta;
    Parameters[P_nc] = Nc;
    Parameters[P_beta] = Beta;
    Parameters[P_nparticles] = NParticles;

}

