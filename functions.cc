#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <random>
#include <mpi.h>

#include "prng.h"

using namespace std;

static inline string int_to_string(const int a){
    ostringstream str;
    str << a;
    return str.str();
}

static inline string float_to_string(const double a){
    ostringstream str;
    str << a;
    return str.str();
}

void Initialize_rho(const int Nx, double *rho, int *positionList, const double *Parameters){
    for(int i=0;i<Parameters[P_nparticles];i++){
        int xi = randSite(Nx);
        positionList[i] = xi;
        rho[xi] += 1;
    }
}

void Read_file(const string fileName, double *rho, const int Nx){
    ifstream input(fileName.c_str(),ios::in|ios::binary);
    input.read((char*)rho,Nx*sizeof(double));
    input.close();
}    

void Save_rho(const string fileName, const double *rho, const int Nx){
    ofstream out(fileName.c_str(),ios::out|ios::binary);
    out.write((char*)rho,Nx*sizeof(double));
    out.close();
}

void Sample_rho(const int Nx, const int NSample, const int NSkip, double *fullRho, double *rhoSample){
    for(int j=0;j<=NSample;j+=NSkip){
        for(int i=0;i<Nx;i++){
            rhoSample[(j/NSkip)*Nx+i] = fullRho[j*Nx+i];
        }
    }
}

void Reset_rho(const int Nx, double *rho, double *fullRho){
    for(int i=0;i<Nx;i++){
        fullRho[i] = rho[i];
    }
}

void Sample_timedata(const int NSample, const int NSkip, double *PhysTime, double *timeSample){
    for(int j=0;j<=NSample;j+=NSkip){
        timeSample[j/NSkip] = PhysTime[j];
    }
}

static inline double Total_mass(const int Nx, const double *rho){
    double sum = 0.;
    for(int i=0;i<Nx;i++){
        sum += rho[i];
    }
    return sum;
}

void Normalize_rho(const int Nx, double *rho, const double *Parameters){
    double norm = Total_mass(Nx,rho);
    norm /= (Parameters[P_nparticles]);
    for(int i=0;i<Nx;i++)
        rho[i] /= norm;
}

// Heaviside function
static inline double kappa(const double rho, const double *Parameters)
{
    if(rho>=Parameters[P_nc]){
        return (pow((Parameters[P_nc]/rho),Parameters[P_eta]));
    }
    else{
        return 1.;
    }
}


// Harmonic Potential

void set_potential(const int Nx, double *Potential, const double *Parameters){
    for(int i=0;i<Nx;i++){
        Potential[i] = pow((i*Parameters[P_dx] - Parameters[P_L]),2);
	Potential[i] *= 0.5;
    }
}

static inline double jumpRate(int CurrentSite, int NextSite, const double *Potential, const double *Parameters){
    double vdiff = Potential[CurrentSite] - Potential[NextSite];
    return (2./(1. + exp(-vdiff*Parameters[P_beta])));
}
   
int Update_density(const int Nx, const int NSkip, int counter, double *rho, int *positionList, double *PhysTime, double *timeSample, double *fullRho, const double *Potential, const double *Parameters){

    int randParticle = randSite(Parameters[P_nparticles]);
    int CurrentSite = positionList[randParticle];
    //cout << "Picked up particle no " << endl; 
    int Direction = randDir();
    int NextSite = CurrentSite + Direction;

    // zero-flux BC

    if ((NextSite>=Nx) || (NextSite<0)){
        NextSite -= Direction;
    }

    double jumpProb = Parameters[P_dt]*jumpRate(CurrentSite,NextSite,Potential,Parameters);
    jumpProb *= kappa(rho[CurrentSite],Parameters);

    double coinToss = randNum();

    if (coinToss<jumpProb){
        //counter += 1;
        rho[CurrentSite] -= 1;
        rho[NextSite] += 1;
        positionList[randParticle] = NextSite;
    }
            
    double tStep = randNum();
    PhysTime[0] -= (0.5*Parameters[P_dt]*log(tStep));

    counter += 1;

    if(counter%NSkip==0){
	    for(int i=0;i<Nx;i++){
	        fullRho[(counter/NSkip)*Nx+i] = rho[i];
	    }
        timeSample[counter/NSkip] = PhysTime[0];
    }

    return counter;
}
