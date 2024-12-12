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
#include "pcg-cpp-0.98/include/pcg_random.hpp"

//Random number generator (pcg64)
//#include "pcg-c-0.94/include/pcg_variants.h"
//#include "pcg-c-0.94/extras/entropy.h"

#include "Parameters.h"
#include "functions.cc"

using namespace std;

/*

Simulation of Lattice model with site-dependent hopping rate using a Rejection Kinetic Monte Carlo (RKMC) algorithm

Author : Aritra Bose, date : 19/12/2023

*/

int *positionList;
double *rho, *fullRho;
//double *localavg, *globalavg;
double *PhysTime, *timeSample;
double *Potential;

int main(int argc, char* argv[])
{
    time_t begin = time(NULL);
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    
    L = 1;
    Nx = 100;
    
    Nt = 1e9;
    dt = 0.5;
    NSample = 1e3;
    NSkip = 1e5;

    k = 1.;
    eta = 1.;

    NParticles = atoi(argv[1]);
    Nc = atoi(argv[2]);

    Beta = 5;

    Set_Parameters();

    rho = new double[Nx];
    //localavg = new double[Nx];
    //globalavg = new double[Nx];
    
    fullRho = new double[Nx*(NSample+1)];
    PhysTime = new double[1];
    timeSample = new double[NSample+1];

    Potential = new double[Nx];
    positionList = new int[NParticles];

    string folderName = "/scratch03.local/abose/lattice_model/dataNewParams/N"+int_to_string(NParticles)+"/Beta"+float_to_string(Beta);
    folderName += "_Nc"+int_to_string(Nc)+"/Run"+int_to_string(myid+1)+"/";  
    //cout << folderName << endl;
    int status = mkdir(folderName.c_str(),S_IRWXU);
    (void)status;

    set_potential(Nx,Potential,Parameters);
    //Save_rho(folderName+"Potential.bin",Potential,Nx);

    Initialize_rho(Nx,rho,positionList,Parameters);
    /*string readFileName = "/scratch03.local/abose/lattice_model/dataNc/otherN/r0.017/run2/N"+int_to_string(NParticles)+"/Beta"+float_to_string(Beta);
    readFileName += "_Nc"+int_to_string(Nc)+"/Run"+int_to_string(myid+1)+"/Density_field10.bin";
    Read_file(readFileName,fullRho,Nx*(NSample+1));
    for(int i=0;i<Nx;i++){
        rho[i] = fullRho[NSample*Nx+i];
    }*/

    MPI_Barrier(MPI_COMM_WORLD);

    Save_rho(folderName+"Density_field_t0.bin",rho,Nx);
    //if(myid==0){
    //MaxRho[0] = rho[Nx/2];
    //}

    int counter = 0;
    int fileNum = 0;

    PhysTime[0] = 0.;
    Reset_rho(Nx,rho,fullRho);
    //for(int i=0;i<Nx;i++){
    //	fullRho[counter*Nx+i] = rho[i];
    //}
    timeSample[0] = 0.;

    for(long int t=1;t<=Nt;t++){
        counter = Update_density(Nx,NSkip,counter,rho,positionList,PhysTime,timeSample,fullRho,Potential,Parameters);

        if(counter/NSkip==NSample){
            fileNum += 1;

	        //cout << "Saving: counter stands at : " << counter*fileNum << endl;

            //Sample_rho(Nx,NSample,NSkip,fullRho,rhoSample);
            //Sample_timedata(NSample,NSkip,PhysTime,timeSample);

            //Save_rho(folderName+"timedata"+int_to_string(fileNum)+".bin",timeSample,NSample+1);
            Save_rho(folderName+"Density_field"+int_to_string(fileNum)+".bin",fullRho,Nx*(NSample+1));
	        //Save_rho(folderName+"Density_field_t"+int_to_string(counter*fileNum)+".bin",rho,Nx);

            timeSample[0] = timeSample[NSample];

            Reset_rho(Nx,rho,fullRho);

            counter = 0;
        }

    }

/*
    fileNum += 1;

    double *restTime = new double[counter+1];
    double *restRho = new double[counter+1];
    
    for(int i=0;i<=counter;i++){
        restTime[i] = PhysTime[i];
        restRho[i] = MaxRho[i];
    }
    
    //Save_rho(folderName+"timedata"+int_to_string(fileNum)+".bin",restTime,counter+1);
    Save_rho(folderName+"rhomax"+int_to_string(fileNum)+".bin",restRho,counter+1);
*/

    Save_rho(folderName+"Density_field_t_final.bin",rho,Nx);

    //delete[] restTime;
    //delete[] restRho;

    //cout << "current counter stands at : " << counter << endl;

    /*Normalize_rho(Nx,localavg,Parameters);
    Save_rho(folderName+"Avg_density.bin",localavg,Nx);*/

    //cout << "mean density (rescaled) at ground state : " << localavg[Nx/2]/Nc << endl;

    MPI_Barrier(MPI_COMM_WORLD);

    /*MPI_Reduce(localavg,globalavg,Nx,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    if(myid==0){
        Normalize_rho(Nx,globalavg,Parameters);
	    cout << "mean density (rescaled) at ground state : " << globalavg[Nx/2]/Nc << endl;
        Save_rho(folderName+"meanrho.bin",globalavg,Nx);
    }*/

    time_t end = time(NULL);

    if(myid==0){
        double Ratio = static_cast<double>(Nc)/static_cast<double>(NParticles);
        cout << "run time : " << (end-begin) << " seconds" << endl;
        cout << "N : " << Total_mass(Nx,rho) << ", Nc = " << Nc << ", Nc/N = " << Ratio << ", beta = " << Beta << endl;
        cout << "System size is : " << Parameters[P_L] << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    delete[] rho;
    //delete[] localavg;
    //delete[] globalavg;
    //delete[] MaxRho;
    delete[] fullRho;
    delete[] PhysTime;
    delete[] timeSample;

    delete[] Potential;
    delete[] positionList;
    
    MPI_Finalize();

    return 0;
}
