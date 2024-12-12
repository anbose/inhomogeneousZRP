#!/bin/bash

# Use bash as shell
#$ -S /bin/bash

# Preserve environment variables
#$ -v $PATH

# Execute from current working directory
#$ -cwd

# Merge standard output and standard error into one file
#$ -j yes

# parallel environment
#$ -pe mvapich2-sam 64

#echo "Got $NSLOTS slots."


#$ -v LD_LIBRARY_PATH=/usr/lmp/intel/lib/intel64:/usr/lmp/intel/mkl/lib/intel64:/usr/lib/x86_64-redhat-linux6E/lib64:/usr/lmp/gmp-6.1.2-gcc-8.2.0/lib:/usr/lmp/gcc-8.2.0/lib64:/usr/lmp/mpc-1.1.0-gcc-8.2.0/lib:/usr/lmp/mpfr-4.0.1-gcc-8.2.0/lib

# Standard name of the job (if none is given on the command line): e.g. "bacterialCells"
#$ -N N9e3_Nc153

# Some diagnostic messages for the output
echo Started: `date`
echo on `hostname`
echo Beta = 20
echo ------------

Np=$"9000"

outputDirectory=$"/scratch03.local/abose/lattice_model/dataNewParams/N"$Np$"/"

Nc=$"153"
#TOld=$"10.508"

mkdir $outputDirectory
mkdir $outputDirectory$"Beta5"$"_Nc"$Nc$"/"

make clean

make PDE

mpirun -n $NSLOTS ./mainRKMC.exe $Np $Nc

make clean

echo ------------

echo Nc = $Nc

echo Stopped: `date`
