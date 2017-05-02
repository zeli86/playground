#!/bin/bash
#PBS -j oe
#PBS -l walltime=11:59:00
#PBS -l nodes=17:ppn=24
#PBS -l feature=mpp
#PBS -A hbp00035
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out

cd 0010
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0015
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0020
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0025
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0030
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0035
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0040
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0045
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0050
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0055
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0060
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0065
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0070
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0075
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0080
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0085
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
cd ../0090
cp ../params.prm . 
aprun -n 24 rt_prop_mpi_cs &
wait
 
