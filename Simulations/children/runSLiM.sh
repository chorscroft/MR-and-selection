#!/bin/bash
#PBS -l walltime=00:10:00
cd $PBS_O_WORKDIR

module load apps/slim-3.4

slim $mySim
