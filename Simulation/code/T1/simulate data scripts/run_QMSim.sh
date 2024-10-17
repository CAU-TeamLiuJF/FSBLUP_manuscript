#!/bin/bash
#SBATCH -J simulate_data_c2021
#SBATCH -p Cnode_all
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4
#SBATCH --time=30:20:00 

module load compiler/intel-compiler/2021.3.0
module load mpi/openmpi/intel/4.0.3
module unload R;source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3
export OMP_NUM_THREADS=4
#source ~/.bash_profile
rm -r r_datasim
./QMSim datasim.prm
