#!/bin/bash
#SBATCH -J simulate_h_fsmat
#SBATCH -o logs/%j__simulate_h_fsmat.out
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p Cnode_all
#SBATCH --mem=250gb



set -e
module unload R
source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3


Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/code/univariate.R 
Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/code/ssGBLUP_scenario2.R
Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/code/ssGBLUP_scenario3.R
Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/code/bayes.R 
Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/code/FSBLUP.R 


echo "calculate done!"
