#!/bin/bash
#SBATCH -J wheat_OF2013
#SBATCH -o logs/%j__wheat_OF2013.out
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p Cnode2
#SBATCH --mem=256gb



set -e

module unload R
source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3


#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/code/data_prepare.R
Rscript /public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/code/univariate.R

echo "calculate done!"
