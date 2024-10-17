#!/bin/bash
#SBATCH -J pig_Mix_fsmat
#SBATCH -o logs/%j__pig_Mix_fsmat.out
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p Cnode_all
#SBATCH --mem=100gb
#SBATCH --exclude=cnode1053,cnode1007


set -e

module unload R
source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3


#Trails=(xxx 'MD2013' 'OB2013' 'H2013' 'SD2013' 'OF2013' 'MD2014' 'OB2014' 'H2014' 'SD2014' 'OF2014' 'MD2015' 'OB2015' 'H2015' 'SD2015' 'OF2015' 'MD2016' 'OB2016' 'H2016' 'SD2016' 'OF2016')
#
#for i in {1..20}; do
#    echo "Trials is :"${Trails[$i]}
#    Rscript data_prepare.R $i
#    wait
#    Rscript univariate.R $i > logs/${Trails[$i]}_miss_0_30.Rout
#    wait
#done



#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/code/CalTime.R
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/code/CalTime_Mix.R 2
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/code/CalTime_Mix_perct_acc.R 4
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/code/CalTime_BL.R 1
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/code/sce2_c2021.R
Rscript /public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/code/ssGOBLUP.R

#args=/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/code/ssGOBLUP.R
#parallel -j 3 Rscript "$args {}" ::: $(seq -w 1 3)


echo "calculate done!"
