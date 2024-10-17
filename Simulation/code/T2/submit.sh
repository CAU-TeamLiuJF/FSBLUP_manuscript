#!/bin/bash
#SBATCH -J simulate_l_fsmat
#SBATCH -o logs/%j__simulate_l_fsmat.out
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p Cnode_all
#SBATCH --mem=100gb



set -e

source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3

#wd=/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code
#for i in {1..3}; do
#    for j in {1..20}; do
##        sbatch "$wd/Extract_genetic_T.sh" $i $j     # adjust T mat
#        echo "Scenario is :$i, rep is :$j"
##        Rscript "$wd/data_prepare.R" $i $j
##        wait
##
#        sbatch "$wd/Mix.sh" $i $j
##
#    done
#done

#echo "univariate.R"
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code/univariate.R
#echo "ssGBLUP_scenario2.R" 
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code/ssGBLUP_scenario2.R
#echo "ssGBLUP_scenario3.R"
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code/ssGBLUP_scenario3.R
#echo "simulate_low_h2.R"
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code/sce1_c2021.R
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code/sce2_c2021.R
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code/sce3_c2021.R
#Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code/bayes.R
Rscript /public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code/ssGOBLUP.R 

#args=/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/code/ssGOBLUP.R 
#parallel -j 3 Rscript "$args {}" ::: $(seq -w 1 3)

echo "calculate done!"
