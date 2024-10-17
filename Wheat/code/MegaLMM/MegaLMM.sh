#!/bin/bash
#SBATCH -J wheat_MegaLMM
#SBATCH -o logs/MegaLMM/%j__MegaLMM.out
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=100gb
#SBATCH -p Cnode_all

set -e

source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3


echo "Trial is :${1} Fold is :${2}"

Rscript /public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/MegaLMM/megalmm_gblup.R ${1} ${2} > logs/MegaLMM/MegaLMM_gblup_${1}_${2}_rep.Rout

#foldid=1
#
#for trialid in {1..20};do
#
#	#Rscript /storage2/xueyh/TempWork/h_matrix_wheat_mix/MegaLMM/prep_data_trial.R ${trialid}
#
#	for foldid in {1..5};do
#
#		echo "Trial is :${trialid} Fold is :${foldid}"
#
#		Rscript /storage2/xueyh/TempWork/h_matrix_wheat_mix/MegaLMM/megalmm_gblup.R ${trialid} ${foldid} > logs/MegaLMM/MegaLMM_gblup_${trialid}_${foldid}_rep.Rout
#
#	done
#done

