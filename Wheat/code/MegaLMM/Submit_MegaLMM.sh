#!/bin/bash
#SBATCH -J wheat_MegaLMM
#SBATCH -o logs/%j__MegaLMM_1:20.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1gb
#SBATCH -t 7-00:00:00

set -e

source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3


for trialid in {1..20}; do

    for foldid in {6..20}; do
				
        sbatch /public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/MegaLMM/MegaLMM.sh ${trialid} ${foldid}
		
    done
	
done

#foldid=1

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

