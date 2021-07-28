#!/bin/bash -l
#########################################################
#SBATCH -J MCRX
#SBATCH --mail-user=rennehan@uvic.ca
#SBATCH --mail-type=ALL 
#SBATCH --account=rrg-babul-ad
#SBATCH -o /home/rennehan/work/sunrise/logs/mcrx-%j.out
#########################################################
#SBATCH --time=1440
##SBATCH --ntasks=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=2200m
#########################################################
#: $SLURM_SUBMIT_DIR:=`pwd`
#echo $SLURM_SUBMIT_DIR
#cd $SLURM_SUBMIT_DIR
#########################################################

cd /home/rennehan/software/sunrize-1.0/bin

./mcrx /home/rennehan/configs/sunrise/mcrx.params

#########################################################
