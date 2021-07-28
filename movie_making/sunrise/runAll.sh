#!/bin/bash -l
#########################################################
#SBATCH -J SUNRISE
#SBATCH --mail-user=rennehan@uvic.ca
#SBATCH --mail-type=ALL 
#SBATCH --account=rrg-babul-ad
#SBATCH -o /home/rennehan/work/sunrise/logs/sunrise-%j.out
#########################################################
#SBATCH --time=0-3:0
##SBATCH --ntasks=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
##SBATCH --mem-per-cpu=2200m
#SBATCH --mem=0
#SBATCH --array=0-43
#########################################################
#: $SLURM_SUBMIT_DIR:=`pwd`
#echo $SLURM_SUBMIT_DIR
#cd $SLURM_SUBMIT_DIR
#########################################################
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#snapDir=/home/rennehan/scratch/EddyDiffusion/MFM/TurbOff/TurbOff_QuickDiceTest/data
#snapDir=/home/rennehan/scratch/EddyDiffusion/MFM/TurbOff/TurbOff_Protocluster/data
snapDir=/home/rennehan/scratch/EddyDiffusion/MFM/TurbOff/TurbOff_Protocluster_HighRes/data

./pipeline.sh $snapDir $SLURM_ARRAY_TASK_ID 

#########################################################
