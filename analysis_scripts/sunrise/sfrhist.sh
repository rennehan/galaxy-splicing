#!/bin/bash -l
#########################################################
#SBATCH -J SfrHist
#SBATCH --mail-user=rennehan@uvic.ca
#SBATCH --mail-type=ALL 
#SBATCH --account=rrg-babul-ad
#SBATCH -o /home/rennehan/work/sunrise/logs/sfrhist-%j.out
#########################################################
#SBATCH --time=180
##SBATCH --ntasks=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=2200m
#########################################################
#: $SLURM_SUBMIT_DIR:=`pwd`
#echo $SLURM_SUBMIT_DIR
#cd $SLURM_SUBMIT_DIR
#########################################################

cd /home/rennehan/software/sunrize-1.0/bin

./sfrhist /home/rennehan/configs/sunrise/sfrhist.params

#########################################################
