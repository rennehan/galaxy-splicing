#!/bin/bash -l
#########################################################
#SBATCH -J ProtoVisualization
#SBATCH --mail-user=rennehan@uvic.ca
#SBATCH --mail-type=ALL 
#SBATCH --account=rrg-babul-ad
#SBATCH -o /home/rennehan/analysis/protocluster/logs/visual-%j.out
#########################################################
#SBATCH --time=2880
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=12000m
#########################################################
#: $SLURM_SUBMIT_DIR:=`pwd`
#echo $SLURM_SUBMIT_DIR
#cd $SLURM_SUBMIT_DIR

python /home/rennehan/analysis/protocluster/scripts/make_movie.py

#########################################################
