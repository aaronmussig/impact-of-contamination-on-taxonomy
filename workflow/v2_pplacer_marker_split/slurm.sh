#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=330G
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_ace

PATH_REFPKG=/home/uqamussi/data/gtdb_r207_bac120.refpkg

PATH_MSA=/home/uqamussi/v2_pplacer_marker_split/msa/msa_$JOB_ID.fasta

PATH_JSON=/scratch/user/uqamussi/v2_pplacer_marker_split_$JOB_ID.json

PATH_STDOUT=/scratch/user/uqamussi/v2_pplacer_marker_split_msa_$JOB_ID.stdout
PATH_STDERR=/scratch/user/uqamussi/v2_pplacer_marker_split_msa_$JOB_ID.stderr

export PATH=$PATH:/home/uqamussi/.conda/envs/gunc-chimeras/bin

srun --job-name=v2_pplacer_marker_split_$JOB_ID pplacer -m wag -j 1 -c /home/uqamussi/data/gtdb_r207_bac120.refpkg -o $PATH_JSON $PATH_MSA > $PATH_STDOUT 2> $PATH_STDERR

mv $PATH_STDOUT /home/uqamussi/v2_pplacer_marker_split/results/pplacer_$JOB_ID.stdout
mv $PATH_STDERR /home/uqamussi/v2_pplacer_marker_split/results/pplacer_$JOB_ID.stderr
mv $PATH_JSON /home/uqamussi/v2_pplacer_marker_split/results/pplacer_$JOB_ID.json


# run as
#export JOB_ID=1_0 && sbatch --job-name pp1_0 main.sh
#export JOB_ID=1_1 && sbatch --job-name pp1_1 main.sh
#export JOB_ID=1_2 && sbatch --job-name pp1_2 main.sh
#export JOB_ID=1_3 && sbatch --job-name pp1_3 main.sh
#export JOB_ID=1_4 && sbatch --job-name pp1_4 main.sh
#export JOB_ID=1_5 && sbatch --job-name pp1_5 main.sh
#export JOB_ID=1_6 && sbatch --job-name pp1_6 main.sh
#export JOB_ID=1_7 && sbatch --job-name pp1_7 main.sh
#export JOB_ID=2_0 && sbatch --job-name pp2_0 main.sh
#export JOB_ID=2_1 && sbatch --job-name pp2_1 main.sh
#export JOB_ID=2_2 && sbatch --job-name pp2_2 main.sh
#export JOB_ID=2_3 && sbatch --job-name pp2_3 main.sh
#export JOB_ID=2_4 && sbatch --job-name pp2_4 main.sh
#export JOB_ID=2_5 && sbatch --job-name pp2_5 main.sh
#export JOB_ID=2_6 && sbatch --job-name pp2_6 main.sh
#export JOB_ID=2_7 && sbatch --job-name pp2_7 main.sh
#export JOB_ID=3_0 && sbatch --job-name pp3_0 main.sh
#export JOB_ID=3_1 && sbatch --job-name pp3_1 main.sh
#export JOB_ID=3_2 && sbatch --job-name pp3_2 main.sh
#export JOB_ID=3_3 && sbatch --job-name pp3_3 main.sh
#export JOB_ID=3_4 && sbatch --job-name pp3_4 main.sh
#export JOB_ID=3_5 && sbatch --job-name pp3_5 main.sh
#export JOB_ID=3_6 && sbatch --job-name pp3_6 main.sh
#export JOB_ID=3_7 && sbatch --job-name pp3_7 main.sh
#export JOB_ID=4_0 && sbatch --job-name pp4_0 main.sh
#export JOB_ID=4_1 && sbatch --job-name pp4_1 main.sh
#export JOB_ID=4_2 && sbatch --job-name pp4_2 main.sh
#export JOB_ID=4_3 && sbatch --job-name pp4_3 main.sh
#export JOB_ID=4_4 && sbatch --job-name pp4_4 main.sh
#export JOB_ID=4_5 && sbatch --job-name pp4_5 main.sh
#export JOB_ID=4_6 && sbatch --job-name pp4_6 main.sh
#export JOB_ID=4_7 && sbatch --job-name pp4_7 main.sh
#export JOB_ID=5_0 && sbatch --job-name pp5_0 main.sh
#export JOB_ID=5_1 && sbatch --job-name pp5_1 main.sh
#export JOB_ID=5_2 && sbatch --job-name pp5_2 main.sh
#export JOB_ID=5_3 && sbatch --job-name pp5_3 main.sh
#export JOB_ID=5_4 && sbatch --job-name pp5_4 main.sh
#export JOB_ID=5_5 && sbatch --job-name pp5_5 main.sh
#export JOB_ID=5_6 && sbatch --job-name pp5_6 main.sh
#export JOB_ID=5_7 && sbatch --job-name pp5_7 main.sh
#export JOB_ID=10_0 && sbatch --job-name pp10_0 main.sh
#export JOB_ID=10_1 && sbatch --job-name pp10_1 main.sh
#export JOB_ID=10_2 && sbatch --job-name pp10_2 main.sh
#export JOB_ID=10_3 && sbatch --job-name pp10_3 main.sh
#export JOB_ID=10_4 && sbatch --job-name pp10_4 main.sh
#export JOB_ID=10_5 && sbatch --job-name pp10_5 main.sh
#export JOB_ID=10_6 && sbatch --job-name pp10_6 main.sh
#export JOB_ID=10_7 && sbatch --job-name pp10_7 main.sh
#export JOB_ID=15_0 && sbatch --job-name pp15_0 main.sh
#export JOB_ID=15_1 && sbatch --job-name pp15_1 main.sh
#export JOB_ID=15_2 && sbatch --job-name pp15_2 main.sh
#export JOB_ID=15_3 && sbatch --job-name pp15_3 main.sh
#export JOB_ID=15_4 && sbatch --job-name pp15_4 main.sh
#export JOB_ID=15_5 && sbatch --job-name pp15_5 main.sh
#export JOB_ID=15_6 && sbatch --job-name pp15_6 main.sh
#export JOB_ID=15_7 && sbatch --job-name pp15_7 main.sh
#export JOB_ID=20_0 && sbatch --job-name pp20_0 main.sh
#export JOB_ID=20_1 && sbatch --job-name pp20_1 main.sh
#export JOB_ID=20_2 && sbatch --job-name pp20_2 main.sh
#export JOB_ID=20_3 && sbatch --job-name pp20_3 main.sh
#export JOB_ID=20_4 && sbatch --job-name pp20_4 main.sh
#export JOB_ID=20_5 && sbatch --job-name pp20_5 main.sh
#export JOB_ID=20_6 && sbatch --job-name pp20_6 main.sh
#export JOB_ID=20_7 && sbatch --job-name pp20_7 main.sh
#export JOB_ID=30_0 && sbatch --job-name pp30_0 main.sh
#export JOB_ID=30_1 && sbatch --job-name pp30_1 main.sh
#export JOB_ID=30_2 && sbatch --job-name pp30_2 main.sh
#export JOB_ID=30_3 && sbatch --job-name pp30_3 main.sh
#export JOB_ID=30_4 && sbatch --job-name pp30_4 main.sh
#export JOB_ID=30_5 && sbatch --job-name pp30_5 main.sh
#export JOB_ID=30_6 && sbatch --job-name pp30_6 main.sh
#export JOB_ID=30_7 && sbatch --job-name pp30_7 main.sh
#export JOB_ID=40_0 && sbatch --job-name pp40_0 main.sh
#export JOB_ID=40_1 && sbatch --job-name pp40_1 main.sh
#export JOB_ID=40_2 && sbatch --job-name pp40_2 main.sh
#export JOB_ID=40_3 && sbatch --job-name pp40_3 main.sh
#export JOB_ID=40_4 && sbatch --job-name pp40_4 main.sh
#export JOB_ID=40_5 && sbatch --job-name pp40_5 main.sh
#export JOB_ID=40_6 && sbatch --job-name pp40_6 main.sh
#export JOB_ID=40_7 && sbatch --job-name pp40_7 main.sh


# sbatch --job-name pp1 main.sh







