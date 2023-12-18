#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=08:20:00
#SBATCH --partition=general
#SBATCH --account=a_ace

PATH_STDOUT=/scratch/user/uqamussi/v2_fasttree_marker_split_msa_$JOB_ID.stdout
PATH_STDERR=/scratch/user/uqamussi/v2_fasttree_marker_split_msa_$JOB_ID.stderr

export PATH=$PATH:/home/uqamussi/.conda/envs/gunc-chimeras/bin
export OMP_NUM_THREADS=24

srun --job-name=v2_fasttree_marker_split_msa_$JOB_ID FastTreeMP -wag /home/uqamussi/v2_fasttree_marker_split/msa/msa_$JOB_ID.faa > $PATH_STDOUT 2> $PATH_STDERR

mv $PATH_STDOUT /home/uqamussi/v2_fasttree_marker_split/results/fasttree_$JOB_ID.tree
mv $PATH_STDERR /home/uqamussi/v2_fasttree_marker_split/results/fasttree_$JOB_ID.stderr

# run as
#export JOB_ID=1 && sbatch --job-name ft1 main.sh
#export JOB_ID=2 && sbatch --job-name ft2 main.sh
#export JOB_ID=3 && sbatch --job-name ft3 main.sh
#export JOB_ID=4 && sbatch --job-name ft4 main.sh
#export JOB_ID=5 && sbatch --job-name ft5 main.sh
#export JOB_ID=10 && sbatch --job-name ft10 main.sh
#export JOB_ID=15 && sbatch --job-name ft15 main.sh
#export JOB_ID=20 && sbatch --job-name ft20 main.sh
#export JOB_ID=30 && sbatch --job-name ft30 main.sh
#export JOB_ID=40 && sbatch --job-name ft40 main.sh

