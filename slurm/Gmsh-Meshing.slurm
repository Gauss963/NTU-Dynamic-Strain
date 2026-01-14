#!/bin/bash
#SBATCH --account=MST113389
#SBATCH --partition=hm112
#SBATCH --job-name=Meshing
#SBATCH --output=/home/gauss112/outs/Meshing-%j.out
#SBATCH --error=/home/gauss112/outs/Meshing-%j.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=85:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b09501028@ntu.edu.tw

echo "[JOB] Meshing started at $(date) on $(hostname) with $SLURM_NTASKS cores."
source /etc/profile.d/modules.sh

module load tools/miniconda3
conda activate gmsh

cd /home/gauss112/NTU-Dynamic-Strain/Gmsh

python3 ./Block-Assembly.py