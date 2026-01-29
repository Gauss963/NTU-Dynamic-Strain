#!/bin/bash
#SBATCH --job-name="Akantu-LSW-50mm"
#SBATCH --partition=cpu-2g
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24-4:30
#SBATCH --chdir=.
#SBATCH --output=/home/r14521220/outs/Akantu-%j.out
#SBATCH --error=/home/r14521220/outs/Akantu-%j.err

echo "[JOB] Job started at $(date) on $(hostname) with $SLURM_NTASKS cores."

DATE_STR=$(date +"%Y-%m-%d-%H-%M-%S")
DEST_DIR="/home/r14521220/work/Akantu-Run/$DATE_STR"

mkdir -p "$DEST_DIR"
cp -r /home/r14521220/NTU-Dynamic-Strain/ "$DEST_DIR/"
echo "[JOB] Job directory $DEST_DIR"
cd "$DEST_DIR/NTU-Dynamic-Strain/Akantu/Velocity-weakening/build/"
mpirun -np 1 ./Velocity-weakening-LSW