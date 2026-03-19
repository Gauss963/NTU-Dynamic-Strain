#!/bin/bash
#SBATCH --account=MST113389
#SBATCH --partition=hm112
#SBATCH --job-name=Akantu-LSW-2D-2.0%dt-20ms
#SBATCH --output=/home/gauss112/outs/Akantu-%j.out
#SBATCH --error=/home/gauss112/outs/Akantu-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=95:59:59
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b09501028@ntu.edu.tw

echo "[JOB] Job started at $(date) on $(hostname) with $SLURM_NTASKS cores."
source /etc/profile.d/modules.sh

module load gcc/11.2.0
module load openmpi/5.0.2

DATE_STR=$(date +"%Y-%m-%d-%H-%M-%S")-${SLURM_JOB_ID}
DEST_DIR="/work1/gauss112/Akantu-Run/$DATE_STR"

mkdir -p "$DEST_DIR"
cp -r /home/gauss112/NTU-Dynamic-Strain/ "$DEST_DIR/"
echo "[JOB] Job directory $DEST_DIR"
cd "$DEST_DIR/NTU-Dynamic-Strain/Akantu/Velocity-weakening/build/"
mpirun -np $SLURM_NTASKS ./Velocity-weakening-LSW