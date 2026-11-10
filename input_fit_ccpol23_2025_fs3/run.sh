#!/bin/bash
#SBATCH --job-name=CCpol23-MPI
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
##SBATCH --exclusive
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=ccpol-biomix-%j.out
#SBATCH --error=ccpol-biomix-%j.err

# Submit from the simulation run directory containing the input files.
# Official guide: https://bioit.dbi.udel.edu/BIOMIX/SLURM-intro.html
# 96 is the largest CPU count listed for the default batch partition.
# Confirm the live configuration with: sinfo -p batch -N -o '%N %c %m %t'
set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:?Submit this script using sbatch}"

# Use the SAME Open MPI installation that compiled the executable.
MPI_PREFIX="${MPI_PREFIX:-${HOME}/local/openmpi-4.1.8}"
if [[ ! -x "${MPI_PREFIX}/bin/mpirun" ]]; then
    echo "ERROR: MPI launcher not found: ${MPI_PREFIX}/bin/mpirun" >&2
    echo "Complete the local Open MPI installation or set MPI_PREFIX." >&2
    exit 1
fi
export PATH="${MPI_PREFIX}/bin:${PATH}"
export LD_LIBRARY_PATH="${MPI_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Override at submission with --export=ALL,EXECUTABLE=/absolute/path/execute
#EXECUTABLE="${EXECUTABLE:-./execute}"
EXECUTABLE="$HOME/MD_ccpol23+/source_W_rf_check_SLR_ravg_rf_restart/execute"
if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "ERROR: Executable missing or not executable: ${EXECUTABLE}" >&2
    echo "Place execute in this run directory or set EXECUTABLE." >&2
    exit 1
fi

printf 'Job: %s | Node: %s | MPI ranks: %s\n' \
    "${SLURM_JOB_ID}" "${SLURM_JOB_NODELIST}" "${SLURM_NTASKS}"
printf 'Working directory: %s\nExecutable: %s\n' "$PWD" "${EXECUTABLE}"
"${MPI_PREFIX}/bin/mpirun" --version

# Single-node MPI launch; no Slurm PMI/PMIx plug-in is required for direct srun.
# CPU binding is left off because the guide does not specify CPU/core topology.
exec "${MPI_PREFIX}/bin/mpirun" \
    -np "${SLURM_NTASKS}" \
    --host "$(hostname):${SLURM_NTASKS}" \
    --map-by slot --bind-to none \
    "${EXECUTABLE}"

