#!/bin/sh
#SBATCH --account=eaton
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --workdir=/moto/eaton/users/slh2181/slurm
#SBATCH --job-name=bcl2fastq
#SBATCH --mail-type=END
#SBATCH --mail-user=sandra.hoffberg@columbia.edu


bcl2fastq --input-dir /moto/eaton/users/slh2181/Pmontana/illumina_runs/126985/Data/Intensities/BaseCalls/ \
-R /moto/eaton/users/slh2181/Pmontana/illumina_runs/126985/ \
--output-dir /moto/eaton/users/slh2181/Pmontana/illumina_runs --use-bases-mask y*,i*,i*,y* --no-lane-splitting --create-fastq-for-index-reads