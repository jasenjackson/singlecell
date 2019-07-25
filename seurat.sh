#!/bin/bash
#SBATCH -A b1042                                        # Allocation
#SBATCH -p genomics                                     # Queue
#SBATCH -t 1:00:00                                     # Walltime/duration of the job
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4800
#SBATCH --mail-user=jasen.jackson@northwestern.edu       # Designate email address for job communications
#SBATCH --mail-type=FAIL                                # Events options are job BEGIN, END, NONE, FAIL, REQUEUE

module load R

Rscript seurat_integrate.R
