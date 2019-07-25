#!/bin/bash
#SBATCH -A b1042                                        # Allocation
#SBATCH -p genomics                                     # Queue
#SBATCH -t 24:00:00                                     # Walltime/duration of the job
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=4800
#SBATCH --mail-user=jasenjackson@northwestern.edu       # Designate email address for job communications
#SBATCH --mail-type=FAIL                                # Events options are job BEGIN, END, NONE, FAIL, REQUEUE

# unload any modules that carried over from your command line session
module load cellranger/3.0.1

cellranger count --id=MCF7_AA_t0_1 \
                 --transcriptome=refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=MCF7_AA_t0
