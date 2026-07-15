#!/bin/bash
#SBATCH --job-name=rgreat
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=rgreat_%j.out
#SBATCH --error=rgreat_%j.err

module purge

# Load R environment (same one you used before)
module load rstudio/v2024.12.1_563-gcc-13.2.0-r-4.5.1-python-3.11.6

echo "Starting rGREAT job..."
date

Rscript /scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/scripts/cluster_peaks_rGREAT.R

echo "Finished rGREAT job..."
date
