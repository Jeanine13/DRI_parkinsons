#!/bin/bash
#SBATCH --job-name=rGREAT_4view
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/logs/rGREAT_4view_%j.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/logs/rGREAT_4view_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu
module purge

module load rstudio/v2024.12.1_563-gcc-13.2.0-r-4.5.1-python-3.11.6
mkdir -p /scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/logs

export R_LIBS_USER=/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5

Rscript /scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/scripts/cluster_peaks_rGREAT_4view.R
