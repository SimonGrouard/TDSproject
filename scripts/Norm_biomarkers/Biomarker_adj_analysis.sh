#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=1:mem=30gb
#PBS -N 1node1core


cd /rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers
module load anaconda3/personal

Rscript Biomarker_adj_analysis.R 