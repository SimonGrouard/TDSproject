#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=1:mem=60gb
#PBS -N 1node1core


cd /rds/general/project/hda_21-22/live/TDS/Group_6/scripts/multivariate_exposures
module load anaconda3/personal

Rscript StabilitySelection_imputed.R