#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=1:mem=40gb
#PBS -N 1node1core


cd /rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/
module load anaconda3/personal


Rscript lasso_exposures.R