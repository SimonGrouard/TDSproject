#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=50gb
#PBS -N 1node1core

cd /rds/general/project/hda_21-22/live/TDS/Group_6/scripts/pls_analysis
module load anaconda3/personal

data_path=/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/

Rscript univariate_exposures_telomere_volcano.R $data_path