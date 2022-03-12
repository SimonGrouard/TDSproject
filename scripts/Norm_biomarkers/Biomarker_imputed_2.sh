#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=60gb
#PBS -N sevnodes1core
#PBS -J 1-24

cd /rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers
module load anaconda3/personal

data_path=/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final
ichunk=$PBS_ARRAY_INDEX
nchunks=24

Rscript Biomarker_imputation.R $data_path $nchunks $ichunk