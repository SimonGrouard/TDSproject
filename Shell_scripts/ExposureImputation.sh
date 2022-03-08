#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=1:mem=50gb
#PBS -N sevnodes1core
#PBS -J 1-24

cd /rds/general/project/hda_21-22/live/TDS/Group_6/scripts
module load anaconda3/personal

data_path=/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded
ichunk=$PBS_ARRAY_INDEX
nchunks=24

Rscript ExposureImputation.R $data_path $nchunks $ichunk