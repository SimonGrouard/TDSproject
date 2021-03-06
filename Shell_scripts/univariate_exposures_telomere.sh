#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=40gb
#PBS -N sevnodes1core
#PBS -J 1-24

cd /rds/general/project/hda_21-22/live/TDS/Group_6/scripts/univariate_exposures
module load anaconda3/personal

data_path=/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/
ichunk=$PBS_ARRAY_INDEX
nchunks=24

Rscript univariate_exposures_telomere_volcano.R $data_path $nchunks $ichunk
