#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -N sevnodes1core
#PBS -J 1-24

cd /rds/general/project/hda_21-22/live/TDS/Group_6/scripts
module load anaconda3/personal

data_path=/rds/general/project/hda_21-22/live/TDS/Group_6/Results_univariate_exposures/
ichunk=$PBS_ARRAY_INDEX
nchunks=24

Rscript ExposureImputation.R $data_path $nchunks $ichunk