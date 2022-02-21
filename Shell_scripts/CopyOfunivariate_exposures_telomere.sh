#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=40gb
#PBS -N 1node1core

cd /rds/general/project/hda_21-22/live/TDS/Group_6/scripts
module load anaconda3/personal

data_path=/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/
start_chunk=580
end_chunk=613

Rscript CopyOfunivariate_exposures_telomere.R $data_path $start_chunk $end_chunk
