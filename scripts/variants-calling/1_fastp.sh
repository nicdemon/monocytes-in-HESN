#!/bin/bash
#SBATCH --job-name=FastP
#SBATCH --account=def-makarenk
#SBATCH --mem=64G
#SBATCH --array=1-13
#SBATCH --output=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.out
#SBATCH --error=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.err
#SBATCH --time=00:20:00

cd $SLURM_SUBMIT_DIR
module load nixpkgs/16.09  intel/2018.3 fastp/0.20.0

#Preset
#Set WD
DATA_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/data
WORK_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/trimmed

#Create outdir
if [ ! -d $WORK_DIR ];then
  mkdir $WORK_DIR;
fi

#List all folders to file
if [ ! -f $WORK_DIR/files_folder_fastq.txt ]; then
  for i in $DATA_DIR/*; do
    echo $i >> $WORK_DIR/files_folder_fastq.txt;
  done
fi

#Set folder list to object
inputFolders=$WORK_DIR/files_folder_fastq.txt
#Isolate each input folder to job ID
folder=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $inputFolders)

for i in $folder/*; do
  out=$WORK_DIR/"$(basename -- $i)"
  fastp -i $i -o $out
done

cd $SLURM_SUBMIT_DIR
echo $file $SLURM_TASK_ID $SLURM_ARRAY_TASK_ID >> _listTask_fastP.txt
