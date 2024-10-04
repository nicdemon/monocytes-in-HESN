#!/bin/bash
#SBATCH --job-name=MarkDups
#SBATCH --account=def-makarenk
#SBATCH --mem=16G
#SBATCH --array=1-13
#SBATCH --output=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.out
#SBATCH --error=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.err
#SBATCH --time=01:00:00

cd $SLURM_SUBMIT_DIR
module load nixpkgs/16.09 sambamba/0.7.1

#Preset
#Set WD
WORK_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/merged
OUT_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/marked

#Create outdir
if [ ! -d $OUT_DIR ];then
  mkdir $OUT_DIR;
fi

#List all files in folder to file
if [ ! -f $WORK_DIR/files_aligned.txt ]; then
  for i in $WORK_DIR/*; do
    echo $i >> $WORK_DIR/files_merged.txt;
  done
fi

#Set file list to object
inputFiles=$WORK_DIR/files_merged.txt
#Isolate each input folder to job ID
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $inputFiles)
#test
#file=$(sed -n "1p" $inputFiles)
out=$OUT_DIR/"$(basename -- $file)"
metrics=$OUT_DIR/"metrics_$(basename -s .bam $file)"

sambamba markdup -t=48 $file $out

cd $SLURM_SUBMIT_DIR
echo $file $SLURM_TASK_ID $SLURM_ARRAY_TASK_ID >> _listTask_markdups.txt
