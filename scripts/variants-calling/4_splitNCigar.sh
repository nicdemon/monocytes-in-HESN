#!/bin/bash
#SBATCH --job-name=SplitNCigar
#SBATCH --account=def-makarenk
#SBATCH --mem=16G
#SBATCH --array=1-13
#SBATCH --output=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.out
#SBATCH --error=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.err
#SBATCH --time=08:00:00

cd $SLURM_SUBMIT_DIR
module load nixpkgs/16.09 gatk/4.1.2.0

#Preset
#Set WD
WORK_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/marked
OUT_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/splitNCigar

#Create outdir
if [ ! -d $OUT_DIR ];then
  mkdir $OUT_DIR;
fi

#List all files in folder to file
if [ ! -f $WORK_DIR/files_splitNcigar.txt ]; then
  for i in $WORK_DIR/*".bam"; do
    echo $i >> $WORK_DIR/files_splitNcigar.txt;
  done
fi

#Set file list to object
inputFiles=$WORK_DIR/files_splitNcigar.txt
#Isolate each input folder to job ID
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $inputFiles)
out=$OUT_DIR/"$(basename -- $file)"
ref=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa

####GATK SPLITNCIGAR VERIFY OPTIONS + OUTPUT FOR NAME

gatk SplitNCigarReads -R $ref -I $file -O $out

cd $SLURM_SUBMIT_DIR
echo $file $SLURM_TASK_ID $SLURM_ARRAY_TASK_ID >> _listTask_splitNcigar.txt
