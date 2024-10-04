#!/bin/bash
#SBATCH --job-name=Genotyping
#SBATCH --account=def-makarenk
#SBATCH --mem=64G
#SBATCH --array=1-13
#SBATCH --output=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.out
#SBATCH --error=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.err
#SBATCH --time=08:00:00

cd $SLURM_SUBMIT_DIR
module load nixpkgs/16.09 gatk/4.1.2.0

#Preset
#Set WD
WORK_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/HaplotypeCaller
OUT_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/Genotyping

#Create outdir
if [ ! -d $OUT_DIR ];then
  mkdir $OUT_DIR;
fi

#List all files in folder to file
if [ ! -f $WORK_DIR/files_genotype.txt ]; then
  for i in $WORK_DIR/*".g.vcf.gz"; do
    echo $i >> $WORK_DIR/files_genotype.txt;
  done
fi

#Set file list to object
inputFiles=$WORK_DIR/files_genotype.txt
#Isolate each input folder to job ID
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $inputFiles)
#output
out=$OUT_DIR/"$(basename -s ".g.vcf.gz" $file).vcf.gz"
ref=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa

gatk --java-options "-Xmx4g" GenotypeGVCFs -R $ref -V $file  -O $out

cd $SLURM_SUBMIT_DIR
echo $file $SLURM_TASK_ID $SLURM_ARRAY_TASK_ID >> _listTask_haplotype.txt
