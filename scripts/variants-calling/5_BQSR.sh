#!/bin/bash
#SBATCH --job-name=BQSR
#SBATCH --account=def-makarenk
#SBATCH --mem=64G
#SBATCH --array=1-13
#SBATCH --output=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.out
#SBATCH --error=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.err
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR
module load nixpkgs/16.09 gatk/4.1.2.0

#Get known variants + index -> interactive node
#mkdir db
#cd db
#module load samtools/1.10 picard/2.20.6 bcftools/1.10.2
#wget ftp://ftp.ensembl.org/pub/release-100/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz
#wget ftp://ftp.ensembl.org/pub/release-100/variation/vcf/homo_sapiens/homo_sapiens_structural_variations.vcf.gz
#gunzip *
#bcftools sort -o /home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_1000GENOMES-phase_3.vcf /home/nicdemon/projects/def-makarenk/nicdemon/db/1000GENOMES-phase_3.vcf
#bcftools sort -o /home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_homo_sapiens_structural_variations.vcf /home/nicdemon/projects/def-makarenk/nicdemon/db/homo_sapiens_structural_variations.vcf
#gatk IndexFeatureFile -F /home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_1000GENOMES-phase_3.vcf
#gatk IndexFeatureFile -F /home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_homo_sapiens_structural_variations.vcf

#Preset
#Set WD
WORK_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/splitNCigar
OUT_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/BQSR

#Create outdir
if [ ! -d $OUT_DIR ];then
  mkdir $OUT_DIR;
fi

#List all files in folder to file
if [ ! -f $WORK_DIR/files_BQSR.txt ]; then
  for i in $WORK_DIR/*".bam"; do
    echo $i >> $WORK_DIR/files_BQSR.txt;
  done
fi

#Set file list to object
inputFiles=$WORK_DIR/files_BQSR.txt
#Isolate each input folder to job ID
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $inputFiles)
out1=$OUT_DIR/"Before_$(basename -s ".bam" $file).table"
out2=$OUT_DIR/"$(basename -- $file)"
out3=$OUT_DIR/"After_$(basename -s ".bam" $file).table"
out4=$OUT_DIR/"Analysis_$(basename -s ".bam" $file).pdf"
ref=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa
sites1=/home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_1000GENOMES-phase_3.vcf
sites2=/home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_homo_sapiens_structural_variations.vcf

gatk BaseRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -I $file -R $ref --known-sites $sites1 --known-sites $sites2 -O $out1 --disable-sequence-dictionary-validation

gatk ApplyBQSR -R $ref -I $file --bqsr-recal-file $out1 -O $out2

gatk BaseRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -I $out2 -R $ref --known-sites $sites1 --known-sites $sites2 -O $out3 --disable-sequence-dictionary-validation

gatk AnalyzeCovariates --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -before $out1 -after $out3 -plots $out4

cd $SLURM_SUBMIT_DIR
echo $file $SLURM_TASK_ID $SLURM_ARRAY_TASK_ID >> _listTask_BQSR.txt
