#!/bin/bash
#SBATCH --job-name=annotation
#SBATCH --account=def-makarenk
#SBATCH --mem=64G
#SBATCH --array=1-13
#SBATCH --output=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.out
#SBATCH --error=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.err
#SBATCH --time=08:00:00

cd $SLURM_SUBMIT_DIR
module load nixpkgs/16.09 intel/2018.3 gatk/4.1.2.0 bcftools/1.10.2

#Convert gff -> bedfile + index
#module load mugqic/bedops/v2.4.35
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.chr_patch_hapl_scaff.annotation.gff3.gz
#(grep ^"#" /home/nicdemon/projects/def-makarenk/nicdemon/db/gencode.v35.chr_patch_hapl_scaff.annotation.gff3; grep -v ^"#" /home/nicdemon/projects/def-makarenk/nicdemon/db/gencode.v35.chr_patch_hapl_scaff.annotation.gff3 | sort -k1,1 -k4,4n) | bgzip > /home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_gencode.v35.chr_patch_hapl_scaff.annotation.gff.gz
#gff2bed --keep-header < /home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_gencode.v35.chr_patch_hapl_scaff.annotation.gff > /home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_gencode.v35.chr_patch_hapl_scaff.annotation.bed
#sed 's/ \+ /\t/g' /home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_gencode.v35.chr_patch_hapl_scaff.annotation.bed > /home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_tab_gencode.v35.chr_patch_hapl_scaff.annotation.bed
#bgzip sorted_tab_gencode.v35.chr_patch_hapl_scaff.annotation.bed
#tabix -p bed project/def-makarenk/nicdemon/db/sorted_tab_gencode.v35.chr_patch_hapl_scaff.annotation.bed.gz

#Preset
#Set WD
WORK_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/Genotyping
OUT_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/annotation

#Create outdir
if [ ! -d $OUT_DIR ];then
  mkdir $OUT_DIR;
fi

#List all files in folder to file
if [ ! -f $WORK_DIR/files_annotation.txt ]; then
  for i in $WORK_DIR/*".vcf.gz"; do
    echo $i >> $WORK_DIR/files_annotation.txt;
  done
fi

#Set file list to object
inputFiles=$WORK_DIR/files_annotation.txt
#Isolate each input folder to job ID
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $inputFiles)
#output
out=$OUT_DIR/"$(basename -- $file)"
ref=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa
#set column headers
columns=CHROM,FROM,-,ID,SCORE,STRAND,SOURCE,TYPE,PHASE,ATTRIBUTES
#Set annotations bed to object
annot=/home/nicdemon/projects/def-makarenk/nicdemon/db/sorted_tab_gencode.v35.chr_patch_hapl_scaff.annotation.bed.gz
header=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/header.txt

#Run BCFtools to annotate variants
bcftools annotate -a $annot -c $columns -h $header -o $out -O z --threads 48 $file

#Index vcf after creation
gatk IndexFeatureFile -F $out

cd $SLURM_SUBMIT_DIR
echo $file $SLURM_TASK_ID $SLURM_ARRAY_TASK_ID >> _listTask_annotation.txt
