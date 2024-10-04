#!/bin/bash
#SBATCH --job-name=Merge
#SBATCH --account=def-makarenk
#SBATCH --mem=64G
#SBATCH --array=1-13
#SBATCH --output=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.out
#SBATCH --error=/home/nicdemon/projects/def-makarenk/nicdemon/scripts/%x_%a.err
#SBATCH --time=01:00:00

cd $SLURM_SUBMIT_DIR
module load nixpkgs/16.09 intel/2018.3 samtools/1.10

#Preset
#Set WD
WORK_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/aligned
OUT_DIR=/home/nicdemon/projects/def-makarenk/nicdemon/results/merged

#Create outdir
if [ ! -d $OUT_DIR ];then
  mkdir $OUT_DIR;
fi

#Set files list to object
#Isolate each input files to job ID
if [ ${SLURM_ARRAY_TASK_ID} == 1 ]; then
  if [ ! -f $WORK_DIR/files_D307TM24.txt ]; then
    for i in $WORK_DIR/"D307TM24"*".bam"; do
      echo $i >> $WORK_DIR/files_D307TM24.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D307TM24.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D307TM24.txt
  out=$OUT_DIR/D307TM24.bam
  rg=$WORK_DIR/rg_D307TM24.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 2 ]; then
  if [ ! -f $WORK_DIR/files_D336TM18.txt ]; then
    for i in $WORK_DIR/*"D336TM18"*".bam"; do
      echo $i >> $WORK_DIR/files_D336TM18.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D336TM18.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D336TM18.txt
  out=$OUT_DIR/D336TM18.bam
  rg=$WORK_DIR/rg_D336TM18.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 3 ]; then
  if [ ! -f $WORK_DIR/files_D377TM12.txt ]; then
    for i in $WORK_DIR/*"D377TM12"*".bam"; do
      echo $i >> $WORK_DIR/files_D377TM12.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D377TM12.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D377TM12.txt
  out=$OUT_DIR/D377TM12.bam
  rg=$WORK_DIR/rg_D377TM12.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 4 ]; then
  if [ ! -f $WORK_DIR/files_D619TM24.txt ]; then
    for i in $WORK_DIR/*"D619TM24"*".bam"; do
      echo $i >> $WORK_DIR/files_D619TM24.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D619TM24.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D619TM24.txt
  out=$OUT_DIR/D619TM24.bam
  rg=$WORK_DIR/rg_D619TM24.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 5 ]; then
  if [ ! -f $WORK_DIR/files_D704TM30.txt ]; then
    for i in $WORK_DIR/*"D704TM30"*".bam"; do
      echo $i >> $WORK_DIR/files_D704TM30.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D704TM30.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D704TM30.txt
  out=$OUT_DIR/D704TM30.bam
  rg=$WORK_DIR/rg_D704TM30.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 6 ]; then
  if [ ! -f $WORK_DIR/files_D714TM27.txt ]; then
    for i in $WORK_DIR/*"D714TM27"*".bam"; do
      echo $i >> $WORK_DIR/files_D714TM27.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D714TM27.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D714TM27.txt
  out=$OUT_DIR/D714TM27.bam
  rg=$WORK_DIR/rg_D714TM27.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 7 ]; then
  if [ ! -f $WORK_DIR/files_D715TM24.txt ]; then
    for i in $WORK_DIR/*"D715TM24"*".bam"; do
      echo $i >> $WORK_DIR/files_D715TM24.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D715TM24.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D715TM24.txt
  out=$OUT_DIR/D715TM24.bam
  rg=$WORK_DIR/rg_D715TM24.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 8 ]; then
  if [ ! -f $WORK_DIR/files_D722TM45.txt ]; then
    for i in $WORK_DIR/*"D722TM45"*".bam"; do
      echo $i >> $WORK_DIR/files_D722TM45.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D722TM45.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D722TM45.txt
  out=$OUT_DIR/D722TM45.bam
  rg=$WORK_DIR/rg_D722TM45.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 9 ]; then
  if [ ! -f $WORK_DIR/files_D749TM24.txt ]; then
    for i in $WORK_DIR/*"D749TM24"*".bam"; do
      echo $i >> $WORK_DIR/files_D749TM24.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D749TM24.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D749TM24.txt
  out=$OUT_DIR/D749TM24.bam
  rg=$WORK_DIR/rg_D749TM24.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 10 ]; then
  if [ ! -f $WORK_DIR/files_D752TM12.txt ]; then
    for i in $WORK_DIR/*"D752TM12"*".bam"; do
      echo $i >> $WORK_DIR/files_D752TM12.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_D752TM12.txt
    done
  fi
  inputFiles=$WORK_DIR/files_D752TM12.txt
  out=$OUT_DIR/D752TM12.bam
  rg=$WORK_DIR/rg_D752TM12.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 11 ]; then
  if [ ! -f $WORK_DIR/files_HLA_003.txt ]; then
    for i in $WORK_DIR/*"HLA_003"*".bam"; do
      echo $i >> $WORK_DIR/files_HLA_003.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_HLA_003.txt
    done
  fi
  inputFiles=$WORK_DIR/files_HLA_003.txt
  out=$OUT_DIR/HLA_003.bam
  rg=$WORK_DIR/rg_HLA_003.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 12 ]; then
  if [ ! -f $WORK_DIR/files_HLA_006.txt ]; then
    for i in $WORK_DIR/*"HLA_006"*".bam"; do
      echo $i >> $WORK_DIR/files_HLA_006.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_HLA_006.txt
    done
  fi
  inputFiles=$WORK_DIR/files_HLA_006.txt
  out=$OUT_DIR/HLA_006.bam
  rg=$WORK_DIR/rg_HLA_006.txt
elif [ ${SLURM_ARRAY_TASK_ID} == 13 ]; then
  if [ ! -f $WORK_DIR/files_HLA_009.txt ]; then
    for i in $WORK_DIR/*"HLA_009"*".bam"; do
      echo $i >> $WORK_DIR/files_HLA_009.txt;
      j="$(basename -s .bam $i)"
      perl -e 'print "\@RG\tID:'$j'\tSM:hs\tLB:'$j'\tPL:Illumina\n"' >> $WORK_DIR/rg_HLA_009.txt
    done
  fi
  inputFiles=$WORK_DIR/files_HLA_009.txt
  out=$OUT_DIR/HLA_009.bam
  rg=$WORK_DIR/rg_HLA_009.txt
fi

samtools merge -b $inputFiles -rh $rg -p --threads 48 $out

cd $SLURM_SUBMIT_DIR
echo $file $SLURM_TASK_ID $SLURM_ARRAY_TASK_ID >> _listTask_merge.txt
