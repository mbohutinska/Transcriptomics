#!/bin/bash
#PBS -l select=1:ncpus=4:mem=32gb:scratch_local=20gb
#PBS -l walltime=3:00:00
#PBS -j oe

# Load necessary modules
module load fastqc
module load trimmomatic

# Define paths
INPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/raw/gen_bas_flow_col_polym_A_arenosa_CP"
OUTPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/trimmed"
ADAPTERS="/storage/brno3-cerit/home/filip_kolar/brno3/600_genomes/adapters/TruSeq3-PE.fa"

SCRATCH_FASTQ="$SCRATCHDIR/fastq_files"
SCRATCH_OUTPUT="$SCRATCHDIR/output"

# Create directories in scratch space
mkdir -p $SCRATCH_FASTQ
mkdir -p $SCRATCH_OUTPUT

# Copy the specific input files to the scratch space
cp $INPUT_DIR/${SAMPLE_ID}_Lib*_L003_R1_001.fastq.gz $SCRATCH_FASTQ
cp $INPUT_DIR/${SAMPLE_ID}_Lib*_L004_R1_001.fastq.gz $SCRATCH_FASTQ
cp $INPUT_DIR/${SAMPLE_ID}_Lib*_L003_R2_001.fastq.gz $SCRATCH_FASTQ
cp $INPUT_DIR/${SAMPLE_ID}_Lib*_L004_R2_001.fastq.gz $SCRATCH_FASTQ

# Move to the scratch directory for processing
cd $SCRATCH_FASTQ

# Combine lane reads for each sample
cat ${SAMPLE_ID}_Lib*_L003_R1_001.fastq.gz ${SAMPLE_ID}_Lib*_L004_R1_001.fastq.gz > ${SAMPLE_ID}_R1_combined.fastq.gz
cat ${SAMPLE_ID}_Lib*_L003_R2_001.fastq.gz ${SAMPLE_ID}_Lib*_L004_R2_001.fastq.gz > ${SAMPLE_ID}_R2_combined.fastq.gz

# FastQC before trimming
fastqc ${SAMPLE_ID}_R1_combined.fastq.gz ${SAMPLE_ID}_R2_combined.fastq.gz -o $SCRATCH_OUTPUT

# Trimmomatic for adapter trimming
trimmomatic PE -threads 4 \
    ${SAMPLE_ID}_R1_combined.fastq.gz ${SAMPLE_ID}_R2_combined.fastq.gz \
    ${SAMPLE_ID}_R1_trimmed_paired.fastq.gz ${SAMPLE_ID}_R1_trimmed_unpaired.fastq.gz \
    ${SAMPLE_ID}_R2_trimmed_paired.fastq.gz ${SAMPLE_ID}_R2_trimmed_unpaired.fastq.gz \
    ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36

# FastQC after trimming
fastqc ${SAMPLE_ID}_R1_trimmed_paired.fastq.gz ${SAMPLE_ID}_R2_trimmed_paired.fastq.gz -o $SCRATCH_OUTPUT

# Move trimmed fastq files to scratch output directory
mv ${SAMPLE_ID}_R1_trimmed_paired.fastq.gz $SCRATCH_OUTPUT
mv ${SAMPLE_ID}_R2_trimmed_paired.fastq.gz $SCRATCH_OUTPUT

# Copy results back to the output directory
cp -r $SCRATCH_OUTPUT/* $OUTPUT_DIR || export CLEAN_SCRATCH=false

echo "Pipeline for sample ${SAMPLE_ID} completed."
