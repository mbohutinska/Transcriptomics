#!/bin/bash
#PBS -l select=1:ncpus=16:mem=68gb:scratch_local=50gb
#PBS -l walltime=6:50:00
#PBS -j oe

# Load necessary modules
module load star

# Define paths
INPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/trimmed"
OUTPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/star"
REFERENCE_DIR="/storage/pruhonice1-ibot/home/holcovam/references/lyrataV2/starRef" 
#SAMPLE_ID="5"

# Go to scratch directory
cd $SCRATCHDIR

# Copy input data to scratch
cp $INPUT_DIR/${SAMPLE_ID}_R1_trimmed_paired.fastq.gz .
cp $INPUT_DIR/${SAMPLE_ID}_R2_trimmed_paired.fastq.gz .
cp -r $REFERENCE_DIR . # Copy genome index to scratch

# Run STAR alignment
STAR --runThreadN 16 \
     --genomeDir starRef \
     --readFilesIn ${SAMPLE_ID}_R1_trimmed_paired.fastq.gz ${SAMPLE_ID}_R2_trimmed_paired.fastq.gz \
     --outFileNamePrefix ${SAMPLE_ID}_ \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts \
     --outSAMstrandField intronMotif \
     --readFilesCommand zcat \
     --twopassMode Basic

# Clean up intermediate files
rm ${SAMPLE_ID}_R1_trimmed_paired.fastq.gz ${SAMPLE_ID}_R2_trimmed_paired.fastq.gz
rm -r starRef # Remove the local genome index copy from scratch
#mkdir ${SAMPLE_ID}
#mv -r * ${SAMPLE_ID}

# Copy results back to output directory
cp -r * $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Pipeline completed."
