#!/bin/bash
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=100gb
#PBS -l walltime=3:00:00
#PBS -j oe

# Load necessary modules
module load py-htseq

SAMPLE_ID="5"

# Define paths
INPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/raw/gen_bas_flow_col_polym_A_arenosa_CP/trial/star"
OUTPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/htseq"
GTF_FILE="/storage/pruhonice1-ibot/home/holcovam/references/lyrataV2/Alyrata_384_v2.1.gene_exonsAGAT.gtf"

# Go to scratch directory
cd $SCRATCHDIR

# Copy input data to scratch
cp $INPUT_DIR/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam .
cp $GTF_FILE .

# Run HTSeq counting
htseq-count -f bam -r pos -s reverse -t exon -i gene_id  --nonunique none ${SAMPLE_ID}_Aligned.sortedByCoord.out.bam Alyrata_384_v2.1.gene_exonsAGAT.gtf > ${SAMPLE_ID}_htseq_counts.txt


#alternatives:
#htseq-count -f bam -r pos -s no -t exon -i gene_id --nonunique none --secondary-alignments ignore --additional-attr gene_name ${SAMPLE}_Aligned.sortedByCoord.out.bam /path/to/annotation.gtf > ${SAMPLE}_htseq_counts.txt

#htseq-count -m intersection-nonempty -s no -f bam -r pos -t exon -i Parent ./${i}*Aligned.sortedByCoord.out.bam /storage/brno12-cerit/home/susnatas/transcriptomics_AA372/X201SC22102893-Z01-F001/01.RawData/Arenosa_v2.gff > /storage/brno12-cerit/home/susnatas/transcriptomics_AA250/readcounts/${i}_read_counts.txt



# Copy results back to output directory
cp ${SAMPLE_ID}_htseq_counts.txt $OUTPUT_DIR || export CLEAN_SCRATCH=false

echo "HTSeq counting completed for sample ${SAMPLE_ID}."
