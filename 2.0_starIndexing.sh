#!/bin/bash
#PBS -l select=1:ncpus=16:mem=70gb:scratch_local=100gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -N 2.0_starIndex

module load star

# Define paths
INPUT="/auto/pruhonice1-ibot/nfs4/home/holcovam/references/lyrataV2/"
OUTPUT="/auto/pruhonice1-ibot/nfs4/home/holcovam/references/lyrataV2/starRef"

# Go to scratch directory
cd $SCRATCHDIR

# Copy necessary input files to scratch
cp $INPUT/alygenomes.fasta .
cp $INPUT/Alyrata_384_v2.1.gene_exonsAGAT.gtf .

# Create a directory for the genome index
mkdir -p genomeDir

# Run STAR to generate the genome index
STAR --runMode genomeGenerate \
     --genomeDir genomeDir \
     --genomeFastaFiles ./alygenomes.fasta \
     --sjdbGTFfile ./Alyrata_384_v2.1.gene_exonsAGAT.gtf \
     --sjdbOverhang 150 \
     --genomeSAindexNbases 12

# Copy the genome index files back to the output directory
cp -r genomeDir/* $OUTPUT || export CLEAN_SCRATCH=false

echo "Done"

