### Petals:
Read count data: normalizedCounts_filtered_240RC_15nonzero_7Oct.csv

Sample data: SampleMetadata.txt

### Leaves:
Read count data: Guillaume_96samples_readCounts.txt

Sample data:  Alignment_96samples_compil_fromGuillaume.xlsx, sample_map...

### R script to directly analyse raw read counts: tkpAnalysis23Feb26.R

### On how to get the raw read counts start by reading the tkpNotes file (attached or copied bellow)

### 1. Trimming and QC

example input
5_Lib_S105_L003_R1_001.fastq.gz
5_Lib_S105_L003_R2_001.fastq.gz
5_Lib_S105_L004_R1_001.fastq.gz
5_Lib_S105_L004_R2_001.fastq.gz
18_Lib_S106_L003_R1_001.fastq.gz
18_Lib_S106_L003_R2_001.fastq.gz
18_Lib_S106_L004_R1_001.fastq.gz
18_Lib_S106_L004_R2_001.fastq.gz

Submit in parallel
INPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/raw/gen_bas_flow_col_polym_A_arenosa_CP"
Extract unique sample IDs from input files
SAMPLE_IDS=($(ls ${INPUT_DIR}/*_R1_001.fastq.gz | sed -E 's|.*/([0-9]+)_.*|\1|' | sort | uniq))

Loop over the array and submit the job for each sample
for SAMPLE in "${SAMPLE_IDS[@]}"; do
    qsub -v SAMPLE_ID=$SAMPLE -N 1_QC_trimming_$SAMPLE 1_QC_trimming.sh
done


### 2. Mapping with STAR

2.0 make genome index: 
qsub auto/pruhonice1-ibot/nfs4/home/holcovam/tkpPetals/2.0_starIndexing.sh


2.1 run mapping in parallel:

example input
18_R2_trimmed_paired_fastqc.html
18_R2_trimmed_paired_fastqc.zip
18_R2_trimmed_paired.fastq.gz

Submit in parallel
INPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/trimmed/"
Extract unique sample IDs from input files
SAMPLE_IDS=($(ls ${INPUT_DIR}/*_R1_trimmed_paired.fastq.gz | sed -E 's|.*/([0-9]+)_.*|\1|' | sort | uniq))

SAMPLE_IDS=("100" "54" "18")
Loop over the array and submit the job for each sample
for SAMPLE in "${SAMPLE_IDS[@]}"; do
    qsub -v SAMPLE_ID=$SAMPLE -N 2_star_$SAMPLE 2.1_star_parallel.sh
done

### 3. get read counts

#### 3.1 directly from STAR
scp holcovam@tilia.metacentrum.cz:/auto/pruhonice1-ibot/home/holcovam/tkpPetals/star/*_ReadsPerGene.out.tab .

R pipeline locally

#### 3.2 htseq-count https://htseq.readthedocs.io/en/release_0.11.1/count.html
My ideal options are already in the star counting pipeline (excude multiple mappers, ignore feature overlap as this not a common case in our reference), no need to refine. my prefered htseq command: 
htseq-count -f bam -r pos -s reverse -t exon -i gene_id  --nonunique none ${SAMPLE_ID}_Aligned.sortedByCoord.out.bam Alyrata_384_v2.1.gene_exonsAGAT.gtf > ${SAMPLE_ID}_htseq_counts.txt
Gave the very same results as star :D


Define input directory and extract sample IDs
INPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/star"
SAMPLE_IDS=($(ls ${INPUT_DIR}/*_Aligned.sortedByCoord.out.bam | sed -E 's|.*/([0-9]+)_.*|\1|' | sort | uniq))

Loop over the sample IDs and submit the jobs
for SAMPLE in "${SAMPLE_IDS[@]}"; do
    qsub -v SAMPLE_ID=$SAMPLE -N 3_htseq_$SAMPLE 3_htseq_parallel.sh
done


