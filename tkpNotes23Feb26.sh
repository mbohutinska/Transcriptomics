################### RNA seq ###################

##### The clean procedure

### 1. Trimming and QC

#example input
5_Lib_S105_L003_R1_001.fastq.gz
5_Lib_S105_L003_R2_001.fastq.gz
5_Lib_S105_L004_R1_001.fastq.gz
5_Lib_S105_L004_R2_001.fastq.gz
18_Lib_S106_L003_R1_001.fastq.gz
18_Lib_S106_L003_R2_001.fastq.gz
18_Lib_S106_L004_R1_001.fastq.gz
18_Lib_S106_L004_R2_001.fastq.gz

### Submit in parallel
INPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/raw/gen_bas_flow_col_polym_A_arenosa_CP"
# Extract unique sample IDs from input files
SAMPLE_IDS=($(ls ${INPUT_DIR}/*_R1_001.fastq.gz | sed -E 's|.*/([0-9]+)_.*|\1|' | sort | uniq))

# Loop over the array and submit the job for each sample
for SAMPLE in "${SAMPLE_IDS[@]}"; do
    qsub -v SAMPLE_ID=$SAMPLE -N 1_QC_trimming_$SAMPLE 1_QC_trimming.sh
done






### 2. Mapping with STAR

## 2.0 make genome index: 
qsub auto/pruhonice1-ibot/nfs4/home/holcovam/tkpPetals/2.0_starIndexing.sh


## 2.1 run mapping in parallel:

#example input
18_R2_trimmed_paired_fastqc.html
18_R2_trimmed_paired_fastqc.zip
18_R2_trimmed_paired.fastq.gz

### Submit in parallel
INPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/trimmed/"
# Extract unique sample IDs from input files
SAMPLE_IDS=($(ls ${INPUT_DIR}/*_R1_trimmed_paired.fastq.gz | sed -E 's|.*/([0-9]+)_.*|\1|' | sort | uniq))

SAMPLE_IDS=("100" "54" "18")
# Loop over the array and submit the job for each sample
for SAMPLE in "${SAMPLE_IDS[@]}"; do
    qsub -v SAMPLE_ID=$SAMPLE -N 2_star_$SAMPLE 2.1_star_parallel.sh
done



### 3. get read counts

## 3.1 directly from STAR
scp holcovam@tilia.metacentrum.cz:/auto/pruhonice1-ibot/home/holcovam/tkpPetals/star/*_ReadsPerGene.out.tab .

R pipeline locally

## 3.2 htseq-count https://htseq.readthedocs.io/en/release_0.11.1/count.html
My ideal options are already in the star counting pipeline (excude multiple mappers, ignore feature overlap as this not a common case in our reference), no need to refine. my prefered htseq command: 
htseq-count -f bam -r pos -s reverse -t exon -i gene_id  --nonunique none ${SAMPLE_ID}_Aligned.sortedByCoord.out.bam Alyrata_384_v2.1.gene_exonsAGAT.gtf > ${SAMPLE_ID}_htseq_counts.txt
Gave the very same results as star :D







## 3.2 using HTSeq

# Define input directory and extract sample IDs
INPUT_DIR="/storage/pruhonice1-ibot/home/holcovam/tkpPetals/star"
SAMPLE_IDS=($(ls ${INPUT_DIR}/*_Aligned.sortedByCoord.out.bam | sed -E 's|.*/([0-9]+)_.*|\1|' | sort | uniq))

# Loop over the sample IDs and submit the jobs
for SAMPLE in "${SAMPLE_IDS[@]}"; do
    qsub -v SAMPLE_ID=$SAMPLE -N 3_htseq_$SAMPLE 3_htseq_parallel.sh
done







##### Notes

#Advantages of STAR:
Speed: STAR is generally much faster than BWA, especially for large RNA-seq datasets.
Splice Junction Detection: STAR is designed for RNA-seq and excels in detecting splice junctions, making it more accurate for transcriptomic mapping.
Multi-mapping: It handles multi-mapping reads better, which is crucial in highly repetitive transcriptomic regions.
Disadvantages of STAR:
High Memory Demand: STAR’s performance comes at the cost of requiring large amounts of RAM, sometimes limiting its use on resource-constrained systems.
Complex Output: The output of STAR can be more complex to interpret compared to BWA, requiring more effort in downstream analysis.
#Advantages of BWA:
Lower Memory Usage: BWA is more efficient in memory usage, making it suitable for smaller computing environments.
Accuracy in Simple Alignment: BWA is highly accurate for general DNA and RNA alignments, especially in simpler transcriptomic datasets where splice junction detection is less critical.
Wide Adoption: BWA is widely used and well-supported for a range of applications, not just RNA-seq.
Disadvantages of BWA:
Slower: BWA is slower than STAR, especially for large RNA-seq datasets.
Splice Junction Handling: BWA is not optimized for detecting splice junctions, making it less accurate for transcriptomic reads with complex splicing patterns.


############# Mapping with STAR
module load cufflinks
grep -v "miRNA_gene" LyV2.gff > LyV2.miRNA_geneRemoved.gff
gffread Alyrata_384_v2.1.gene_exons.gff3 -T -o Alyrata_384_v2.1.gene_exons.gtf
cut -f3 LyV2_23Sep.gtf | sort | uniq

gffread Alyrata_384_v2.1.gene_exons.gff3 -T -o Alyrata_384_v2.1.gene_exons.gtf --keep-genes --force-exons --include-utr --add-ids

gffread Alyrata_384_v2.1.gene_exons.gff3 -F -T -o Alyrata_384_v2.1.gene_exons.gtf 
cut -f3 Alyrata_384_v2.1.gene_exons.gtf  | sort | uniq
# Always problem with only CDS and exons!


# This finally sorted it out!
module load conda-modules   # load (default) module
conda activate agat
agat_convert_sp_gff2gtf.pl --gff Alyrata_384_v2.1.gene_exons.gff3 -o Alyrata_384_v2.1.gene_exonsAGAT.gtf --relax
cut -f3 Alyrata_384_v2.1.gene_exonsAGAT.gtf  | sort | uniq

# CDS
# exon
# five_prime_UTR
# gene
# mRNA
# three_prime_UTR


### make index: better use scratch in the future, takes a long time!
module load star
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ./alygenomes.fasta --sjdbGTFfile ./Alyrata_384_v2.1.gene_exonsAGAT.gtf --sjdbOverhang 150 --genomeSAindexNbases 12


### run STAR
# Mahnaz
INPUT="put you fastq files' directory"
INPUT1"put the directory of your index file" ### This index you produced it once and use it for your mapping of all samples. 
OUTPUT="the output directory"

cd $SCRATCHDIR
cp -r $INPUT/*.fastq.gz .
cp -r $INPUT1/genomeDir .

STAR --genomeDir ./genomeDir --runThreadN 16 --readFilesCommand zcat --readFilesIn *_R1_paired.fq.gz *_R2_paired.fq.gz --quantMode GeneCounts

mkdir -p $OUTPUT
rm -rf *_R[1,2]_paired.fq.gz genomeDir
cp -r * $OUTPUT || export CLEAN_SCRATCH=false

# chatgpt
STAR --runThreadN 8 \
     --genomeDir /path/to/genomeDir \
     --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
     --outFileNamePrefix sample_ \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts \
     --outSAMstrandField intronMotif \
     --readFilesCommand zcat \
     --twopassMode Basic

--quantMode TranscriptomeSAM GeneCounts: Outputs transcript-level BAM files and gene count files for downstream analysis (DE and DS).
--twopassMode Basic: Improves accuracy by aligning twice, which helps in detecting novel splice junctions.
This workflow ensures that your mapping process is ready for both differential expression and differential splicing analysis.


#CARLOS
STAR --twopassMode Basic --sjdbFileChrStartEnd $parsed_list --chimSegmentMin ${read_length}/3 --alignSJDBoverhangMin 3 --alignIntronMin 70 --alignIntronMax 562000 --alignMatesGapMax 562000 --limitSjdbInsertNsj 2000000 --outSAMtype BAM SortedByCoordinate --runThreadN $nthreads --genomeDir $genome_outdir --outFileNamePrefix ${mapout_dir}/ --readFilesPrefix ${data_dir}/ --readFilesCommand gunzip -c --readFilesIn ${sample_r1} ${sample_r2}


### STAR output
The `readPerGene.out.tab` file generated by STAR contains crucial information about how the reads from your RNA-seq data map to the reference genome. Here's how to interpret the columns in this output:

### Column Breakdown

1. **First Column**: Gene Identifier
   - This usually represents the gene ID, such as `AL1G10010.v2.1`, which indicates the name of the gene based on the annotation.

2. **Second Column**: Number of Reads Assigned to the Gene
   - This column indicates the total number of reads that were uniquely mapped to this specific gene. 

3. **Third Column**: Number of Reads Assigned to the Gene (for Multiple Mappings)
   - This column shows the number of reads that are considered multimapping, meaning they can map to multiple locations. These are generally reads that don't provide definitive evidence for the expression of any specific gene.

4. **Fourth Column**: Number of Reads Assigned to the Gene (including Unassigned)
   - This column represents the total reads assigned to the gene, combining both unique and multimapping reads. In some contexts, it might also include reads that are not confidently assigned to any features.

### Interpretation of Specific Lines

- **N_unmapped**: Total number of reads that did not map to any feature in the genome (6,313,116 reads).
  
- **N_multimapping**: Total number of reads that mapped to multiple locations (2,931,999 reads).
  
- **N_noFeature**: Reads that did not map to any gene (651,529 reads). The high number here suggests that a significant proportion of reads could not be associated with any known features in your reference genome.

- **N_ambiguous**: Reads that are ambiguous, meaning they could map to multiple genes but are not uniquely assignable (781,743 reads).

### Gene-Specific Data

- For example, the gene `AL1G10010.v2.1` has:
  - **81926** unique reads,
  - **41** multimapping reads,
  - **81912** total reads including those that are either uniquely assigned or ambiguously mapped.

### Final Thoughts

This output is essential for evaluating the quality of your mapping and understanding the distribution of reads across genes. If there are a significant number of unmapped or multimapping reads, it might indicate issues with the reference genome or the quality of the RNA-seq data. 

Additionally, when conducting downstream analyses like differential expression, you would typically focus on the unique reads to gauge the expression levels of genes, while considering multimapping reads may provide insights into gene family expressions or similar scenarios.






###### DE featureCounts, HTSeq or STAR directly??
#Yes, it’s possible to use the output of STAR directly for differential expression (DE) analysis in edgeR, as STAR can generate gene-level counts if `--quantMode GeneCounts` is specified in the command. These counts are stored in a file named `ReadsPerGene.out.tab`, which contains raw counts for each gene across different alignment categories (such as unstranded, stranded forward, and stranded reverse). These raw counts can be input into edgeR for downstream DE analysis.

#However, the advantage of adding a **featureCounts** step lies in flexibility and precision:

#1. **Control over counting options**: featureCounts offers greater control over how multi-mapping reads, exon-overlapping reads, and strand-specific reads are handled. You can define custom parameters to improve the accuracy of gene quantification.

#2. **Handling of gene annotation updates**: featureCounts allows you to specify more detailed or updated annotation files in GTF/GFF format, which can be useful when the provided annotation isn't ideal.

#3. **Additional statistics**: featureCounts provides detailed summary statistics, including how many reads are uniquely mapped, multi-mapping, and unassigned, which can be useful for quality control.

#In short, while STAR provides convenient and quick count generation, featureCounts offers more detailed control and precision. If your project demands fine-tuning or quality control, featureCounts is usually preferred for counting.

### Mahnaz:
#I haven't worked with featurecounts and I am not familiar with it, I always use STAR counting or Htseq. But for the differential analysis, edgeR and Deseq2 are two of the best packages which Deseq2 outperformed (I have this from some literature and forums).
Deseq2 is good for different designs (factorial designs, having batch effect and ...). I always use Deseq2, but edgeR is good as well.
I am not with rMATS and HDST either.
here is the manual of Htseq which is very straightforward.
And here is the link explaining how it counts the aligned reads (which I was trying to tell you today).
For running htseq, I ran it interactively and I did not save the command I just wrote what parameters I add besides the default parameters. I asked Susnata to send you the command.
htseq-count -m intersection-nonempty -s no -f bam -r pos -t exon -i Parent ./${i}*Aligned.sortedByCoord.out.bam /storage/brno12-cerit/home/susnatas/transcriptomics_AA372/X201SC22102893-Z01-F001/01.RawData/Arenosa_v2.gff > /storage/brno12-cerit/home/susnatas/transcriptomics_AA250/readcounts/${i}_read_counts.txt

#My trial:
htseq-count -f bam -r pos -s reverse -t exon -i gene_id  --nonunique none ${SAMPLE_ID}_Aligned.sortedByCoord.out.bam Alyrata_384_v2.1.gene_exonsAGAT.gtf > ${SAMPLE_ID}_htseq_counts.txt


#For your setup, an ideal HTSeq command should address the following considerations: paired-end reads, sorting of BAM files by coordinate, exclusion of multimapped reads, and compatibility with the annotation GTF file. Here’s a suggested command:

htseq-count -f bam -r pos -s no -t exon -i gene_id --nonunique none --secondary-alignments ignore --additional-attr gene_name ${SAMPLE}_Aligned.sortedByCoord.out.bam /path/to/annotation.gtf > ${SAMPLE}_htseq_counts.txt

Explanation:
-f bam: Specifies that the input files are in BAM format.
-r pos: Indicates that the BAM file is sorted by position (sortedByCoord from STAR).
-s no: No strand-specific information (adjust if your experiment requires strand-specific counting).
-t exon: Specifies that only exonic features will be counted.
-i gene_id: The GTF attribute to group reads by (gene-level analysis, assuming gene_id is the identifier in your GTF file).
--nonunique none: Excludes multimapping reads, as your goal is to avoid counting them.
--secondary-alignments ignore: Ignores secondary alignments, which can be produced for multimapping reads.
--additional-attr gene_name: This adds gene names alongside gene IDs in the output, which might be helpful in downstream analysis.
This command should be compatible with your STAR results and genome characteristics. HTSeq will count reads mapped to genes based on exonic regions, while excluding multimapped reads, making the output cleaner and more reliable for DESeq2. It’s also optimized for plant genomes, where some degree of repetition might affect read mapping.


### RESULTS: htseq my optional = STAR counting
			starReverse	htseq
AL1G10010.v2.1	81926	41	81912	81912
AL1G10020.v2.1	141	6	138	138
AL1G10030.v2.1	5	0	5	5
AL1G10040.v2.1	146	42	136	136pridat
AL1G10050.v2.1	1424	0	1456	1456
AL1G10060.v2.1	306	0	306	306
AL1G10070.v2.1	744	4	743	743
AL1G10080.v2.1	1028	4	1028	1028
AL1G10090.v2.1	4703	3	4707	4707
AL1G10102.v2.1	223	2	223	223
AL1G10104.v2.1	4803	2	4804	4804
AL1G10110.v2.1	4127	5	4125	4125
AL1G10120.v2.1	15	3	15	15
AL1G10130.v2.1	1014	0	1014	1014
AL1G10140.v2.1	2324	13	2311	2311
AL1G10150.v2.1	1	1	0	0








################## RANDOM NOTES FOR BETY
######## DS rMATS
To use rMATS for differential splicing (DS) analysis, here’s a general workflow and some key considerations:

### Steps for DS Analysis with rMATS:

1. **Prepare the Input Data**:
   - **BAM Files**: You’ll need aligned BAM files from your STAR mapping. These should be sorted by coordinates (as you’ve already done with the `--outSAMtype BAM SortedByCoordinate` option).
   - **GTF File**: The GTF file used during the STAR alignment will also be needed for rMATS.

2. **Run rMATS**:
   - **Input Data**: Organize your BAM files into groups based on your experimental design (e.g., condition A vs. condition B).
   - **rMATS Command**:
     ```bash
     rMATS-turbo --b1 <file_with_conditionA_BAMs> --b2 <file_with_conditionB_BAMs> \
                 --gtf <your_GTF_file> --od <output_directory> \
                 --nthread <number_of_threads> --readLength <read_length>
     ```
     - `--b1`: File listing BAM paths for condition A.
     - `--b2`: File listing BAM paths for condition B.
     - `--gtf`: The GTF file you used during STAR mapping.
     - `--readLength`: The length of your reads (important for accurate splicing detection).
     - `--nthread`: Number of threads to use.

3. **Review Output**:
   - rMATS provides results for different types of splicing events:
     - SE: Skipped Exon
     - MXE: Mutually Exclusive Exon
     - A5SS: Alternative 5' Splice Site
     - A3SS: Alternative 3' Splice Site
     - RI: Retained Intron

   Each event type will have separate result files with p-values, FDR values, and splicing inclusion levels.

### Key Points to Consider:
- **Read Length**: Ensure that the `--readLength` matches your actual sequencing read length for accuracy.
- **Strandedness**: If your data is strand-specific, specify it during the rMATS run (`--libType` option).
- **Multi-mapping Reads**: rMATS uses all reads, including multi-mapping reads, which may inflate splicing event counts if not handled properly. It can tolerate these, but carefully consider whether it fits your analysis goals.

If you need any help setting up the parameters or interpreting the rMATS results, let me know!




### Analysis with WCGNA FOR JACHYM
## Nice new tutorial: https://bigomics.ch/blog/introduction-to-wgcna-and-its-applications-in-gene-correlation-network-analysis/

file:///home/aa/pigmentation/tkpAnalysis/petalRNASeq/WCGNA/blockwiseTOM-block.5.RData
file:///home/aa/pigmentation/tkpAnalysis/petalRNASeq/WCGNA/blockwiseTOM-block.4.RData
file:///home/aa/pigmentation/tkpAnalysis/petalRNASeq/WCGNA/blockwiseTOM-block.3.RData
file:///home/aa/pigmentation/tkpAnalysis/petalRNASeq/WCGNA/blockwiseTOM-block.2.RData
file:///home/aa/pigmentation/tkpAnalysis/petalRNASeq/WCGNA/blockwiseTOM-block.1.RData



To perform co-expression analysis using WGCNA (Weighted Gene Co-expression Network Analysis), you typically need gene expression data in the form of a gene expression matrix. Here are the key inputs and details you will need for WGCNA:

#1. Gene Expression Matrix
Rows: Each row represents a gene.
Columns: Each column represents a sample (e.g., biological replicates, time points, conditions).
Values: Each entry in the matrix is the normalized expression level of a gene in a sample (e.g., log-transformed FPKM, TPM, or counts from RNA-seq after normalization).
Format: The gene expression matrix can be in a text format like .csv, .txt, or .tsv. The data needs to be normalized, and batch effects should be corrected if needed. Rows should be labeled by gene IDs or gene names.

#2. Phenotypic Data (Optional but Recommended)
A matrix or table that contains phenotypic information for each sample (e.g., treatment group, time points, or any external traits you want to associate with co-expression modules). Each row should correspond to a sample in the gene expression matrix, and columns should contain the traits or conditions.
Format: This is usually in the form of a .csv, .txt, or .xlsx file.

#3. Soft-thresholding Power Selection
Before proceeding with network construction, WGCNA requires the selection of a soft-thresholding power. This can be done using functions in WGCNA, which evaluate different powers and suggest the best fit based on scale-free topology.
Summary of Input Files
Gene expression matrix: Normalized expression levels (e.g., counts, FPKM, TPM).
Phenotypic data (optional): Information linking samples to traits or conditions.
Soft-thresholding power: To be selected during the analysis.

