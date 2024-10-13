# Bioinformatics_Analysis
Exome Sequencing Analysis Pipeline  A pipeline for paired-end exome data analysis from SRA, including QC, alignment (BWA), variant calling (GATK), and annotation. All steps, commands, and results are documented with runtime metrics.

# Exome Analysis Project

1)exome sequencing (SRR30834370)         
2)exome sequencing (SRR30834374)

## Data Download

**Sample 1 (SRR30834370)**

prefetch --output-directory /home/vaishnavi/Bioinfo_Analysis/Samples SRR30834370
fastq-dump --split-files /home/vaishnavi/Bioinfo_Analysis/Samples/SRR30834370/SRR30834370.sra

**Sample 2 (SRR30834374)**

prefetch --output-directory /home/vaishnavi/Bioinfo_Analysis/Samples SRR30834374
fastq-dump --split-files /home/vaishnavi/Bioinfo_Analysis/Samples/SRR30834374/SRR30834374.sra


**Subsetting and Calculating Statistics
Downsampling:**

seqkit head -n 10000 sample1_R1.fastq > data/sample1_R1_subset.fastq
seqkit head -n 10000 sample1_R2.fastq > data/sample1_R2_subset.fastq
seqkit head -n 10000 sample2_R1.fastq > data/sample2_R1_subset.fastq
seqkit head -n 10000 sample2_R2.fastq > data/sample2_R2_subset.fastq


**FASTQ statistics:**

seqkit stats data/sample1_R1_subset.fastq data/sample1_R2_subset.fastq data/sample2_R1_subset.fastq data/sample2_R2_subset.fastq > results/fastq_statistics.txt


**Quality Control (QC)**

fastqc data/sample1_R1_subset.fastq -o results/
fastqc data/sample1_R2_subset.fastq -o results/
fastqc data/sample2_R1_subset.fastq -o results/
fastqc data/sample2_R2_subset.fastq -o results/


**Genome Download and Indexing**

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip hg38.fa.gz

bwa index hg38.fa


**Alignment**

**Align the paired-end reads:**

bwa mem -t 8 hg38.fa sample1_R1_subset.fastq sample1_R2_subset.fastq > sample1_aligned.sam

bwa mem -t 8 hg38.fa sample2_R1_subset.fastq  sample2_R2_subset.fastq  > sample2_aligned.sam


**Converting SAM to BAM and Sorting**

samtools view -Sb sample1_aligned.sam | samtools sort -o sample1_sorted.bam

samtools view -Sb sample2_aligned.sam | samtools sort -o sample2_sorted.bam


**BAM Indexing**

samtools index sample1_sorted.bam

samtools index sample2_sorted.bam


**Alignment statistics**

samtools flagstat sample1_sorted.bam > sample1_alignment_stats.txt                         **or**           samtools stats sample1_sorted.bam > sample1_stats.txt

samtools flagstat sample2_sorted.bam > sample2_alignment_stats.txt                   ** or**               samtools stats sample2_sorted.bam > sample2_stats.txt


**Mark Duplicates**

**Sample1:**

samtools markdup is an efficient way to identify and mark duplicate reads in BAM files:

samtools sort -n sample1_sorted.bam -o sample1_name_sorted.bam

samtools fixmate -m sample1_name_sorted.bam sample1_fixmate.bam

samtools sort sample1_fixmate.bam -o sample1_positions_sorted.bam

samtools markdup sample1_positions_sorted.bam sample1_markdup.bam

samtools flagstat sample1_markdup.bam > sample1_markdup_stats.txt

samtools index sample1_markdup.bam

**Sample2:**

samtools sort -n sample2_sorted.bam -o sample2_name_sorted.bam

samtools fixmate -m sample2_name_sorted.bam sample2_fixmate.bam

samtools sort sample2_fixmate.bam -o sample2_positions_sorted.bam

samtools markdup sample2_positions_sorted.bam sample2_markdup.bam

samtools flagstat sample2_markdup.bam > sample2_markdup_stats.txt

samtools index sample2_markdup.bam


**Variant Calling with GATK**

Prepare the reference genome:

gatk CreateSequenceDictionary -R data/hg38.fa

samtools faidx data/hg38.fa

**Adding or Replacing Read Groups with GATK**

**Sample1:**

gatk AddOrReplaceReadGroups \
    -I sample1_sorted.bam \
    -O sample1_rg.bam \
    --RGID sample1 \
    --RGLB lib1 \
    --RGPL illumina \
    --RGPU unit1 \
    --RGSM sample1

**Sample2:**

gatk AddOrReplaceReadGroups \
    -I sample2_sorted.bam \
    -O sample2_rg.bam \
    --RGID sample2 \
    --RGLB lib2 \
    --RGPL illumina \
    --RGPU unit2 \
    --RGSM sample2


**HaplotypeCaller**

**Sample1:**

./gatk HaplotypeCaller \
    -R hg38.fa \
    -I ../sample1_rg.bam \
    -O ../sample1_variants.vcf \
    --emit-ref-confidence GVCF \
   --sample-name sample1

**Sample2:**

./gatk HaplotypeCaller \
    -R hg38.fa \
    -I ../sample2_rg.bam \
    -O ..sample2_variants.vcf \
    --emit-ref-confidence GVCF \
    --sample-name sample2

**SNP and Indel Separation**

**Sample1:**

awk 'BEGIN {OFS="\t"} {if ($5 ~ /,/) $5 = split($5, a, ",") ? a[1] : $5; print}' sample1_variants.vcf | grep -v "NON" > sample1_variants_out.vcf

**Sample2:**

awk 'BEGIN {OFS="\t"} {if ($5 ~ /,/) $5 = split($5, a, ",") ? a[1] : $5; print}' sample2_variants.vcf | grep -v "NON" > sample2_variants_out.vcf

**Filtering SNPs & Filtering Indels**

**Sample1:**

awk 'BEGIN {OFS="\t"} /^#/ {print > "sample1_variants_out_snps.vcf"} length($4) == 1 && length($5) == 1 {print > "sample1_variants_out_snps.vcf"}' sample1_variants_out.vcf

awk 'BEGIN {OFS="\t"} /^#/ {print > "sample1_variants_out_indels.vcf"} length($4) != 1 || length($5) != 1 {print > "sample1_variants_out_indels.vcf"}' sample1_variants_out.vcf

**Sample2:**

awk 'BEGIN {OFS="\t"} /^#/ {print > "sample2_variants_out_snps.vcf"} length($4) == 1 && length($5) == 1 {print > "sample2_variants_out_snps.vcf"}' sample2_variants_out.vcf

awk 'BEGIN {OFS="\t"} /^#/ {print > "sample2_variants_out_indels.vcf"} length($4) != 1 || length($5) != 1 {print > "sample2_variants_out_indels.vcf"}' sample2_variants_out.vcf


**Comparative Analysis**

**#Common Variants:**

awk 'BEGIN {OFS="\t"} FNR==NR && !/^#/ {key=$1":"$2":"$4":"$5; sample1[key]=$0; next} /^#/ {print > "common_variants.vcf"} !/^#/ {key=$1":"$2":"$4":"$5; if (key in sample1) {print sample1[key], $NF >> "common_variants.vcf"}}' sample1_variants_out.vcf sample2_variants_out.vcf

**Common SNPs:**

awk 'BEGIN {OFS="\t"} FNR==NR && !/^#/ {key=$1":"$2":"$4":"$5; sample1[key]=$0; next} /^#/ {print > "common_snps.vcf"} !/^#/ {key=$1":"$2":"$4":"$5; if (key in sample1) {print sample1[key], $NF >> "common_snps.vcf"}}' sample1_variants_out_snps.vcf sample2_variants_out_snps.vcf

**Common Indels:**

awk 'BEGIN {OFS="\t"} FNR==NR && !/^#/ {key=$1":"$2":"$4":"$5; sample1[key]=$0; next} /^#/ {print > "common_indels.vcf"} !/^#/ {key=$1":"$2":"$4":"$5; if (key in sample1) {print sample1[key], $NF >> "common_indels.vcf"}}' sample1_variants_out_indels.vcf sample2_variants_out_indels.vcf

**Variant Annotation**

for i in $(find . -type f -name "*.vcf"); do vep -i $i -o ${i/.vcf/_vep.vcf} --cache --dir_cache /mnt/d/Downloads/VEP_database/ --refseq --fasta hg38.fa --everything --format vcf; done

