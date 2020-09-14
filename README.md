# TCGA Microbe Analysis

## Start from bam file

- STAR's analyzing pipeline  
- https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
- Bam file generated in steps 4 is the downloaded bam files
  - --outSAMunmapped Within was specified, and this bam files actually contains all information in bam files
- Reference used TCGA
  - https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
  - hg38 reference genome
  - ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
  
  - Sequence decoy
    - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz
  - Virus decoy
    - https://gdc.cancer.gov/files/public/file/GRCh83.d1.vd1_virus_decoy.txt


```bash
# STAR-2.4.2a

### For users with access to the ICGC pipeline:

python star_align.py \
--genomeDir <star_index_path> \
--FastqFileIn <input_fastq_path> \
--workDir <work_dir> \
--out <output_bam> \
--genomeFastaFiles <reference> \
--runThreadN 8 \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--limitBAMsortRAM 0 \
--alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--twopass1readsN -1 \
--sjdbOverhang 100 \
--outSAMstrandField intronMotif \
--outSAMunmapped Within

### For users without access to the ICGC pipeline:

### Step 1: Building the STAR index.*

STAR
--runMode genomeGenerate
--genomeDir <star_index_path>
--genomeFastaFiles <reference>
--sjdbOverhang 100
--sjdbGTFfile <gencode.v22.annotation.gtf>
--runThreadN 8

### Step 2: Alignment 1st Pass.

STAR
--genomeDir <star_index_path>
--readFilesIn <fastq_left_1>,<fastq_left2>,... <fastq_right_1>,<fastq_right_2>,...
--runThreadN <runThreadN>
--outFilterMultimapScoreRange 1
--outFilterMultimapNmax 20
--outFilterMismatchNmax 10
--alignIntronMax 500000
--alignMatesGapMax 1000000
--sjdbScore 2
--alignSJDBoverhangMin 1
--genomeLoad NoSharedMemory
--readFilesCommand <bzcat|cat|zcat>
--outFilterMatchNminOverLread 0.33
--outFilterScoreMinOverLread 0.33
--sjdbOverhang 100
--outSAMstrandField intronMotif
--outSAMtype None
--outSAMmode None

### Step 3: Intermediate Index Generation.

STAR
--runMode genomeGenerate
--genomeDir <output_path>
--genomeFastaFiles <reference>
--sjdbOverhang 100
--runThreadN <runThreadN>
--sjdbFileChrStartEnd <SJ.out.tab from previous step>

### Step 4: Alignment 2nd Pass.

STAR
--genomeDir <output_path from previous step>
--readFilesIn <fastq_left_1>,<fastq_left2>,... <fastq_right_1>,<fastq_right_2>,...
--runThreadN <runThreadN>
--outFilterMultimapScoreRange 1
--outFilterMultimapNmax 20
--outFilterMismatchNmax 10
--alignIntronMax 500000
--alignMatesGapMax 1000000
--sjdbScore 2
--alignSJDBoverhangMin 1
--genomeLoad NoSharedMemory
--limitBAMsortRAM 0
--readFilesCommand <bzcat|cat|zcat>
--outFilterMatchNminOverLread 0.33
--outFilterScoreMinOverLread 0.33
--sjdbOverhang 100
--outSAMstrandField intronMotif
--outSAMattributes NH HI NM MD AS XS
--outSAMunmapped Within
--outSAMtype BAM SortedByCoordinate
--outSAMheaderHD @HD VN:1.4
--outSAMattrRGline <formatted RG line provided by wrapper>
```

## Get Unmapped Reads

```bash
# -f 13: both mate should unmapped
# -F 0x200: reads should pass quality filter
# samtools sort -n : sort by read id
samtools view -h -f 13 -F 0x200 input.bam | samtools sort -n | samtools view -b > unmapped.bam
```
## Get Reads Aligned to Viral Decoy

```bash
bin/get-viral-reads.py -ib input.bam -ii input.bam.bai -ob viral-decoy.bam
```

## Quantify reads aligned to viral decoy

```bash
bin/count_reads.py count_transcript -i viral-decoy.bam -o viral-decoy.counts.txt
```

## Bam to fastq

```bash
# bedtools bamtofastq -i unmapped.bam -fq unmapped_1.fastq -fq2 unmapped_2.fastq 
# Do not automatically compress output fastq file 
# samtools fastq unmapped.bam [-N] -1 unmapped_1.fastq.gz -2 unmapped_2.fastq.gz 
# Generate exactly same result of bedtools bamtofastq if -N is specified, or read mate (/1 or /2 after read id) will not automatically appended

## Convert bam file of unmapped reads and viral decoy to fastq
samtools fastq viral-decoy.bam -1 viral-decoy_1.fastq.gz -2 viral-decoy_2.fastq.gz 
samtools fastq unmapped.bam -1 unmapped_1.fastq.gz -2 unmapped_2.fastq.gz 

## Concatenate unmapped reads and viral decoy
zcat viral-decoy_1.fastq.gz unmapped_1.fastq.gz | gzip -c > exogenous_1.fastq.gz
zcat viral-decoy_2.fastq.gz unmapped_2.fastq.gz | gzip -c > exogenous_2.fastq.gz
```

## Classification

- Two database 
  - /BioII/lulab_b/jinyunfan/data/minikraken2_v1_8GB
    - Downloaded from ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v1_8GB_201904.tgz
  - /BioII/lulab_b/jinyunfan/data/standard-db
    - kraken2's standard database, plus fungi database

```bash
kraken2 --db {db}  --threads {threads} --unclassified-out {unclassified}  --report {report} --paired --use-names  exogenous_1.fastq.gz exogenous_2.fastq.gz
```

