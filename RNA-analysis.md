## Notes for RNA analysis

### Starting point
- Copy the scripts, genome  and tools directory in /BioII/lulab_b/jinyunfan/projects/diagnostic-marker/RNA-analysis into an empty folder
- Prepare required input files
  - cleaned fastq files (no adapter, remove spike in if any, remove potential UniVec contamination, remove rRNA)
  - genome aligned bam files, sorted by coordinate, with duplication removed
    - Unaligned reads of STAR mapping in exSEEK may out of order due to a bug in STAR_2.5.3a_modified
    - use scripts/checkFastqPairing.py to check whether circRNA_1.fastq.gz and circRNA_2.fastq.gz is out of order
    - If its the case, may rerun mapping, and fix unmapped reads of each step with repair.sh in bbmap ( see Chimeric RNA quantification section )
  - reads unmapped to circRNA


### Add reads group to bam file
  - gatk require these fields
  ```{bash}
  gatk AddOrReplaceReadGroups --java-options -Xmx4G \
        --INPUT {input.bam} --OUTPUT {output.bam} -SO  coordinate \
        --RGLB library --RGPL illumina --RGPU HiSeq2000 --RGSM {sample_id} > {log} 2>&1
  samtools index {output.bam}
  ```

### SNP calling
  - Split Reads by "N" in cigar string, input is bam file with reads group
  - Get red of spliced junctions
  ```{bash}
  gatk SplitNCigarReads --java-options -Xmx4G --input {input.bam} \
  --output {output.bam} \
  --create-output-bam-index -R genome/fasta/hg38.fa \
  --tmp-dir tmp 2> {log}
  ```
  - Here the base quality recalibration step is skipped, as it is time consuming and makes little difference
  - Run HaplotypeCaller, using bam output of SplitNCigarReads as input
  ```{bash}
  gatk HaplotypeCaller --java-options -Xmx4G -R genome/fasta/hg38.fa \
        -I {input.bam} -O {output.vcf} \
        --tmp-dir tmp 2>  {log}
  ``` 
  - Filter the vcf file
  ```{bash}
  gatk VariantFiltration --java-options -Xmx4G \
        -R genome/fasta/hg38.fa -V {input.vcf} \
        -window 35 -cluster 3 \
        --filter-name FS20 -filter "FS > 20.0" \
        --filter-name QD2 -filter "QD < 2.0" \
        --filter-name DP10 -filter "DP < 10.0" \
        --filter-name QUAL20 -filter "QUAL < 20.0" -O {output.vcf} 2> {log} 
  ```
  
### RNA editing analysis
  - long RNA editing analysis
  - Input: bam with read group
  ```{bash}
  gatk ASEReadCounter --java-options -Xmx4G --input {input.bam} --variant genome/vcf/REDIportal.vcf.gz  --reference genome/fasta/hg38.fa --output-format TABLE 2> {log} | gzip -c > {output.count}
  ```

### ASE analysis
  - Replace the vcf file in RNA editing analysis with COSMIC and dbSNP
  ```{bash}
  gatk ASEReadCounter --java-options -Xmx4G --input {input.bam} --variant genome/vcf/dbSNP.vcf.gz.vcf.gz --reference genome/fasta/hg38.fa --output-format TABLE 2> {log} | gzip -c > {output.dbSNP-count}
  gatk ASEReadCounter --java-options -Xmx4G --input {input.bam} --variant genome/vcf/COSMIC.vcf.gz --reference genome/fasta/hg38.fa --output-format TABLE 2> {log} | gzip -c > {output.COSMIC-count}
  ```


### Alternative Splicing Analysis
  - rMATS
  - Get bam path for position and negative set, see output/test/splicing/bam-path for example
  - According to https://groups.google.com/g/rmats-user-group/c/qq5ryHPWrIc, --libType is not actually used if using bam of STAR alignment as input, this parameter actually passed to tophat for alignment
    - That means for bam input or STAR alignment, different --libType generate same result
  - For each run, positive and negative samples should be specified. For very large sample size, the program may crash.
  - Not friendly for large scale analysis or samples with mautiple class
  ```{bash}
  python2 /BioII/lulab_b/jinyunfan/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 {pos_path} --b2 {neg_path} --gtf genome/gtf/gencode.v27.annotation.gtf --od {outdir} -t paired  --libType fr-firststrand --readLength 150 
  python2 /BioII/lulab_b/jinyunfan/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 output/test/splicing/bam-path/PDAC.txt  --b2 output/test/splicing/bam-path/HD.txt --gtf genome/gtf/gencode.v27.annotation.gtf --od output/test/splicing/rmats -t paired  --libType fr-firststrand --readLength 150
  ``` 
  - Summarize the result, assign a unique ID to each alternative splicing event
  - The input pos and neg sample ids should same as that provided for rMATS
  - splicing_type if one of MXE A3SS A5SS SE RI 
  ```{bash}
  python3 scripts/summarize-splicing.py --input {rmats_outdir}/{splicing_type}.MATS.JC.txt  --outdir {outdir} --type {splicing_type} --method JC  --pos {pos}  --neg  {neg}
  #python3 scripts/summarize-splicing.py --input output/test/splicing/rmats/SE.MATS.JC.txt  --type SE --method JC --pos metadata/test/PDAC.txt --neg metadata/test/HD.txt  --outdir output/test/splicing/matrix
  ``` 


### APA
  - The DaPar software asks for wig input, while actually it needs a bedgraph , suggested by its example
    - See /BioII/lulab_b/jinyunfan/software/dapars-modified/DaPars_Test_Dataset/Condition_B_chrX.wig
  - To generate the annotation
  ```{bash}
  # Prepare transcript information in bed12 format
  # /BioII/lulab_b/jinyunfan/software/gffread-0.11.4.Linux_x86_64/gffread
  gffread genome/gtf/gencode.v27.annotation.gtf --bed > genome/bed/gencode.v27.annotation.bed
  cat genome/bed/gencode.v27.annotation.bed | cut -f 1,2,3,4,5,6,7,8,9,10,11,12 > genome/bed/gencode.v27.annotation.12.bed
  # Preapre transcript to gene mapping
  cat genome/gtf/gencode.v27.annotation.gtf | awk 'BEGIN{OFS="\t";FS="\t";}$3=="transcript"{print $9}' | awk 'BEGIN{OFS="\t";FS=";";}{split($1,gene," ");split($2,tx," ");print tx[2],gene[2]}' | sed 's/"//g' > genome/tx2gene.txt
  # Generate 3 prime UTR annotation. Use python2 here
  python /BioII/lulab_b/jinyunfan/software/dapars-modified/src/DaPars_Extract_Anno.py -b genome/bed/gencode.v27.annotation.12.bed -s genome/tx2gene.txt -o genome/bed/gencode.v27.annotation.UTR3p.bed
  ```
  - To generate the "wig" file:
  ```{bash}
  bedtools genomecov -ibam {input.bam} -bg -split | sort-bed - > {output.wig}
  ``` 
  - Follow the documentation to run the script
  ```{bash}
  # Prepare config file
  scripts/prepareDaParConfig.py -p metadata/test/PDAC.txt -n metadata/test/HD.txt -i output/test/APA/wig -o output/test/APA/dapar --config output/test/APA/config.txt #python3 
  # Run DaPar
  python2 /BioII/lulab_b/jinyunfan/software/dapars-modified/src/DaPars_main.py output/test/APA/config.txt
  # Summarize result
  scripts/parseDapar.py -i output/test/APA/dapar/result_All_Prediction_Results.txt -c output/test/APA/config.txt -l output/test/APA/matrix/long.txt -s output/test/APA/matrix/short.txt -p output/test/APA/matrix/PDUI.txt
  ```


### Alternative promoter analysis
  - Quantify expression level of all transcript isoforms of a gene
  - lib type: 
    - If paired end, ISF for forward stranded, ISR for reverse stranded, IU for unstranded
    - Also has an automatic strandness detection feature: --libType A
    - See https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype for detail
  ```{bash}
  #salmon 1.3.0
  salmon quant --threads {threads} --libType {libtype} -i genome/salmon-index -1 {input.fastq1} -2 {input.fastq2} --validateMappings --gcBias -o {outdir}
  ```
  - Summarize salmon output
  ```{bash}
  scripts/prepareTxQuantMat.py -i output/test/salmon -o output/test/TPM-by-tx.txt
  ```
  - Assign expression of transcript with similar TSS to same promoter
  ```{bash}
  scripts/getPromoterActivity.py -i output/test/TPM-by-tx.txt -o output/test/TPM-by-promoter.txt 
  ``` 


### TE element analysis
  - Analysis with salmon TE
  ```{bash}
  #/BioII/lulab_b/jinyunfan/anaconda3/envs/quantification/bin/python
  tools/SalmonTE/SalmonTE.py quant --reference=hs --num_threads=6 --outpath=output/{sample_id}  {fastq1} {fastq2}
  ``` 


### Chimeric RNA quantification
  - Align reads that unaligned to circRNA to curated jucntion sequences with STAR
  ```{bash}
  STAR --genomeDir genome/chimera-star-index \
            --readFilesIn {input.reads1} {input.reads2} \
            --runThreadN 4 \
            --outFileNamePrefix {output_prefix} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --readFilesCommand gzip -d -c \
            --outSAMmultNmax 1 \
            --seedPerWindowNmax 20
  # STAR --genomeDir genome/chimera-star-index \
  #           --readFilesIn data/test/unmapped/SRR10080507_1.fastq.gz data/test/unmapped/SRR10080507_2.fastq.gz \
  #           --runThreadN 4 \
  #          --outFileNamePrefix output/test/chimeric-RNA/SRR10080507/ \
  #           --outSAMtype BAM Unsorted \
  #           --outReadsUnmapped Fastx \
  #           --readFilesCommand gzip -d -c \
  #           --outSAMmultNmax 1 \
  #           --seedPerWindowNmax 20
  # To fix the bug of STAR_2.5.3a_modified 
  /BioII/lulab_b/jinyunfan/software/bbmap/repair.sh in={output_prefix}/Unmapped.out.mate1 in2={output_prefix}/Unmapped.out.mate1 out=unmapped_1.fastq.gz out2=unmapped_2.fastq.gz overwrite=t 
  ```


### Metagenomic classification of unmapped sequences
  - Input: reads unaligned to curated chimeras
  - Database: /Share2/home/lulab/jinyunfan/data/kraken2db/standard-db
  ```{bash}
  kraken2 --db {input.database} --paired --threads {threads} --unclassified-out {params.unclassified} --report {output.report}  --use-names  {input.fastq1} {input.fastq2}  >  {output.assignment} 2> {log} 
  ```
