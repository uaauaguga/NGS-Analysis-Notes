# Notes for human genome annotation

- See `/BioII/lulab_b/jinyunfan/data/human-genome-annotations`

## Human genome assembly

- sequences in hg38 assembly from different sources should be identical, but sequence ids can be different

- make sure that sequence id in reference genome is consistent with that in gff3/gtf file
  - you can simply use both reference genome and gff3 file from gencode
  - this repo <https://github.com/dpryan79/ChromosomeMappings> provide some chromosome id cross mappings

```bash
# ebi-gencode
# here use ebi-gencode for further processing
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz

# other sources
# ucsc
# wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
# ensembl
# 1 2 3 instead of chr1 chr2 chr3 ...
# wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# ncbi
# ncbi refseq sequence id, instead of chr1 chr2 chr3 ...
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
```

## Download gene annotations

### Gencode

```{bash}
# download gencode annotation
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
```
- gencode annotation
  - gene id: gene_id, `^ENSG\d{11}\.\d$` 
  - gene biotype annotate in field 9: 
  - protein_coding	19955
  - lncRNA	16888
  - miRNA	1879
  - misc_RNA	2212
  - Mt_rRNA	2
  - Mt_tRNA	22
  - Pseudogene
    - pseudogene	15
    - transcribed_processed_pseudogene	502
    - transcribed_unitary_pseudogene	143
    - transcribed_unprocessed_pseudogene	950
    - translated_processed_pseudogene	2
    - translated_unprocessed_pseudogene	1
    - polymorphic_pseudogene	49
    - processed_pseudogene	10163
    - unitary_pseudogene	98
    - unprocessed_pseudogene	2614
  - rRNA	47
  - rRNA_pseudogene	497
  - IG genes
    - IG_C_gene	14
    - IG_C_pseudogene	9
    - IG_D_gene	37
    - IG_J_gene	18
    - IG_J_pseudogene	3
    - IG_pseudogene	1
    - IG_V_gene	145
    - IG_V_pseudogene	187
  - TR genes
    - TR_C_gene	6
    - TR_D_gene	4
    - TR_J_gene	79
    - TR_J_pseudogene	4
    - TR_V_gene	106
    - TR_V_pseudogene	33
  - vault_RNA	1
  - ribozyme	8
  - scaRNA	49
  - scRNA	1
  - snoRNA	943
  - snRNA	1901
  - sRNA	5
  - TEC	1056

  - tx id: transcript_id, `^ENST\d{11}\.\d$`

### refseq 

- Use rRNA annotation in refseq

```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
zcat source/annotation/GRCh38_latest_genomic.gff.gz > gff3/refseq.GRCh38_latest_genomic.gff3 
```

### mitranscriptome 

- Seems not under active maintenance

- **The raw gtf is on hg19 coordinate**

```bash
  wget http://mitranscriptome.org/download/mitranscriptome.gtf.tar.gz
```

- gene id: gene_id, `^G\d{6}$` 
- gene biotype: tcat
  - lncrna: 63615
  - mixed_read_through: 4188
  - protein_coding: 19922
  - pseudogene: 7502
  - tucp: 3755

- tx id: transcript_id, `^T\d{6}$`

### tRNA annotation

- tRNA-SE scan result provided by gencode

```bash
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.tRNAs.gtf.gz

# http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/
```

### miRNA annotation

- miRbase

```bash
wget --recursive --no-parent ftp://mirbase.org/pub/mirbase/22.1/
```

- human miRNA gff: `source/annotation/mirbase.org/pub/mirbase/22.1/genomes/hsa.gff3`

- primary transcript
  - id: `MI\d{7}`
  - name: `hsa-mir-\d{4}` 
- mature miRNA 
  - derived from primary transcript
  - five prime arm
  - id: `MIMAT\d{7}`
    - name: `hsa-miR-\d{4}-5p `
  - three prime arm
    - id: `MIMAT\d{7}`
    - name: `hsa-miR-\d{4}-3p`

- piRNA annotation
  
  - piRBase <http://www.regulatoryrna.org/database/piRNA/download.html>
  
```bash
wget http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/fasta/hsa.fa.gz
wget http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/data/piR_hsa.txt.gz
wget http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/bed/hsa.bed.gz
```

### circRNA annotation

- circatlas
  
- http://circatlas.biols.ac.cn/
  
- circbase
  - http://www.circbase.org/cgi-bin/downloads.cgi

```bash
  wget http://www.circbase.org/download/hscore/hscores_human_gencode27.tar.gz
  wget http://www.circbase.org/download/hsa_hg19_circRNA.bed
  wget http://www.circbase.org/download/hsa_hg19_circRNA.txt
  wget http://www.circbase.org/download/hg19_circID_to_name.txt
  wget http://www.circbase.org/download/human_hg19_circRNAs_putative_spliced_sequence.fa.gz
```

### repeats
- UCSC table browser RepeatMasker track
  - Go to UCSC table browser <http://genome.ucsc.edu/cgi-bin/hgTables>
  - set group to Repeats, track to RepeatMasker
  - type file name to save as
  - click "get output"


## Processing of annotations

- Make sure the version of genome annotation is consistent with reference genome

- Get mRNA from gencodev38 gff3 file

```bash
cat gff3/gencode.v38.annotation.gff3 | awk 'BEGIN{FS="\t";OFS="\t";}$0~"#"{print;next }$0~"gene_type=protein_coding"{print}' > gff3/gencode.v38.mRNA.gff3
```

- Get mRNA trasncript

```bash
gffread -g fasta/hg38/genome.fa -s fasta/hg38/chrom.size -W -M -F -G -A -O -E -w fasta/transcript/gencode.v38.mRNA.fa -d fasta/transcript/gencode.v38.mRNA.fa.collapsed.info gff3/gencode.v38.mRNA.gff3 
```

- Prepare mRNA bed12 format transcript annotation 

```bash
 gffread -W -M -F -G -A -E --bed gff3/gencode.v38.mRNA.gff3 > bed12/gencode.v38.mRNA.bed
```


- Prepare mRNA gene/exon/intron bed file

```bash
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.mRNA.gff3 -b bed/gencode.v38.mRNA.gene.bed -f gene
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.mRNA.gff3 -b bed/gencode.v38.mRNA.exon.bed -f exon 
  bedtools subtract -s -a bed/gencode.v38.mRNA.gene.bed -b bed/gencode.v38.mRNA.exon.bed > bed/gencode.v38.mRNA.intron.bed 
```

- Get lncRNA from gencodev38 gff3 file

```bash
cat gff3/gencode.v38.annotation.gff3 | awk 'BEGIN{FS="\t";OFS="\t";}$0~"#"{print;next }$0~"gene_type=lncRNA"{print}' > gff3/gencode.v38.lncRNA.gff3
```

- Prepare mRNA bed12 format transcript annotation

```bash
gffread -W -M -F -G -A -E --bed gff3/gencode.v38.lncRNA.gff3 > bed12/gencode.v38.lncRNA.bed 
```

- Get transcript sequences with reference genome and gff3 annotation

```bash
gffread -g fasta/hg38/genome.fa -s fasta/hg38/chrom.size -W -M -F -G -A -O -E -w fasta/transcript/gencode.v38.lncRNA.fa -d fasta/transcript/gencode.v38.lncRNA.fa.collapsed.info gff3/gencode.v38.lncRNA.gff3 
```

- Get gene/exon/intron bed file of lncRNA

```bash
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.lncRNA.gff3 -b bed/gencode.v38.lncRNA.gene.bed -f gene
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.lncRNA.gff3 -b bed/gencode.v38.lncRNA.exon.bed -f exon
  bedtools subtract -s -a bed/gencode.v38.lncRNA.gene.bed -b bed/gencode.v38.lncRNA.exon.bed > bed/gencode.v38.lncRNA.intron.bed
```

- Convert mitranscriptome to hg38 coordinate

```bash
  zcat source/annotation/mitranscriptome.gtf/mitranscriptome.v2.gtf.gz > gtf/mitranscriptome.v2.hg19.gtf
  wget  http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
  CrossMap.py gff source/annotation/hg19ToHg38.over.chain.gz gtf/mitranscriptome.v2.hg19.gtf gtf/mitranscriptome.v2.hg38.gtf
```

- Processing lncRNA in mitranscriptome

```bash
# get lncRNA gtf
cat gtf/mitranscriptome.v2.hg38.gtf  | grep 'tcat "lncrna"' > gtf/mitranscriptome.v2.hg38.lncRNA.gtf

# get gene intervals
cat gtf/mitranscriptome.v2.hg38.lncRNA.gtf | grep $'\ttranscript\t' | awk 'BEGIN{FS="\t";OFS="\t";}{split($9,attrs,";");split(attrs[5],fields," ");tx_id=fields[2];gsub(/"/, "", tx_id);print $1,$4-1,$5,tx_id,".",$7;}' > bed/mitranscriptome.v2.hg38.lncRNA.gene.bed 

# get exon intervals
cat gtf/mitranscriptome.v2.hg38.lncRNA.gtf | grep $'\texon\t' | awk 'BEGIN{FS="\t";OFS="\t";}{split($9,attrs,";");split(attrs[6],tx_fields," ");tx_id=tx_fields[2];gsub(/"/, "", tx_id); split(attrs[1],exon_fields," ");exon_id=exon_fields[2];gsub(/"/, "", exon_id); print $1,$4-1,$5,tx_id"-"exon_id,".",$7;}' > bed/mitranscriptome.v2.hg38.lncRNA.exon.bed

# get intron intervals
bedtools subtract -s -a bed/mitranscriptome.v2.hg38.lncRNA.gene.bed -b bed/mitranscriptome.v2.hg38.lncRNA.exon.bed > bed/mitranscriptome.v2.hg38.lncRNA.intron.bed 
```


- Genomic tRNA 
  - gencode tRNA prediction
```bash
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.tRNAs.gff3 -b bed/gencode.v38.tRNA.bed -f tRNA --name ID,gene_type
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa -bed bed/gencode.v38.tRNA.bed -s  | sed 's/::/ /g' > fasta/transcript/gencode.v38.tRNA.fa  
```

- Mt_tRNA
  - 22 mt tRNA, each gene has one tx, each tx has one exon
```bash
  cat gff3/gencode.v38.annotation.gff3 | grep 'gene_type=Mt_tRNA' > gff3/gencode.v38.Mt_tRNA.gff3
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.Mt_tRNA.gff3 -b bed/gencode.v38.Mt_tRNA.bed -f gene 
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa -bed bed/gencode.v38.Mt_tRNA.bed -s  | sed 's/::/ /g' > fasta/transcript/gencode.v38.Mt_tRNA.fa  
```

- Genomic rRNA
  - gencode v38 rRNA annotation is incomplete (no 16S and 28S)
  - Here we use refseq

```bash
  # NC_012920.1 is refseq id of mitochondrial
  cat gff3/refseq.GRCh38_latest_genomic.gff3 | awk 'BEGIN{FS="\t";OFS="\t";}$0~"^#"{print;next;}$3=="rRNA"{print}' | grep -v 'NC_012920.1' > gff3/refseq.rRNA.gff3
  scripts/gff3-to-bed.py -g gff3/refseq.rRNA.gff3 -f rRNA -b bed/refseq.rRNA.bed -n Name,product
```

  - Get rRNA sequences

```bash
  wget -P source https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_RefSeq2UCSC.txt
  wget -P source https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_UCSC2gencode.txt
  scripts/rename-chrom.py -i bed/refseq.rRNA.bed -o bed/refseq.rRNA.UCSC.seq_id.bed  -n source/GRCh38_RefSeq2UCSC.txt
  # skip NW_018654708.1 and NW_021160023.1 
  scripts/rename-chrom.py -i bed/refseq.rRNA.UCSC.seq_id.bed -o bed/refseq.rRNA.gencode.seq_id.bed  -n source/GRCh38_UCSC2gencode.txt
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa -bed bed/refseq.rRNA.gencode.seq_id.bed -s | sed 's/::/ /g' > fasta/transcript/refseq.rRNA.fa 
  
```

- Mt rRNA
  - Only 2 mt rRNA
  - gencode and refseq annotation are identical 

```bash
  cat gff3/gencode.v38.annotation.gff3 | awk 'BEGIN{FS="\t";OFS="\t";}$0~"#"{print;next;}$0~"gene_type=Mt_rRNA;"{print}' > gff3/gencode.v38.Mt_rRNA.gff3
  scripts/gff3-to-bed.py -g gff3/gencode.v38.Mt_rRNA.gff3 -b bed/gencode.v38.Mt_rRNA.bed -f gene 
  # refseq has identical mt rRNA annotation
  # cat gff3/refseq.GRCh38_latest_genomic.gff3 | awk 'BEGIN{FS="\t";OFS="\t";}$3=="rRNA"{print}' | grep 'NC_012920.1' > gff3/refseq.Mt_rRNA.gff3 
```
  - get mt rRNA sequence

  ```bash
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa -bed bed/gencode.v38.Mt_rRNA.bed -s | sed 's/::/ /g' > fasta/transcript/gencode.v38.Mt_rRNA.fa 
  ```

- miRNA
  - miRbase

  ```bash
  cp source/annotation/mirbase.org/pub/mirbase/22.1/genomes/hsa.gff3 gff3/mirbase.22.1.hsa.gff3
  cat gff3/mirbase.22.1.hsa.gff3 | awk 'BEGIN{FS="\t";OFS="\t";}$0~"#"{print;next }$3=="miRNA_primary_transcript"{print}' > gff3/mirbase.22.1.pre-miRNA.gff3
  scripts/gff3-to-bed.py -g gff3/mirbase.22.1.pre-miRNA.gff3 -b bed/mirbase.22.1.pre-miRNA.bed -n Name -f miRNA_primary_transcript 
  cat gff3/mirbase.22.1.hsa.gff3 | awk 'BEGIN{FS="\t";OFS="\t";}$0~"#"{print;next }$3=="miRNA"{print}' > gff3/mirbase.22.1.mature-miRNA.gff3
  scripts/gff3-to-bed.py -g gff3/mirbase.22.1.mature-miRNA.gff3 -b bed/mirbase.22.1.mature-miRNA.bed -n Name -f miRNA 
  ```
  
  - Get miRNA sequence

```bash
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa -bed bed/mirbase.22.1.pre-miRNA.bed -s | sed 's/::/ /g' > fasta/transcript/mirbase.22.1.pre-miRNA.fa 
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa -bed bed/mirbase.22.1.mature-miRNA.bed -s | sed 's/::/ /g' > fasta/transcript//mirbase.22.1.mature-miRNA.fa
```

- piRNA
  
  - piRbase annotation

```bash
  zcat source/annotation/piRBase2.0/hsa.bed.gz >  bed/piRBase2.0.bed
  zcat source/annotation/piRBase2.0/hsa.fa.gz > fasta/transcript/piRBase2.0.fa   
```


- snRNA

```bash
  cat gff3/gencode.v38.annotation.gff3  | awk 'BEGIN{FS="\t";OFS="\t";}$0~"^#"{print;next;}$9~"gene_type=snRNA;"{print}' > gff3/gencode.v38.snRNA.gff3
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.snRNA.gff3 -b bed/gencode.v38.snRNA.bed -f gene
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa  -bed bed/gencode.v38.snRNA.bed -s | sed 's/::/ /g' > fasta/transcript/gencode.v38.snRNA.fa 
```


- snoRNA

```bash
  cat gff3/gencode.v38.annotation.gff3  | awk 'BEGIN{FS="\t";OFS="\t";}$0~"^#"{print;next;}$9~"gene_type=snoRNA;"{print}' > gff3/gencode.v38.snoRNA.gff3
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.snoRNA.gff3 -b bed/gencode.v38.snoRNA.bed -f gene
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa  -bed bed/gencode.v38.snoRNA.bed -s | sed 's/::/ /g' > fasta/transcript/gencode.v38.snoRNA.fa 
```

- srpRNA
  - Annotated as misc_RNA in gencode

```bash
  cat gff3/gencode.v38.annotation.gff3 | awk 'BEGIN{FS="\t";OFS="\t";}$0~"^#"{print;next;}$9~/gene_name=RN7SL.+/{print}' > gff3/gencode.v38.srpRNA.gff3
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.srpRNA.gff3 -b bed/gencode.v38.srpRNA.bed -f gene
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa  -bed bed/gencode.v38.srpRNA.bed -s | sed 's/::/ /g' > fasta/transcript/gencode.v38.srpRNA.fa 
```

- Y RNA
  - Annotated as misc_RNA in gencode

```bash
  cat gff3/gencode.v38.annotation.gff3 | awk 'BEGIN{FS="\t";OFS="\t";}$0~"^#"{print;next;}$9~/gene_name=Y_RNA;/{print}' > gff3/gencode.v38.Y_RNA.gff3
  scripts/gff3-to-bed.py -g  gff3/gencode.v38.Y_RNA.gff3 -b bed/gencode.v38.Y_RNA.bed -f gene
  bedtools getfasta -name+ -fi fasta/hg38/genome.fa  -bed bed/gencode.v38.Y_RNA.bed -s | sed 's/::/ /g' > fasta/transcript/gencode.v38.Y_RNA.fa 
```

- circRNA
  - circbase: putative back-spliced sequences: `fasta/circRNA/cricbase.hg19.fa`
  - circatlas
  ```bash
  # Convert back spliced sequence to fasta format
  cat human_sequence_v2.0.txt | awk 'BEGIN{FS="\t";}$3!="partial"{print ">"$2;print $3;}' >  fasta/circRNA/circatlas.v2.0.fa
  ```

## Notes on file format

- gtf and gff3 format
  - 1 based coordinate (similar to array in R/matlab), closed interval [start, end]
  - different lines / records have a hierarchical structure
    - Each genome features can access its parent by attribute in field 9
    - gencode gff3 annotation
      - there is a gene -> transcript -> exon hierarchy
      - For protein coding gene, the hierarchy is gene -> transcript -> (exon,five_prime_UTR,start_codon,CDS,stop_codon,three_prime_UTR)
    - mirbase gff3 annotation
      - miRNA_primary_transcript -> miRNA hierarchy

- bed format and bed12 format
  - 0 based coordinate (similar to array in C/python), half open interval [start,end)
  - different lines are independent
  - bed format describe general genomic interval
  - bed12 format is a specific case of bed format
    - each line could describe the structure of a transcript
    - genomic positions of exons that make up this trasncript
    - where is CDS (for mRNA transcript)

- For more information, see
  - <http://genome.ucsc.edu/FAQ/FAQformat.html>
  - <https://bedtools.readthedocs.io/en/latest/content/general-usage.html>