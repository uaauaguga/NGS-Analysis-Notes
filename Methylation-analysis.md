## Methylation Analysis

### Experimental Methods
  - Capture and sequencing
    - Capture methylated fragments with antibody, the read abundance directionly reflects level of methylation
  - Bisulfite conversion and sequencing
    - Bisulfite convert unmethylated C to U in both strand of DNA, while leaved methylated C unaffected
    - Converted C will sequenced as T

### Reads mappers designed for mapping BS-sequencing data
  - The two complementary DNA strands are no longer complementary after bisulfite treatment, if not all Cs were methylated
    - Directional protocol: Only the top strand and bottom strand are sequenced
    - Undirectional protocol: Converted top strand, strand complementary to converted top strand, converted bottom strand, strand complementary to converted bottom strand are all sequenced with similar likelihood
    - Some spacialized protocol (MCTA-Seq for example) actually selectively sequencing strand complementary to converted top strand and strand complementary to converted bottom strand
  - 3-letters mappers
    - Get C to T converted reference genome
    - Get sequence complementary to genome sequence, convert C to T
    - Build two genome index
    - Directional protocol
      - For each read, convert all C to T
      - Map converted reads to two converted genome, respectively, determine its genome location as the best alignment
    - Undirectional protocol
      - For each read, get its C to T converted sequence and G to A converted sequence
      - Run four way mapping
        - C to T converted reads to C to T converted reference genome
        - C to T converted reads to C to T converted complementary genome 
        - G to A converted reads to C to T converted reference genome
        - G to A converted reads to C to T converted complementary genome
      - Determine genome location by the best alignment
    - Many software implement all of these features <https://academic.oup.com/bib/article/17/6/938/2606438>
    - Bismark is a popular choice <https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html>
      - The user guide is very helpful

### Pipelines
  - MEDIP data analysis: `snakefiles/MEDIP.snakemake`
  - Bisulfite sequencing data analysis: `snakefiles/BS-seq-*e.snakemake`
      - See `config` directory for example configuration
  
