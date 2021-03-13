## cfDNA sequencing data analysis
### Environment setting
- Add `/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin` to `$PATH`
- Install manta for structural variation calling
  ```bash
  conda create -n manta-env
  conda activate manta-env
  conda install -c bioconda manta
  ``` 
- Run `snakefiles/cfDNA-pipeline.snakemake` with in manta-env
- `/BioII/lulab_b/jinyunfan/projects/diagnostic-marker/cfDNA/Snakefile` is the same file
