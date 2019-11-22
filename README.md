# csp2

## prerequisites

1. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

2. Set up your conda installation channels properly.

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. Install [snakemake](https://snakemake.readthedocs.io/en/stable/)

```bash
conda install snakemake
```

4. *OPTIONAL:* [Install singularity](https://sylabs.io/singularity/)

## getting started

1. Clone or download the repo.
2. Modify `config.yaml` to include your data.
3. From the repo repository run:
```bash
snakemake --configfile config.yaml --cores 2 --keep-remote --directory work/
```
This will run the workflow using at most 2 threads in the `work/` directory.
Additionally it will keep any remote resources downloaded during the run in case
debugging is necessary.

It is recommended to run with one of the following flags to take care of dependency
management automatically:
```bash
--use-singularity

# or

--use-conda

```

## features

Paths in your config file can be local paths, remote URIs including FTP,
google cloud, amazon aws (if you have appropriately configured amazon or google
accounts). The most frequent use-case would be providing a path to a fastq
on ENA's ftp server. If you have trouble with this note that ENA has technical
difficulty on occasion, so it may make sense to try again later before trying
more intense debugging.

Multiple fastqs can be automatically concatenated per sample - just be aware that
the order for paired-end fastqs must be conserved for each read pair.

Input files are listed in the config file under "input_files." The same input
file can be used for multiple ChIP libraries where appropriate, so you can specify
the name of each input library for each ChIP.

```
 JURKAT_H3K27AC_01:
   phred: 33
   input: JURKAT_WCE_01 # <- this should be the name of entry in 'input_files'
    ... etc ...
```

The template config.yaml file provided in the repo points to bowtie2 indices
and genome files hosted in a private google cloud bucket, which I will likely keep
private for a while because it costs me money to ingress/egress data.
However the exact same indices are publically available in compressed archives on the
[bowtie2 website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
The corresponding required fasta files can be derived from the indices with the following
commands:

```bash
bowtie2-inspect {index_prefix} | bgzip -c -i > genome.fa.gz
samtools faidx genome.fa.gz
```

Once you have extracted the indices and reconsituted the genome fasta just modify
the config file to point to your local copies.


## to get help

To print a description of every rule:

```bash
snakemake -l
```

To print a walkthrough of a google cloud run:

```bash
snakemake help
```
