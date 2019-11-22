from os.path import realpath
from os.path import split as pathsplit
import subprocess
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import sys

# Block annoying warnings
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# META
__author__ = "Matt Lawlor"

# SETUP
shell.executable("/bin/bash")
GS = GSRemoteProvider()
FTP = FTPRemoteProvider()
S3 = S3RemoteProvider()
HTTP = HTTPRemoteProvider()

# PARAMS
RUNS = {**config.get("samples",None), **config.get("input_samples",None)}
BT2_MAX_ISIZE = config.get("BT2_MAX_ISIZE",800)
MACS2_SHIFT_SIZE = config.get("MACS2_SHIFT_SIZE",-100)
MACS2_EXTENSION = config.get("MACS2_EXTENSION",200)
MAPQ_CUTOFF = config.get("MAPQ_CUTOFF",30)
GENOME_SIZE = config.get("GENOME_SIZE",None)

# DETERMINE REMOTE OR LOCAL RESOURCE
def determine_resource(path):
    if "gs://" in path:
         return GS.remote(path.replace("gs://",""))
    elif "ftp://" in path:
         return FTP.remote(path)
    elif "s3://" in path:
         return S3.remote(path.replace("s3://",""))
    elif "http://" in path:
         return HTTP.remote(path.replace("http://",""))
    elif "https://" in path:
         return HTTP.remote(path.replace("https://",""))
    else:
        return path

# REMOTE RESOURCES
bt2_idx_paths = [determine_resource(y) for y in config.get("BT2_FILES",None)]
genome_fa = determine_resource(config.get("GENOME_FA",None))
genome_fai = determine_resource(config.get("GENOME_FAI",None))
genome_gzi = determine_resource(config.get("GENOME_GZI",None))
genome_chroi_names = config.get("GENOME_CHR",None).split(",")
genome_bl = determine_resource(config.get("GENOME_BL",None))


rule target:
    """
    Generate peaks and tracks from all samples in the manifest.
    """
    input:
        expand("{s}_peaks.narrowPeak",s=config.get("samples",None)),
        expand("{s}_summits.bed",s=config.get("samples",None)),
        expand("{s}.rpkm.bw",s=config.get("samples",None)),

# ------------------------------------------------------------------------------
# Preproc
# ------------------------------------------------------------------------------

rule concat_fqs:
    """
    Concatenate fastqs for samples with multiple fqs.
    """
    input:
        lambda wc: [determine_resource(x) for x in RUNS[wc.s]["fastq"][wc.end]]
    output:
        temp("fastq/{s}_{end}.fq.gz")
    shell:
        "cat {input} > {output}"

def get_fqs_for_trim(x):
    if (len(RUNS[x]["fastq"].keys()) == 1):
        return ["fastq/{s}_r1.fq.gz"]
    else:
        return ["fastq/{s}_r1.fq.gz", "fastq/{s}_r2.fq.gz"]

def get_fqs_for_aln(x):
    if (len(RUNS[x]["fastq"].keys()) == 1):
        return ["fastq/{s}_r1.trimmed.fq.gz"]
    else:
        return ["fastq/{s}_r1.trimmed.fq.gz", "fastq/{s}_r2.trimmed.fq.gz"]

def get_proper_ended_fastp_call(x):
    fqs = get_fqs_for_trim(x)
    call = ""
    if len(fqs) == 1:
        call = "--in1 {r1}".format(r1=fqs[0].format(s=x))
    else:
        call = "--in1 {r1} --in2 {r2}".format(r1=fqs[0].format(s=x), r2=fqs[1].format(s=x))
    if RUNS[x].get("phred",33) == 64:
        call += " --phred64"
    return(call)



def get_proper_ended_fastp_out(x):
    fqs = get_fqs_for_aln(x)
    if len(fqs) == 1:
        return "--out1 {r1}".format(r1=fqs[0].format(s=x))
    else:
        return "--out1 {r1} --out2 {r2}".format(r1=fqs[0].format(s=x), r2=fqs[1].format(s=x))

rule trim_se:
    """
    Trim single-end fastqs with fastp.
    """
    input:
        fq = lambda wc: get_fqs_for_trim(wc.s),
    output:
        r1 = "fastq/{s}_r1.trimmed.fq.gz",
        html = "fastq/{s}_fastp.html",
        json = "fastq/{s}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.s),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.s)
    conda:
        "envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp {params.call_in} "
        "{params.call_out} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.s}_fastp"

rule trim_pe:
    """
    Trim paired-end fastqs with fastp.
    """
    input:
        fq = lambda wc: get_fqs_for_trim(wc.samp)
    output:
        r1 = "fastq/{s}_r1.trimmed.fq.gz",
        r2 = "fastq/{s}_r2.trimmed.fq.gz",
        html = "fastq/{s}_fastp.html",
        json = "fastq/{s}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.s),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.s)
    conda:
        "envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp {params.call_in} "
        "{params.call_out} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.s}_fastp"

def get_proper_ended_aln_call(x):
    fqs = get_fqs_for_aln(x)
    if len(fqs) == 1:
        return "-U {r1}".format(r1=fqs[0].format(s=x))
    else:
        return "-1 {r1} -2 {r2}".format(r1=fqs[0].format(s=x), r2=fqs[1].format(s=x))


#http://biolearnr.blogspot.com/2017/11/snakemake-using-inputoutput-values-in.html
rule align_bt2:
    """
    Align reads with bowtie2.
    """
    input:
        fqs = lambda wc: rules.trim_se.output if (len(RUNS[wc.s]["fastq"].keys()) == 1) else rules.trim_pe.output,
        idx = bt2_idx_paths,
    output:
        temp("{s}.raw.sam")
    conda:
        "envs/bowtie2.yaml"
    singularity:
        "docker://quay.io/biocontainers/bowtie2:2.3.5--py37he860b03_0"
    threads:
        4
    params:
        idx_pfx = config.get("BT2_IDX_PFX",None),
        reads = lambda wc: get_proper_ended_aln_call(wc.s)
    shell:
        "bowtie2 --phred33 -p {threads} " # fastp always returns phred33
        "--no-discordant --no-unal -k 1 "
        "-x {params.idx_pfx} {params.reads} -S {output}"

rule raw_to_cram:
    """
    Convert alignments in sam format to cram.
    """
    input:
        crm = rules.align_bt2.output,
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        temp("{s}.raw.cram")
    threads:
        4
    conda:
        "envs/samtools.yaml"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    shell:
        "samtools sort -O cram -@ {threads} "
        "--output-fmt-option lossy_names=1,level=9,store_md=0,store_nm=0 "
        "-o {output} --reference {input.fa} {input.crm}"

rule blacklist_filter_reads:
    """
    Remove blacklistlisted reads.
    """
    input:
        crm=rules.raw_to_cram.output,
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
        bl=genome_bl,
        crai="{s}.raw.cram.crai"
    output:
        temp("{s}.bl.cram")
    conda:
        "envs/bedtools.yaml"
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.29.0--hc088bd4_3"
    threads:
        1
    shell:
        "CRAM_REFERENCE={input.fa} "
        "bedtools intersect -v -a {input.crm} -b {input.bl} > {output}"

rule fix_mate_info:
    """
    Update mate info in aux tags and coord sort.
    """
    input:
        crm="{s}.bl.nsrt.cram",
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        temp("{s}.fixm.cram")
    conda:
        "envs/samtools.yaml"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    threads: 1
    shell:
        "samtools fixmate -m --reference {input.fa} {input.crm} - | "
        "samtools sort -O cram "
        "--output-fmt-option lossy_names=1,level=9,store_md=0,store_nm=0 "
        "-o {output} --reference {input.fa}"

rule clean_reads:
    """
    Filter reads by MAPQ, canonical chromosomes, pcr dups.
    """
    input:
        crm="{s}.fixm.cram",
        crai="{s}.fixm.cram.crai",
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        "{s}.clean.cram"
    params:
        chr=genome_chroi_names,
        mapq=MAPQ_CUTOFF,
    conda:
        "envs/samtools.yaml"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    threads:
        1
    shell:
        "samtools view -u -q {params.mapq} "
        "--reference {input.fa} {input.crm} {params.chr} | "
        "samtools markdup -r "
        "--output-fmt-option lossy_names=1,level=9,store_md=0,store_nm=0 "
        "--reference {input.fa} - {output}"

# ------------------------------------------------------------------------------
# Peak calling
# ------------------------------------------------------------------------------

def get_alns_for_sample(x):
    alns = ["{s}.clean.bam".format(s=x)]

    if config["samples"][x].get("input",None):
        alns = alns + ["{s}.clean.bam".format(s=config["samples"][x].get("input",None))]

    return alns

def get_proper_macs2_call(x):
    alns = get_alns_for_sample(x)
    call = "-t {ip} ".format(ip=alns[0])
    if len(alns) == 2:
        call += "-c {bg} ".format(bg = alns[1])

    if len(RUNS[x]["fastq"].keys()) != 2:
        call += "-f BAM"
    else:
        call += "-f BAMPE"

    return call



rule call_peaks:
    """
    Call peaks using macs2 and shifting the reads to center on the cut site.
    """
    input:
        lambda wc: get_alns_for_sample(wc.s)
    output:
        sum="{s}_summits.bed",
        np="{s}_peaks.narrowPeak"
    params:
        gs=GENOME_SIZE,
        call = lambda wc: get_proper_macs2_call(wc.s)
    conda:
        "envs/macs2.yaml"
    singularity:
        "docker://quay.io/biocontainers/macs2:2.2.5--py37h516909a_0"
    shadow:
        "shallow"
    threads:
        2
    shell:
        "macs2 callpeak {params.call} " # TODO right now this only takes the 1st read in the pair, which is find for now...
        "--keep-dup all "
        "-g {params.gs} -n {wildcards.s} --call-summits; "

# ------------------------------------------------------------------------------
# viz
# ------------------------------------------------------------------------------

rule make_bigwigs:
    """
    Generate rpkm-normalized bigwigs from filtered bams.
    """
    input:
        crm="{s}.clean.cram",
        crai="{s}.clean.cram.crai"
    output:
        "{s}.rpkm.bw"
    threads:
        4
    conda:
        "envs/deeptools.yaml"
    singularity:
        "docker://quay.io/biocontainers/deeptools:3.3.1--py_0"
    shell:
        "bamCoverage -b {input.crm} "
        "-p {threads} "
        "--outFileName {output} "
        "--outFileFormat bigwig "
        "--extendReads 200 "
        "--binSize 50 --smoothLength 150 "
        "--verbose --normalizeUsing RPKM"

# ------------------------------------------------------------------------------
# Generics
# ------------------------------------------------------------------------------

rule index_cram:
    """
    Create a *.crai for fast random access of your cram files.
    """
    input:
        "{file}.cram"
    output:
        temp("{file}.cram.crai")
    conda:
        "envs/samtools.yaml"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    shell:
        "samtools index {input}"

rule index_bam:
    """
    Create a *.bai for fast random access of your cram files.
    """
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    conda:
        "envs/samtools.yaml"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    shell:
        "samtools index {input}"


rule cram_to_bam:
    """
    A generic rule for converting cram to bam when downstream tools require bam.
    """
    input:
        cram="{file}.cram",
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        temp("{file}.bam")
    conda:
        "envs/samtools.yaml"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    shell:
        "samtools view -b {input.cram} -o {output} -T {input.fa}"

rule nsort_cram:
    """
    Generic rule for sorting a cram by read name.
    """
    input:
        crm="{file}.cram",
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        temp("{file}.nsrt.cram")
    conda:
        "envs/samtools.yaml"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    shell:
        "samtools sort -n -O cram {input.crm} --reference {input.fa} -o {output} "
        "--output-fmt-option lossy_names=1,level=9,store_md=0,store_nm=0"

# ------------------------------------------------------------------------------
# HELP
# ------------------------------------------------------------------------------

rule help:
    shell:
        """
        echo '
        ===== google cloud =====

        # install prereqs
        pip install kubernetes
        gcloud components install kubectl

        # set up cluster variables
        CLUSTER_NAME=snk-cl2
        NODES=8
        ZONE=us-central1-a
        REMOTE=GS
        PREFIX=archibald
        MACHINE_TYPE=n1-standard-2

        # initialize cluster
        gcloud container clusters create $CLUSTER_NAME \
            --num-nodes=$NODES \
            --scopes storage-rw \
            --machine-type=$MACHINE_TYPE \
            --zone $ZONE

        # register cluster info
        gcloud container clusters get-credentials $CLUSTER_NAME --zone $ZONE


        snakemake --kubernetes --use-conda \
            --default-remote-provider $REMOTE \
            --default-remote-prefix $PREFIX \
            --latency-wait 300 \
            --jobs 8 \
            --verbose \
            --debug-dag

        # shut down your cluster
        gcloud container clusters delete $CLUSTER_NAME --zone $ZONE

        '
        """
