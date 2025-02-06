configfile: "config.yaml"
# SAMPLES = ["A", "B", "C"]
rule all:
    input:
        "plots/quals.svg"
rule download_data:
    output:
        "snakemake-tutorial-data-v5.4.5.tar.gz"
    shell:
        "curl -L -o {output} https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz"
rule extract_data:
    input:
        "snakemake-tutorial-data-v5.4.5.tar.gz"
    output:
        "data/genome.fa"
    shell:
        """
        mkdir -p data
        tar -xzf {input} -C data --strip-components 2
        """
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
