
# Snakefile
import glob
import os
import pandas as pd

######################
# Generate paths to R1 and R2 reads
######################
SAMPLES = pd.read_csv(config['SAMPLE_NAMES'], header=None)[0].tolist()
os.makedirs(os.path.join(config['WORKING_DIR'], "logs"), exist_ok=True)

def read1(sample):
    return os.path.join(config['READ_DIR'], f"{sample}", f"{sample}_1.fq.gz")

def read2(sample):
    return os.path.join(config['READ_DIR'], f"{sample}" ,f"{sample}_2.fq.gz")

######################
# Rules
######################
rule all:
    input:
        map_reads_markers = expand(os.path.join(config['WORKING_DIR'], "map_reads", "map_reads_{sample}.done"), sample=SAMPLES),
        filter_mapped_reads_markers = expand(os.path.join(config['WORKING_DIR'], "filter_mapped_reads", "filter_{sample}.done"), sample=SAMPLES),
        bam_to_coverage_markers = expand(os.path.join(config['WORKING_DIR'], "coverage", "bam_to_coverage_{sample}.done"), sample=SAMPLES)

rule map_reads:
    conda:
        os.path.join(config['WORKFLOW_DIR'], "envs", "bwa_env.yaml")
    input:
        r1=lambda wc: read1(wc.sample),
        r2=lambda wc: read2(wc.sample)
    output:
        marker=os.path.join(config['WORKING_DIR'], "map_reads", "map_reads_{sample}.done")
    threads: 
        12
    resources:
        mem=10000
    params:
        ref=config['REFERENCE'],
        map_reads_dir=os.path.join(config['WORKING_DIR'], "map_reads"),
        bam=os.path.join(config['WORKING_DIR'], "map_reads", "{sample}.bam")
    log:
        os.path.join(config['WORKING_DIR'], "logs", "map_reads_{sample}.log")
    shell:
        """
        mkdir -p $(dirname {params.bam})

        echo "Mapping reads for {wildcards.sample}" > {log} 2>&1
        echo "R1: {input.r1}" >> {log}
        echo "R2: {input.r2}" >> {log}

        bwa mem -a -t {threads} {params.ref} {input.r1} {input.r2} 2>> {log} | \
            samtools view -@ {threads} -bh - > {params.bam} 2>> {log}

        if [[ -s {params.bam} ]]; 
        then
            echo "BAM file created for {wildcards.sample}" >> {log}
            touch {output.marker}
        else
            echo "Creation of BAM file FAILED for {wildcards.sample}" >> {log}
            exit 1
        fi
        """

rule filter_mapped_reads:
    conda:
        os.path.join(config['WORKFLOW_DIR'], "envs", "bwa_env.yaml")
    input:
        marker=os.path.join(config['WORKING_DIR'], "map_reads", "map_reads_{sample}.done")
    output:
        marker=os.path.join(config['WORKING_DIR'], "filter_mapped_reads", "filter_{sample}.done")
    threads: 
        4
    resources:
        mem=12000
    params:
        bam=os.path.join(config['WORKING_DIR'], "map_reads", "{sample}.bam"),
        filtbam=os.path.join(config['WORKING_DIR'], "filter_mapped_reads", "{sample}.filtered.bam")
    log:
        os.path.join(config['WORKING_DIR'], "logs", "filter_mapped_reads_{sample}.log")
    shell:
        """

        mkdir -p $(dirname {params.filtbam}) >> {log}

        # Keep only reads with:
        #  - identity >= 0.95  (1 - NM / aligned_M)
        #  - coverage >= 0.8   (aligned_M / read_length)
        samtools view -@ {threads} -h {params.bam} 2>> {log} | \
        awk 'BEGIN{{OFS="\\t"}}
             /^@/ {{print; next}}
             {{
                 cigar=$6; nm=0; aln=0;

                 # extract NM:i tag (edit distance)
                 for(i=12;i<=NF;i++) if($i ~ /^NM:i:/) {{nm=substr($i,6); break}}

                 # sum all M operations in CIGAR as aligned length
                 tmp=cigar;
                 while (match(tmp, /[0-9]+M/)) {{
                     mlen=substr(tmp,RSTART,RLENGTH-1);
                     aln+=mlen;
                     tmp=substr(tmp,RSTART+RLENGTH);
                 }}

                 # read length from SEQ field
                 qlen=length($10);

                 if (aln>0 && qlen>0) {{
                     pid = 1 - nm/aln;       # identity
                     cov = aln/qlen;         # fraction of read aligned
                     if (pid >= 0.95 && cov >= 0.80) print;
                 }}
             }}' 2>> {log} | \
        samtools view -@ {threads} -bh - 2>> {log} | \
        samtools sort -@ {threads} -o {params.filtbam} - 2>> {log}

        if [[ -s {params.filtbam} ]]; 
        then
            samtools index -@ {threads} {params.filtbam} >> {log}
            touch {output.marker}
        else
            echo "Filtered BAM not created for {wildcards.sample}" >> {log}
            exit 1
        fi
        """

rule bam_to_coverage:
    conda:
        os.path.join(config['WORKFLOW_DIR'], "envs", "python_env.yaml")
    input:
        marker=os.path.join(config['WORKING_DIR'], "filter_mapped_reads", "filter_{sample}.done")
    output:
        marker=os.path.join(config['WORKING_DIR'], "coverage", "bam_to_coverage_{sample}.done")
    threads: 
        12
    resources:
        mem=8000
    params:
        coverage_script=os.path.join(config['WORKFLOW_DIR'], "calc_gene_coverage_from_bam.py"),
        filtbam=os.path.join(config['WORKING_DIR'], "filter_mapped_reads", "{sample}.filtered.bam"),
        coverage=os.path.join(config['WORKING_DIR'], "coverage", "{sample}.coverage.tsv")
    log:
        os.path.join(config['WORKING_DIR'], "logs", "bam_to_coverage_{sample}.log")
    shell:
        """
        if [[ -s {params.filtbam} ]];
        then
            echo "Filtering of bam file complete"
            echo "Calculating coverage statistics"

            python {params.coverage_script} -i {params.filtbam} -o {params.coverage} -t {threads}
        else
            echo "Filtering of bam file failed for {wildcards.sample}"
        fi

        if [[ -s {params.coverage} ]];
        then
            echo "Calculation of coverage from bam file complete."
            touch {output.marker}
        else
            echo "Calculation of coverage from bam file FAILED"
        fi
        """
