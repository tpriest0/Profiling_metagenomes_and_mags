import glob
import os
import pandas as pd

######################
# Generate paths to R1 and R2 reads
######################
SAMPLES = pd.read_csv(config['SAMPLE_NAMES'], header=None)[0].tolist()
os.makedirs(os.path.join(config['WORKING_DIR'], "logs"), exist_ok=True)

######################
# Rules
######################
rule all:
    input:
        pred_gene_marker = expand(os.path.join(config['WORKING_DIR'], "genes", "pred_genes_{sample}.done"), sample=SAMPLES),
        gene_catalog_marker = expand(os.path.join(config['WORKING_DIR'], "gene_catalog", "build_index.done")),        
        scgs_marker = expand(os.path.join(config['WORKING_DIR'], "gene_catalog", "scgs", "extract_scgs.done"))
        
rule pred_genes:
    conda:
        os.path.join(config['WORKFLOW_DIR'], "envs", "predict_genes_env.yaml")
    input:
        contigs=os.path.join(config['CONTIGS_DIR'], "{sample}.contigs.fa")
    output:
        marker=os.path.join(config['WORKING_DIR'], "genes", "pred_genes_{sample}.done"),
        genes_aa=os.path.join(config['WORKING_DIR'], "genes", "{sample}.genes.aa.fa"),
        genes_nt=os.path.join(config['WORKING_DIR'], "genes", "{sample}.genes.nt.fa"),
        genes_gff=os.path.join(config['WORKING_DIR'], "genes", "{sample}.genes.gff"),
        genes_stats=os.path.join(config['WORKING_DIR'], "genes", "{sample}.genes.stats.tsv")
    threads: 
        16
    resources:
        mem=500
    params:
        input_script=os.path.join(config['WORKFLOW_DIR'], "predict_genes_chunked.py")
    log:
        os.path.join(config['WORKING_DIR'], "logs", "pred_genes_{sample}.log")
    shell:
        """
        if [[ ! -s {output.genes_aa} ]] && [[ -s {input.contigs} ]];
        then
            echo "Predicting genes for {wildcards.sample}"
            python {params.input_script} \
                -i {input.contigs} \
                -a {output.genes_aa} \
                -n {output.genes_nt} \
                -g {output.genes_gff} \
                -s {output.genes_stats} \
                -t {threads}
        else
            echo "Input contigs not found for {wildcards.sample}"
        fi

        if [[ -s {output.genes_aa} ]];
        then
            echo "Genes predicted for {wildcards.sample}"
            touch {output.marker}
        else
            echo "Gene prediction failed for {wildcards.sample}"
        fi
        """


def get_aa_files(wildcards):
    """Return list of all predicted amino acid gene FASTAs."""
    return [
        os.path.join(config['WORKING_DIR'], "genes", f"{s}.genes.aa.fa")
        for s in SAMPLES
    ]

def get_nt_files(wildcards):
    """Return list of all predicted nucleotide gene FASTAs."""
    return [
        os.path.join(config['WORKING_DIR'], "genes", f"{s}.genes.nt.fa")
        for s in SAMPLES
    ]

rule concat_genes:
    input:
        aa=get_aa_files,
        nt=get_nt_files
    output:
        marker=os.path.join(config['WORKING_DIR'], "gene_catalog", "concat_genes.done"),
    threads: 
        1
    resources:
        mem=2000
    params:
        all_aa=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.aa.fa")),
        all_nt=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.nt.fa"))
    log:
        os.path.join(config['WORKING_DIR'], "logs", "concat_genes.log")
    shell:
        """
        cat {input.aa} > {params.all_aa}
        cat {input.nt} > {params.all_nt}

        if [[ -s {params.all_aa} ]] && [[ -s {params.all_nt} ]];
        then
            echo "Concatenation of gene files complete"
            touch {output.marker}
        else
            echo "Concatenation of gene files FAILED"
            echo "Exiting..."
            exit 1
        fi
        """

rule build_catalog:
    conda:
        os.path.join(config['WORKFLOW_DIR'], "envs", "mmseqs2_env.yaml")
    input:
        marker=os.path.join(config['WORKING_DIR'], "gene_catalog", "concat_genes.done")
    output:
        marker=os.path.join(config['WORKING_DIR'], "gene_catalog", "build_catalog.done")
    threads: 
        50
    resources:
        mem=1000
    params:
        gene_catalog_dir=os.path.join(config['WORKING_DIR'], "gene_catalog"),
        out_prefix=(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"])),
        all_aa=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.aa.fa")),
        all_nt=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.nt.fa"))
    log:
        os.path.join(config['WORKING_DIR'], "logs", "build_catalog.log")
    shell:
        """
        if [[ -s {params.all_nt} ]];
        then
            echo "Building gene catalog"
            mmseqs easy-linclust {params.all_nt} {params.out_prefix} {params.gene_catalog_dir}/temp \
                --min-seq-id 0.95 \
                -c 0.9 \
                --threads {threads} \
                --split-memory-limit 40G
        else
            echo "Concatenated gene files not found"
        fi

        touch {output.marker}
        """

rule extract_nucleotide_reps:
    input:
        marker=os.path.join(config['WORKING_DIR'], "gene_catalog", "build_index.done")
    output:
        marker=os.path.join(config['WORKING_DIR'], "gene_catalog", "extract_reps.done")
    threads: 
        1
    resources:
        mem=1000
    params:
        out_prefix=(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"])),
        reps_aa=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.reps.aa.fa")),
        reps_nt=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.reps.nt.fa")),
        all_aa=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.aa.fa"))
    log:
        os.path.join(config['WORKING_DIR'], "logs", "extract_nucleotide_reps.log")
    shell:
        """
        if [[ -s {params.out_prefix}-linclust_rep_seq.fasta ]];
        then
            rm {params.out_prefix}-linclust_all_seqs.fasta

            mv {params.out_prefix}-linclust_rep_seq.fasta {params.reps_nt}

            sed -i 's/ //' {params.reps_nt}

            grep '>' {params.reps_nt} | grep -A 1 -wFf - {params.all_aa} | grep -v '^--' > {params.reps_aa}
        else
            echo "MMSeqs gene clustering failed"
        fi

        touch {output.marker}
        """

rule build_index:
    conda:
        os.path.join(config['WORKFLOW_DIR'], "envs", "bwa_env.yaml")
    input:
        marker=os.path.join(config['WORKING_DIR'], "gene_catalog", "extract_reps.done")
    output:
        marker=os.path.join(config['WORKING_DIR'], "gene_catalog", "build_index.done")
    threads: 
        1
    resources:
        mem=40000
    params:
        reps_nt=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.reps.nt.fa"))
    log:
        os.path.join(config['WORKING_DIR'], "logs", "build_index.log")
    shell:
        """
        echo "Building BWA index of {params.ref}" >> {log}

        bwa index {params.ref} 2>> {log}

        touch {output.marker}
        """

rule extract_scgs:
    conda:
         os.path.join(config['WORKFLOW_DIR'], "envs", "fetchmgs_env.yaml")
    input:
        marker=os.path.join(config['WORKING_DIR'], "gene_catalog", "extract_reps.done")
    output:
        marker=os.path.join(config['WORKING_DIR'], "gene_catalog", "scgs", "extract_scgs.done")
    threads: 
        50
    resources:
        mem=500
    params:
        scgs_dir=os.path.join(config['WORKING_DIR'], "gene_catalog", "scgs"),
        out_prefix=config["OUT_PREFIX"],
        reps_aa=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.reps.aa.fa")),
        reps_nt=expand(os.path.join(config['WORKING_DIR'], "gene_catalog", config["OUT_PREFIX"] + ".genes.reps.nt.fa")),
    log:
        os.path.join(config['WORKING_DIR'], "logs", "extract_scgs.log")
    shell:
        """
        echo "Extracting single copy marker genes using FetchMGs"

        if [[ ! -d {params.scgs_dir} ]];
        then
            mkdir -p {params.scgs_dir}
        fi

        fetchMGs extraction {params.reps_aa} gene {params.scgs_dir} -d {params.reps_nt} -t {threads}

        if [[ -s {params.scgs_dir}/{params.out_prefix}.genes.reps.aa.fa.fetchMGs.faa ]] && [[ -s {params.scgs_dir}/{params.out_prefix}.genes.reps.aa.fa.fetchMGs.fna ]] && [[ -s {params.scgs_dir}/{params.out_prefix}.genes.reps.aa.fa.fetchMGs.scores ]];
        then
            grep '>' {params.scgs_dir}/{params.out_prefix}.genes.reps.aa.fa.fetchMGs.faa | sed 's/>//' | sed 's/.COG/\tCOG/' | sed '1i Gene_name\tCOG' > {params.scgs_dir}/COG_marker_to_gene_name_mapping.tsv
            sed -i 's/.COG.*//' {params.scgs_dir}/{params.out_prefix}.genes.reps.aa.fa.fetchMGs.faa
            sed -i 's/.COG.*//' {params.scgs_dir}/{params.out_prefix}.genes.reps.aa.fa.fetchMGs.fna

        else
            echo "FetchMGs failed!"
        fi

        touch {output.marker}
        """
