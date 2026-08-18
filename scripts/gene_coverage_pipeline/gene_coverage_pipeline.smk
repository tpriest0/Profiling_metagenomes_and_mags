# Snakefile
import glob
import os
import pandas as pd

######################
# Generate paths to R1 and R2 reads
######################
SAMPLES = pd.read_csv(config['SAMPLE_NAMES'], header=None)[0].tolist()
os.makedirs(os.path.join(config['WORKING_DIR'], "logs"), exist_ok=True)

R1_TAGS = ["_1", "_R1"]
R2_TAGS = ["_2", "_R2"]
EXTS = [".fa.gz", ".fq.gz", ".fasta.gz", ".fastq.gz"]


def find_read(sample, read_dir, tags):
    """
    Find exactly one read file for a given sample and read direction.
    Raises a clear error if zero or multiple matches are found.
    """
    sample_dir = os.path.join(read_dir, sample)
    candidates = []

    for tag in tags:
        for ext in EXTS:
            pattern = os.path.join(sample_dir, f"{sample}{tag}*{ext}")
            candidates.extend(glob.glob(pattern))

    candidates = sorted(set(candidates))

    if len(candidates) == 0:
        raise ValueError(
            f"[ERROR] No reads found for sample '{sample}' with tags {tags} "
            f"in {sample_dir}"
        )
    if len(candidates) > 1:
        raise ValueError(
            f"[ERROR] Multiple read files found for sample '{sample}' with tags {tags}:\n"
            + "\n".join(candidates)
        )

    return candidates[0]


def read1(sample):
    return find_read(sample, config["READ_DIR"], R1_TAGS)


def read2(sample):
    return find_read(sample, config["READ_DIR"], R2_TAGS)

######################
# Rules
######################
rule all:
    input:
        build_index_marker = os.path.join(config['OUTPUT_DIR'], "gene_catalog", "build_index.done"),
        map_reads_markers = expand(os.path.join(config['OUTPUT_DIR'], "map_reads", "map_reads_{sample}.done"), sample=SAMPLES),
        filter_mapped_reads_markers = expand(os.path.join(config['OUTPUT_DIR'], "filter_mapped_reads", "filter_{sample}.done"), sample=SAMPLES),
        bam_to_coverage_markers = expand(os.path.join(config['OUTPUT_DIR'], "coverage", "bam_to_coverage_{sample}.done"), sample=SAMPLES),
        normalise_coverage_markers = expand(os.path.join(config['OUTPUT_DIR'], "coverage_normalised", "normalise_{sample}.done"), sample=SAMPLES)

rule bwa_index:
    conda:
        "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/tpriest/projects/giulia_project/scripts/gene_coverage_pipeline/coverage_env.yml"
    input:
        ref=config['REFERENCE']
    output:
        amb=config['REFERENCE'] + ".amb",
        ann=config['REFERENCE'] + ".ann",
        bwt=config['REFERENCE'] + ".bwt",
        pac=config['REFERENCE'] + ".pac",
        sa=config['REFERENCE'] + ".sa",
        marker=os.path.join(config['OUTPUT_DIR'], "gene_catalog", "build_index.done")
    resources:
        mem=32000
    log:
        os.path.join(config['OUTPUT_DIR'], "logs", "bwa_index.log")
    shell:
        """
        bwa index {input.ref}
        """

rule map_reads:
    conda:
        os.path.join(config['WORKFLOW_DIR'], "envs", "bwa_env.yaml")
    input:
        r1=lambda wc: read1(wc.sample),
        r2=lambda wc: read2(wc.sample),
        ref=config['REFERENCE'],
        bwt=config['REFERENCE'] + ".bwt",
        marker=os.path.join(config['OUTPUT_DIR'], "gene_catalog", "build_index.done")
    output:
        marker=os.path.join(config['OUTPUT_DIR'], "map_reads", "map_reads_{sample}.done")
    threads: 
        12
    resources:
        mem=6000
    params:
        map_reads_dir=os.path.join(config['OUTPUT_DIR'], "map_reads"),
        bam=os.path.join(config['OUTPUT_DIR'], "map_reads", "{sample}.bam")
    log:
        os.path.join(config['OUTPUT_DIR'], "logs", "map_reads_{sample}.log")
    shell:
        """
        mkdir -p {params.map_reads_dir}

        echo "Mapping reads for {wildcards.sample}" >> {log}
        echo "R1: {input.r1}" >> {log}
        echo "R2: {input.r2}" >> {log}

        if [[ -s {input.r1} ]] && [[ -s {input.r2} ]];
        then

            bwa mem -a -t {threads} {input.ref} {input.r1} {input.r2} &>> {log} | \
            samtools view -@ {threads} -bh - > {params.bam} &>> {log}
        
        else
            echo "Could not find input reads for {wildcards.sample}."
            exit 1
        fi

        if [[ -s {params.bam} ]];
        then
            echo "Bam file created for {wildcards.sample}" >> {log}
            touch {output.marker}
        else
            echo "Creation of bam file FAILED for {wildcards.sample}" >> {log}
            exit 1
        fi
        """

rule filter_mapped_reads:
    conda:
        os.path.join(config['WORKFLOW_DIR'], "envs", "bwa_env.yaml")
    input:
        marker=os.path.join(config['OUTPUT_DIR'], "map_reads", "map_reads_{sample}.done")
    output:
        marker=os.path.join(config['OUTPUT_DIR'], "filter_mapped_reads", "filter_{sample}.done")
    threads: 
        2
    resources:
        mem=12000
    params:
        bam=os.path.join(config['OUTPUT_DIR'], "map_reads", "{sample}.bam"),
        filtbam=os.path.join(config['OUTPUT_DIR'], "filter_mapped_reads", "{sample}.filtered.bam")
    log:
        os.path.join(config['OUTPUT_DIR'], "logs", "filter_mapped_reads_{sample}.log")
    shell:
        """
        mkdir -p $(dirname {params.filtbam})

        # Keep only reads with:
        #  - identity >= 0.95  (1 - NM / aligned_M)
        #  - coverage >= 0.8   (aligned_M / read_length)

        echo "Filtering read alignments to retain those with >=0.95 identity and >=0.8 horizontal coverage for {wildcards.sample}" >> {log}

        samtools view -@ {threads} -h {params.bam} &>> {log} | \
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
        samtools view -@ {threads} -bh - &>> {log} | \
        samtools sort -@ {threads} -o {params.filtbam} - &>> {log}

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
        os.path.join(config['WORKFLOW_DIR'], "envs", "bwa_env.yaml")
    input:
        marker=os.path.join(config['OUTPUT_DIR'], "filter_mapped_reads", "filter_{sample}.done")
    output:
        marker=os.path.join(config['OUTPUT_DIR'], "coverage", "bam_to_coverage_{sample}.done")
    threads: 
        1
    resources:
        mem=40000
    params:
        bam=os.path.join(config['OUTPUT_DIR'], "map_reads", "{sample}.bam"),
        filtbam=os.path.join(config['OUTPUT_DIR'], "filter_mapped_reads", "{sample}.filtered.bam"),
        raw_coverage=os.path.join(config['OUTPUT_DIR'], "coverage", "{sample}.samtools_coverage.tsv"),
        coverage=os.path.join(config['OUTPUT_DIR'], "coverage", "{sample}.coverage.tsv")
    log:
        os.path.join(config['OUTPUT_DIR'], "logs", "bam_to_coverage_{sample}.log")
    shell:
        """
        ### Use samtools coverage command to return coverage statistics for each reference
        ### We also use the -d 0, to include even zero-covered positions, and we manually
        ### specify the --ff fields to include, so that we can remove the default SECONDARY,
        ### which ignores secondary alignments

        if [[ -s {params.filtbam} ]];
        then
            echo "Filtering of read alignments complete"

            echo "Calculating coverage statistics based on filtered read alignments" >> {log}

            samtools coverage -d 0 --ff UNMAP,QCFAIL,DUP {params.filtbam} -o {params.raw_coverage} &>> {log}

            awk 'BEGIN{{OFS="\\t"}}
                NR==1 {{print "Gene","Gene_length","Num_bases_covered","Prop_bases_covered","Mean_depth","Mean_depth_per_kbp"; next}}
                {{gene_length=$3-$2+1; print $1,gene_length,$5,$6,$7,$7*(1000/gene_length)}}' \
                {params.raw_coverage} > {params.coverage} &>> {log}
        else
            echo "Filtered bam file not found" >> {log}
            echo "Exiting..." >> {log}
            exit 1
        fi

        if [[ -s {params.coverage} ]];
        then
            echo "Calculation of coverage statistics for {wildcards.sample} complete" >> {log}
            echo "Cleaning up working directory"
            rm -rf {params.bam} {params.raw_coverage}
            touch {output.marker}
        else
            echo "Calculation of coverage statistics for {wildcards.sample} FAILED" >> {log}
            exit 1
        fi
        """

rule normalise_coverage_by_scgs:
    conda:
        os.path.join(config['WORKFLOW_DIR'], "envs", "bwa_env.yaml")
    input:
        marker=os.path.join(config['OUTPUT_DIR'], "coverage", "bam_to_coverage_{sample}.done"),
        scg_mapping=config['SCG_MAPPING']
    output:
        marker=os.path.join(config['OUTPUT_DIR'], "coverage_normalised", "normalise_{sample}.done")
    threads:
        1
    resources:
        mem=8000
    params:
        breadth_threshold=95,
        coverage=os.path.join(config['OUTPUT_DIR'], "coverage", "{sample}.coverage.tsv"),
        normalised_coverage=os.path.join(config['OUTPUT_DIR'], "coverage_normalised", "{sample}.coverage.normalised.tsv")
    log:
        os.path.join(config['OUTPUT_DIR'], "logs", "normalise_coverage_{sample}.log")
    shell:
        """
        echo "Estimating the number of genomes sequenced and calculating genome-normalised gene coverage for {wildcards.sample}" >> {log}

        if [[ -s {params.coverage} ]] && [[ -s {input.scg_mapping} ]];
        then
            awk -v OFS="\\t" -v thresh={params.breadth_threshold} '
            BEGIN {{
                # Load single-copy marker gene -> COG mapping (global, one file
                # for the whole gene catalogue -- not sample-specific)
                while ((getline line < "{input.scg_mapping}") > 0) {{
                    split(line, f, "\\t")
                    if (f[1] != "Gene_name") scg_cog[f[1]] = f[2]
                }}
                close("{input.scg_mapping}")
            }}
            NR==1 {{ next }}  # skip coverage file header
            $4 >= thresh {{
                # Task 1: retain only genes with Prop_bases_covered >= 95
                n++
                gene[n]=$1; genelen[n]=$2; covbases[n]=$3; propcov[n]=$4; meandepth[n]=$5; mdpk[n]=$6

                # Task 2: sum length-normalised coverage per marker COG
                # (only genes surviving the filter above contribute)
                if ($1 in scg_cog) {{
                    cog_sum[scg_cog[$1]] += $6
                }}
            }}
            END {{
                # Task 3: Number_of_genomes_sequenced = mean of the per-COG
                # sums, denominator = number of marker COGs actually detected
                # (observed-only -- undetected COGs are simply not counted)
                n_cogs = 0; total = 0
                for (c in cog_sum) {{ total += cog_sum[c]; n_cogs++ }}
                genomes_sequenced = (n_cogs > 0) ? total / n_cogs : "NA"

                # Task 5: write filtered table with Gene_length retained and
                # the new normalised column added
                print "Gene","Gene_length","Num_bases_covered","Prop_bases_covered","Mean_depth","Mean_depth_per_kbp","Mean_depth_per_genome"
                for (i=1; i<=n; i++) {{
                    # Task 4: normalise every retained gene's coverage by
                    # Number_of_genomes_sequenced
                    copy_num = (genomes_sequenced != "NA" && genomes_sequenced > 0) ? mdpk[i] / genomes_sequenced : "NA"
                    print gene[i], genelen[i], covbases[i], propcov[i], meandepth[i], mdpk[i], copy_num
                }}
            }}' {params.coverage} > {params.normalised_coverage} &>> {log}
        else
            echo "Coverage table for {wildcards.sample} not found" >> {log}
            echo "Exiting..."
            exit 1
        fi

        if [[ -s {params.normalised_coverage} ]];
        then
            echo "Genome-normalised coverage table created for {wildcards.sample}"
            touch {output.marker}
        else
            echo "Genome-normalised coverage table FAILED for {wildcards.sample}"
            exit 1
        fi
        """
