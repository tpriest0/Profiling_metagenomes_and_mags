#!/bin/bash
#SBATCH --ntasks=1                 # Total number of tasks (processes)
#SBATCH --cpus-per-task=48           # CPUs per task (small number per process)
#SBATCH --time=1-24:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --tmp=2000                        # per node!!
#SBATCH --job-name=fetchmgs
#SBATCH --output=02output/logs/fetchmgs.out
#SBATCH --error=02output/logs/fetchmgs.err
#SBATCH --partition=sunagawa

#####################################

#####
# CONFIGURATION
#####
eval "$(conda shell.bash hook)"
conda activate fetchmgs_env

SCRIPTS_DIR="/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/tpriest/projects/bledina_cadagno/01code/gene_coverage_pipeline"
OUTPUT_DIR="/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/tpriest/projects/bledina_cadagno/02output/scgs"
INPUT_GENES_AA="/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/tpriest/projects/bledina_cadagno/02output/gene_catalog/CadagnoSediment_2024.genes.reps.aa.fa"
INPUT_GENES_NT="/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/tpriest/projects/bledina_cadagno/02output/gene_catalog/CadagnoSediment_2024.genes.reps.nt.fa"

#####################################

#####
# RUN SCGS SEARCH
#####
if [[ -s "$INPUT_GENES_AA" ]]; 
then
    echo "[$(date)] Input file $INPUT_GENES found..."
    echo "[$(date)] Running FetchMGs"

    fetchMGs extraction "$INPUT_GENES_AA" gene "$OUTPUT_DIR" -d "$INPUT_GENES_NT" -t "$SLURM_CPUS_PER_TASK"

else
    echo "[$(date)] ERROR: Input file $INPUT_GENES not found or empty" >&2
    exit 1
fi

if [[ -s "$OUTPUT_DIR/CadagnoSediment_2024.genes.reps.scgs.aa.fa.fetchMGs.faa" ]]; 
then
    echo "[$(date)] FetchMGs run complete"
    echo "[$(date)] Reformatting outputs..."

    grep '>' "$OUTPUT_DIR/CadagnoSediment_2024.genes.reps.scgs.aa.fa.fetchMGs.faa" | \
    sed 's/>//' | sed 's/.COG/\tCOG/' | \
    sed '1i Gene_name\tCOG' > "$OUTPUT_DIR/gene_to_cog_name.tsv" 

    sed -i 's/.COG.*//' "$OUTPUT_DIR/CadagnoSediment_2024.genes.reps.scgs.aa.fa.fetchMGs.faa"
    sed -i 's/.COG.*//' "$OUTPUT_DIR/CadagnoSediment_2024.genes.reps.scgs.aa.fa.fetchMGs.fna"

else
    echo "[$(date)] FetchMGs FAILED"
    echo "[$(date)] Exiting..."
    exit 1
fi

