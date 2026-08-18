#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="Generate KEGG/PFAM functional coverage & expression profiles (genes + MAGs).")
parser.add_argument("-g", "--metag_cov", required=True, help="Normalized MetaG coverage of gene cluster representatives.")
parser.add_argument("-t", "--metat_cov", required=True, help="Normalized MetaT coverage of gene cluster representatives.")
parser.add_argument("-s", "--sample_pairs", required=True, help="MetaG↔MetaT sample mapping file with columns: MetaG_sample, MetaT_sample")
parser.add_argument("-a", "--annotations", required=True, help="Table mapping Gene_name → KEGG_ko (at least these two columns)")
parser.add_argument("-m", "--mapping_file", required=True, help="Mapping of gene cluster representative → MAG_name")
parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
parser.add_argument("-p", "--out_prefix", required=True, help="Output file prefix")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

#####
# Load inputs
#####
print("[INFO] Loading input tables...")
metaG_reps = pd.read_csv(args.metag_cov, sep="\t")
metaT_reps = pd.read_csv(args.metat_cov, sep="\t")
sample_pairs = pd.read_csv(args.sample_pairs,    sep="\t")
annotations  = pd.read_csv(args.annotations, sep="\t")
map_df = pd.read_csv(args.mapping_file, sep="\t")

# Ensure required columns exist
for col in ["Sample","Gene_name","Coverage_per_cell"]:
    if col not in metaG_reps.columns: raise ValueError(f"MetaG file missing column: {col}")
    if col not in metaT_reps.columns: raise ValueError(f"MetaT file missing column: {col}")
for col in ["MetaG_sample","MetaT_sample"]:
    if col not in sample_pairs.columns: raise ValueError(f"sample_pairs file missing column: {col}")

#####
# Add MAG mapping for genes
##### 
metaG_reps = metaG_reps.merge(
    map_df, left_on="Gene_name", right_on="Gene_cluster_representative", how="left"
)[["Sample","Gene_name","Coverage_per_cell","MAG_name"]]

metaT_reps = metaT_reps.merge(
    map_df, left_on="Gene_name", right_on="Gene_cluster_representative", how="left"
)[["Sample","Gene_name","Coverage_per_cell","MAG_name"]]

# See calc_norm_cov_and_exp-mags.py for why this fillna is required: a
# left-join leaves unmatched (genuinely unbinned) genes with MAG_name =
# NaN, which groupby() silently drops rather than reporting as their own
# category -- this would otherwise cause MG_prop_of_sample_function /
# MT_prop_of_sample_function to not sum to 1.0 across MAGs, with no
# indication why.
metaG_reps["MAG_name"] = metaG_reps["MAG_name"].fillna("Unbinned")
metaT_reps["MAG_name"] = metaT_reps["MAG_name"].fillna("Unbinned")

metaG_reps = metaG_reps.drop_duplicates(subset=["Sample","Gene_name","MAG_name"])
metaT_reps = metaT_reps.drop_duplicates(subset=["Sample","Gene_name","MAG_name"])

#####
# Process annotation table
#####
kegg = annotations[["Gene_name","KEGG_ko"]].dropna()
pfam = annotations[["Gene_name","PFAM_accession"]].dropna()
metabolic = annotations[["Gene_name","METABOLIC_hmm"]].dropna()

#####
# Calculate coverages and expression of functions at the sample and MAG level
#####

def compute_profiles(func_map: pd.DataFrame, func_col: str, label: str):
    # Attach function to gene rows
    metaG_f = metaG_reps.merge(func_map, on="Gene_name", how="inner")
    metaT_f = metaT_reps.merge(func_map, on="Gene_name", how="inner")

    ### Community function coverage profile: coverage of functions in metaG and metaT (sum of per cell gene coverages per function×sample)
    
    funcG = (metaG_f.groupby(["Sample", func_col], as_index=False)
             .agg(MG_coverage_per_cell=("Coverage_per_cell","sum")))
    funcT = (metaT_f.groupby(["Sample", func_col], as_index=False)
             .agg(MT_coverage_per_cell=("Coverage_per_cell","sum")))

    funcG.to_csv(os.path.join(args.output_dir, f"{args.out_prefix}.genes_reps.{label}.MG_cov_normalised.tsv"), sep="\t", index=False)
    funcT.to_csv(os.path.join(args.output_dir, f"{args.out_prefix}.genes_reps.{label}.MT_cov_normalised.tsv"), sep="\t", index=False)

    ### Community function expression profile: expression of functions (for samples with both MetaG and MetaT)
    
    expr = (funcG.merge(sample_pairs, left_on="Sample", right_on="MetaG_sample", how="inner")
                 .merge(funcT.rename(columns={"Sample":"MetaT_sample"}),
                        on=[func_col,"MetaT_sample"], how="inner"))
    # Compute log2 only when both > 0
    expr["log2_expr"] = np.where(
        (expr["MG_coverage_per_cell"] > 0) & (expr["MT_coverage_per_cell"] > 0),
        np.log2(expr["MT_coverage_per_cell"] / expr["MG_coverage_per_cell"]),
        np.nan
    )
    expr = expr[[func_col,"MetaG_sample","MetaT_sample","MG_coverage_per_cell","MT_coverage_per_cell","log2_expr"]]
    expr.to_csv(os.path.join(args.output_dir, f"{args.out_prefix}.genes_reps.{label}.expression_profile.tsv"), sep="\t", index=False)

    ### MAG function coverage profile: coverage of MAG functions in metaG and metaT in relation to whole community
    
    # Sum gene coverages within each MAG × sample × function
    magG = (metaG_f.groupby(["MAG_name","Sample",func_col], as_index=False)
            .agg(MG_coverage_per_cell=("Coverage_per_cell","sum")))
    magT = (metaT_f.groupby(["MAG_name","Sample",func_col], as_index=False)
            .agg(MT_coverage_per_cell=("Coverage_per_cell","sum")))

    # For proportions, we need sample-level totals per function (across MAGs).
    # These are exactly funcG/funcT; join to compute per-MAG share.
    magG = magG.merge(funcG.rename(columns={"MG_coverage_per_cell":"MG_coverage_per_cell_total"}),
                      on=["Sample",func_col], how="left")
    magG["MG_prop_of_sample_function"] = np.where(
        magG["MG_coverage_per_cell_total"] > 0,
        magG["MG_coverage_per_cell"] / magG["MG_coverage_per_cell_total"],
        np.nan
    )

    magT = magT.merge(funcT.rename(columns={"MT_coverage_per_cell":"MT_coverage_per_cell_total"}),
                      on=["Sample",func_col], how="left")
    magT["MT_prop_of_sample_function"] = np.where(
        magT["MT_coverage_per_cell_total"] > 0,
        magT["MT_coverage_per_cell"] / magT["MT_coverage_per_cell_total"],
        np.nan
    )

    magG_out = magG[["MAG_name","Sample",func_col,"MG_coverage_per_cell","MG_prop_of_sample_function"]]
    magT_out = magT[["MAG_name","Sample",func_col,"MT_coverage_per_cell","MT_prop_of_sample_function"]]

    magG_out.to_csv(os.path.join(args.output_dir, f"{args.out_prefix}.mags.{label}.MG_cov_normalised.tsv"), sep="\t", index=False)
    magT_out.to_csv(os.path.join(args.output_dir, f"{args.out_prefix}.mags.{label}.MT_cov_normalised.tsv"), sep="\t", index=False)

    ### MAG function expression profile: coverage of MAG functions in metaG and metaT in relation to whole community
    
    # Pair via MetaG_sample→MetaT_sample; then bring MAG-level MT function coverage from magT_out
    expr_mag = (magG_out.merge(sample_pairs, left_on="Sample", right_on="MetaG_sample", how="inner")
                       .merge(magT_out.rename(columns={"Sample":"MetaT_sample"}),
                              on=["MAG_name",func_col,"MetaT_sample"], how="inner"))

    # Compute log2 only when both > 0
    expr_mag["log2_expr"] = np.where(
        (expr_mag["MG_coverage_per_cell"] > 0) & (expr_mag["MT_coverage_per_cell"] > 0),
        np.log2(expr_mag["MT_coverage_per_cell"] / expr_mag["MG_coverage_per_cell"]),
        np.nan
    )

    # Keep both MG and MT proportion columns (each relative to its own sample-level function total)
    expr_mag = expr_mag[[
        "MAG_name", func_col, "MetaG_sample", "MetaT_sample",
        "MG_coverage_per_cell", "MG_prop_of_sample_function",
        "MT_coverage_per_cell", "MT_prop_of_sample_function",
        "log2_expr"
    ]]
    expr_mag.to_csv(os.path.join(args.output_dir, f"{args.out_prefix}.mags.{label}.expression_profile.tsv"), sep="\t", index=False)

    print(f"[INFO] {label}: functions={funcG[func_col].nunique()}, MAGs={magG_out['MAG_name'].nunique()}")

# --------- Run for KEGG & PFAM ----------
compute_profiles(kegg, "KEGG_ko", "KEGG")
compute_profiles(pfam, "PFAM_accession", "PFAM")
# compute_profiles(metabolic, "METABOLIC_hmm", "METABOLIC")

print("\n[INFO] ✅ Functional profiling completed.")
