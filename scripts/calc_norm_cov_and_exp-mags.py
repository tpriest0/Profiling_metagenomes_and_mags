import pandas as pd
import argparse
import numpy as np
import os

#####
# Command-line arguments
######
parser = argparse.ArgumentParser(description="Expand normalized gene representative coverage to all genes and MAGs, then compute expression profiles consistently.")
parser.add_argument("-g", "--metag_cov", required=True, help="Normalized MetaG coverage of gene cluster representatives")
parser.add_argument("-t", "--metat_cov", required=True, help="Normalized MetaT coverage of gene cluster representatives")
parser.add_argument("-n", "--sample_pairs", required=True, help="File mapping MetaG_sample to MetaT_sample")
parser.add_argument("-m", "--mapping_file", required=True, help="Mapping of representative -> member -> contig -> MAG")
parser.add_argument("-s", "--scgs_file", required=True, help="Mapping of representative -> single copy gene COGs")
parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
parser.add_argument("-p", "--out_prefix", required=True, help="Output file prefix")
args = parser.parse_args()

#####
# Load inputs
#####
print("[INFO] Loading input tables...")
metaG_reps = pd.read_csv(args.metag_cov, sep="\t")
metaT_reps = pd.read_csv(args.metat_cov, sep="\t")
sample_pairs = pd.read_csv(args.sample_pairs, sep="\t")
map_df = pd.read_csv(args.mapping_file, sep="\t")
scg_map = pd.read_csv(args.scgs_file, sep="\t")

# Consistency check
print(f"[INFO] MetaG samples: {metaG_reps['Sample'].nunique()} | MetaT samples: {metaT_reps['Sample'].nunique()}")
print(f"[INFO] Loaded mapping file: {len(map_df)} mappings")

#####
# Add MAG mapping for genes
##### 
metaG_reps = metaG_reps.merge(
    map_df, left_on="Gene_name", right_on="Gene_cluster_representative", how="left"
)[["Sample","Gene_name","Coverage_per_cell","MAG_name"]]

metaT_reps = metaT_reps.merge(
    map_df, left_on="Gene_name", right_on="Gene_cluster_representative", how="left"
)[["Sample","Gene_name","Coverage_per_cell","MAG_name"]]

# Genes with no match in map_df (left join) get MAG_name = NaN, not the
# string "Unbinned". Without this fill, NaN keys are silently dropped by
# groupby() below, and the "MAG_name == Unbinned" checks further down
# never match anything -- genuinely unbinned genes would vanish from the
# output entirely rather than being reported as their own category.
metaG_reps["MAG_name"] = metaG_reps["MAG_name"].fillna("Unbinned")
metaT_reps["MAG_name"] = metaT_reps["MAG_name"].fillna("Unbinned")

metaG_reps = metaG_reps.drop_duplicates(subset=["Sample","Gene_name","MAG_name"])
metaT_reps = metaT_reps.drop_duplicates(subset=["Sample","Gene_name","MAG_name"])
    
    
####################
# Calculate metagenome coverage of MAGs
####################

# STEP 1 — Attach SCG type (COG) to coverage table
metaG_scg = metaG_reps.merge(scg_map, on="Gene_name", how="inner")
metaT_scg = metaT_reps.merge(scg_map, on="Gene_name", how="inner")

# STEP 2 — Separate MAG SCGs vs Unbinned SCGs
metaG_scg_mag = metaG_scg[metaG_scg["MAG_name"] != "Unbinned"]
metaG_scg_unb = metaG_scg[metaG_scg["MAG_name"] == "Unbinned"]

metaT_scg_mag = metaT_scg[metaT_scg["MAG_name"] != "Unbinned"]
metaT_scg_unb = metaT_scg[metaT_scg["MAG_name"] == "Unbinned"]

# STEP 3 — MAG-level mean SCG coverage (normal MAGs)

metaG_mag = (
    metaG_scg_mag.groupby(["Sample","MAG_name"], as_index=False)
                 .agg(Mean_coverage_MG=("Coverage_per_cell","mean"))
)

metaT_mag = (
    metaT_scg_mag.groupby(["Sample","MAG_name"], as_index=False)
                 .agg(Mean_coverage_MT=("Coverage_per_cell","mean"))
)


# STEP 4 — UNBINNED abundance (sum per COG then mean across COGs)

# MetaG
unbinned_G = (
    metaG_scg_unb.groupby(["Sample","COG"], as_index=False)
                 .agg(sum_cov=("Coverage_per_cell","sum"))
)

unbinned_G = (
    unbinned_G.groupby("Sample", as_index=False)
              .agg(Mean_coverage_MG=("sum_cov","mean"))
)
unbinned_G["MAG_name"] = "Unbinned"

# MetaT
unbinned_T = (
    metaT_scg_unb.groupby(["Sample","COG"], as_index=False)
                 .agg(sum_cov=("Coverage_per_cell","sum"))
)

unbinned_T = (
    unbinned_T.groupby("Sample", as_index=False)
              .agg(Mean_coverage_MT=("sum_cov","mean"))
)
unbinned_T["MAG_name"] = "Unbinned"


# STEP 5 — Combine MAGs + Unbinned, compute totals, compute proportions

# MetaG
metaG_mag_full = pd.concat([metaG_mag, unbinned_G], ignore_index=True)

metaG_totals = (
    metaG_mag_full.groupby("Sample", as_index=False)
                  .agg(total=("Mean_coverage_MG","sum"))
)

metaG_mag_full = metaG_mag_full.merge(metaG_totals, on="Sample")
metaG_mag_full["Proportion_of_cells_MG"] = (
    metaG_mag_full["Mean_coverage_MG"] / metaG_mag_full["total"]
)
metaG_mag_full = metaG_mag_full.drop(columns=["total"])


# MetaT
metaT_mag_full = pd.concat([metaT_mag, unbinned_T], ignore_index=True)

metaT_totals = (
    metaT_mag_full.groupby("Sample", as_index=False)
                  .agg(total=("Mean_coverage_MT","sum"))
)

metaT_mag_full = metaT_mag_full.merge(metaT_totals, on="Sample")
metaT_mag_full["Proportion_of_cells_MT"] = (
    metaT_mag_full["Mean_coverage_MT"] / metaT_mag_full["total"]
)
metaT_mag_full = metaT_mag_full.drop(columns=["total"])

# Save results
metaG_mag_full.to_csv(
    os.path.join(args.output_dir, f"{args.out_prefix}.mags.MG_cov_normalised.tsv"),
    sep="\t", index=False
)

metaT_mag_full.to_csv(
    os.path.join(args.output_dir, f"{args.out_prefix}.mags.MT_cov_normalised.tsv"),
    sep="\t", index=False
)

#####
# Compute expression of MAGs
#####
expr_mags = pd.merge(
    metaG_mag_full.merge(sample_pairs, left_on="Sample", right_on="MetaG_sample"),
    metaT_mag_full.rename(columns={"Sample":"MetaT_sample"}),
    on=["MAG_name","MetaT_sample"],
    how="inner"
)

mask = (expr_mags["Proportion_of_cells_MG"] > 0) & (expr_mags["Proportion_of_cells_MT"] > 0)
expr_mags["log2_expr"] = np.nan
expr_mags.loc[mask, "log2_expr"] = np.log2(
    expr_mags.loc[mask, "Proportion_of_cells_MT"] /
    expr_mags.loc[mask, "Proportion_of_cells_MG"]
)

# Export
if not expr_mags.empty:
    print(f"[INFO] Expression profile of MAGs generated.")
    expr_mags.to_csv(os.path.join(args.output_dir, f"{args.out_prefix}.mags.expression_profile.tsv"), sep="\t", index=False)
