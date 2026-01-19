import pandas as pd
import argparse
import numpy as np
import os

#####
# Command-line arguments
#####
parser = argparse.ArgumentParser(
    description="Compute MAG-level normalized coverage and expression using SCG median logic (sum per COG, median across COGs). "
                "If MetaT inputs are not provided, MT/expression outputs are skipped."
)
parser.add_argument("-g", "--metag_cov", required=True, help="Normalized MetaG coverage of gene cluster representatives")
parser.add_argument("-t", "--metat_cov", default=None,
                    help="Normalized MetaT coverage of gene cluster representatives. If omitted, MT/expression is skipped.")
parser.add_argument("-n", "--sample_pairs", default=None,
                    help="File mapping MetaG_sample to MetaT_sample. Required only if expression is computed.")
parser.add_argument("-m", "--mapping_file", required=True, help="Mapping of representative -> member -> contig -> MAG")
parser.add_argument("-s", "--scgs_file", required=True, help="Mapping of representative -> single copy gene COGs")
parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
parser.add_argument("-p", "--out_prefix", required=True, help="Output file prefix")
parser.add_argument("--epsilon", type=float, default=1e-6, help="Small value to avoid log(0)")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

#####
# Load inputs and run checks
#####
print("[INFO] Loading input tables...")
metaG_reps = pd.read_csv(args.metag_cov, sep="\t")
map_df = pd.read_csv(args.mapping_file, sep="\t")
scg_map = pd.read_csv(args.scgs_file, sep="\t")

# Optional MetaT/sample pairs
do_metat = False
metaT_reps = None
sample_pairs = None

if args.metat_cov:
    metaT_reps = pd.read_csv(args.metat_cov, sep="\t")
    do_metat = True

if do_metat and args.sample_pairs:
    sample_pairs = pd.read_csv(args.sample_pairs, sep="\t")
elif do_metat and not args.sample_pairs:
    print("[WARN] MetaT coverage provided but --sample_pairs not provided. Expression will be skipped.")

# Basic checks
required_cov_cols = {"Sample", "Gene_name", "Coverage_per_cell"}
if not required_cov_cols.issubset(metaG_reps.columns):
    raise ValueError(f"[ERROR] MetaG coverage file missing columns {required_cov_cols}. Found: {set(metaG_reps.columns)}")

if do_metat and (not required_cov_cols.issubset(metaT_reps.columns)):
    raise ValueError(f"[ERROR] MetaT coverage file missing columns {required_cov_cols}. Found: {set(metaT_reps.columns)}")

if not {"Gene_cluster_representative", "MAG_name"}.issubset(map_df.columns):
    raise ValueError(
        "[ERROR] mapping_file must contain columns: Gene_cluster_representative, MAG_name. "
        f"Found: {list(map_df.columns)}"
    )

if do_metat and (sample_pairs is not None):
    if not {"MetaG_sample", "MetaT_sample"}.issubset(sample_pairs.columns):
        raise ValueError(
            "[ERROR] sample_pairs file must contain columns: MetaG_sample, MetaT_sample. "
            f"Found: {list(sample_pairs.columns)}"
        )

print(f"[INFO] MetaG samples: {metaG_reps['Sample'].nunique()}")
if do_metat:
    print(f"[INFO] MetaT samples: {metaT_reps['Sample'].nunique()}")
else:
    print("[INFO] MetaT not provided: skipping MT/expression.")

print(f"[INFO] Loaded mapping file: {len(map_df)} mappings")

#####
# Attach MAG mapping for representative genes
#####
metaG_reps = metaG_reps.merge(
    map_df, left_on="Gene_name", right_on="Gene_cluster_representative", how="left"
)[["Sample", "Gene_name", "Coverage_per_cell", "MAG_name"]]

metaG_reps = metaG_reps.drop_duplicates(subset=["Sample", "Gene_name", "MAG_name"])

if do_metat:
    metaT_reps = metaT_reps.merge(
        map_df, left_on="Gene_name", right_on="Gene_cluster_representative", how="left"
    )[["Sample", "Gene_name", "Coverage_per_cell", "MAG_name"]]

    metaT_reps = metaT_reps.drop_duplicates(subset=["Sample", "Gene_name", "MAG_name"])

#####
# Attach SCG COG to representative gene table
#####
metaG_scg = metaG_reps.merge(scg_map, on="Gene_name", how="inner")
if do_metat:
    metaT_scg = metaT_reps.merge(scg_map, on="Gene_name", how="inner")
    
#####
# Helper: MAG abundance from SCGs using "sum per COG then median across COGs"
#####
def mag_abundance_from_scgs(df_scg, cov_col, out_col):
    """
    df_scg must contain: Sample, MAG_name, COG, cov_col
    Returns: Sample, MAG_name, out_col
    """
    # Sum SCG coverage within each COG (handles multiple genes per COG)
    per_cog = (
        df_scg.groupby(["Sample", "MAG_name", "COG"], as_index=False)
              .agg(cog_sum=(cov_col, "sum"))
    )

    # Median across COGs -> robust MAG estimate consistent with your per-sample anchor logic
    mag = (
        per_cog.groupby(["Sample", "MAG_name"], as_index=False)
               .agg(**{out_col: ("cog_sum", "median")})
    )
    return mag

#####
# Calculate MAG coverage (MetaG / MetaT)
#####

### MetaG

# Split MAG vs Unbinned (MetaG)
metaG_scg_mag = metaG_scg[metaG_scg["MAG_name"] != "Unbinned"]
metaG_scg_unb = metaG_scg[metaG_scg["MAG_name"] == "Unbinned"]

# MAGs (binned): sum per COG then median across COGs
metaG_mag = mag_abundance_from_scgs(metaG_scg_mag, "Coverage_per_cell", "MAG_cov_MG")

# Unbinned: sum per COG then median across COGs
unbinned_G = (
    metaG_scg_unb.groupby(["Sample", "COG"], as_index=False)
                 .agg(cog_sum=("Coverage_per_cell", "sum"))
)
unbinned_G = (
    unbinned_G.groupby("Sample", as_index=False)
              .agg(MAG_cov_MG=("cog_sum", "median"))
)
unbinned_G["MAG_name"] = "Unbinned"

# Combine MAGs + Unbinned
metaG_mag_full = pd.concat([metaG_mag, unbinned_G], ignore_index=True)

# Totals and proportions per sample
metaG_totals = metaG_mag_full.groupby("Sample", as_index=False).agg(total=("MAG_cov_MG", "sum"))
metaG_mag_full = metaG_mag_full.merge(metaG_totals, on="Sample", how="left")
metaG_mag_full["Proportion_of_cells_MG"] = metaG_mag_full["MAG_cov_MG"] / metaG_mag_full["total"]
metaG_mag_full = metaG_mag_full.drop(columns=["total"])

# Save MG results
metaG_mag_full.to_csv(
    os.path.join(args.output_dir, f"{args.out_prefix}.mags.MG_cov_normalised.tsv"),
    sep="\t", index=False
)
print("[INFO] Wrote MetaG MAG normalized coverage table.")

# Save MG results
metaG_mag_full.to_csv(
    os.path.join(args.output_dir, f"{args.out_prefix}.mags.MG_cov_normalised.tsv"),
    sep="\t", index=False
)
print("[INFO] Wrote MetaG MAG normalized coverage table.")

### MetaT, if provided
metaT_mag_full = None
if do_metat:
    metaT_scg_mag = metaT_scg[metaT_scg["MAG_name"] != "Unbinned"]
    metaT_scg_unb = metaT_scg[metaT_scg["MAG_name"] == "Unbinned"]

    metaT_mag = mag_abundance_from_scgs(metaT_scg_mag, "Coverage_per_cell", "MAG_cov_MT")

    unbinned_T = (
        metaT_scg_unb.groupby(["Sample", "COG"], as_index=False)
                     .agg(cog_sum=("Coverage_per_cell", "sum"))
    )
    unbinned_T = (
        unbinned_T.groupby("Sample", as_index=False)
                  .agg(MAG_cov_MT=("cog_sum", "median"))
    )
    unbinned_T["MAG_name"] = "Unbinned"

    metaT_mag_full = pd.concat([metaT_mag, unbinned_T], ignore_index=True)

    metaT_totals = metaT_mag_full.groupby("Sample", as_index=False).agg(total=("MAG_cov_MT", "sum"))
    metaT_mag_full = metaT_mag_full.merge(metaT_totals, on="Sample", how="left")
    metaT_mag_full["Proportion_of_cells_MT"] = metaT_mag_full["MAG_cov_MT"] / metaT_mag_full["total"]
    metaT_mag_full = metaT_mag_full.drop(columns=["total"])

    metaT_mag_full.to_csv(
        os.path.join(args.output_dir, f"{args.out_prefix}.mags.MT_cov_normalised.tsv"),
        sep="\t", index=False
    )
    print("[INFO] Wrote MetaT MAG normalized coverage table.")

#####
# Compute expression of MAGs (optional)
#####
if do_metat and (sample_pairs is not None):
    expr_mags = pd.merge(
        metaG_mag_full.merge(sample_pairs, left_on="Sample", right_on="MetaG_sample"),
        metaT_mag_full.rename(columns={"Sample": "MetaT_sample"}),
        on=["MAG_name", "MetaT_sample"],
        how="inner"
    )

    mask = (expr_mags["Proportion_of_cells_MG"] > 0) & (expr_mags["Proportion_of_cells_MT"] > 0)
    expr_mags["log2_expr"] = np.nan
    expr_mags.loc[mask, "log2_expr"] = np.log2(
        (expr_mags.loc[mask, "Proportion_of_cells_MT"] + args.epsilon) /
        (expr_mags.loc[mask, "Proportion_of_cells_MG"] + args.epsilon)
    )

    if not expr_mags.empty:
        expr_out = os.path.join(args.output_dir, f"{args.out_prefix}.mags.expression_profile.tsv")
        expr_mags.to_csv(expr_out, sep="\t", index=False)
        print(f"[INFO] Wrote MAG expression profile: {expr_out}")
    else:
        print("[WARN] No MAG expression pairs produced (check sample_pairs and MAG overlap).")
else:
    if do_metat and (sample_pairs is None):
        print("[INFO] Skipping MAG expression (missing --sample_pairs).")
    elif not do_metat:
        print("[INFO] Skipping MAG expression (no MetaT coverage provided).")