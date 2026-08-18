import pandas as pd
import os
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import glob

#####
# Command-line arguments
#####
parser = argparse.ArgumentParser(description="Normalize coverage by SCG mean and compute expression for shared samples.")
parser.add_argument("-g", "--metag_dir", required=True, help="Directory containing metagenomic coverage files (*.tsv)")
parser.add_argument("-t", "--metat_dir", required=True, help="Directory containing metatranscriptomic coverage files (*.tsv)")
parser.add_argument("-n", "--name_matching", required=True, help="File containing metagenome to metatranscriptome name matching")
parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
parser.add_argument("-p", "--output_prefix", required=True, help="Output file prefix")
parser.add_argument("-s", "--scgs_list", required=True, help="Path to SCG list mapping Gene_name to COG")
parser.add_argument("-b", "--breadth_threshold", type=float, default=95.0,
                     help="Minimum horizontal coverage (Prop_bases_covered, on a 0-100 scale) required to retain a gene (default: 95.0)")
parser.add_argument("--threads", type=int, default=4, help="Number of threads")
args = parser.parse_args()

#####
# Load single-copy marker gene list
#####
scg_map = pd.read_csv(args.scgs_list, sep="\t")
required_scg_cols = {"Gene_name", "COG"}
missing_scg_cols = required_scg_cols - set(scg_map.columns)
if missing_scg_cols:
    raise ValueError(f"SCG list ({args.scgs_list}) is missing required column(s): {missing_scg_cols}")

#####
# Process coverage files: filter by breadth, normalize by SCG mean
#####
def process_coverage_file(cov_file, scg_map, breadth_threshold):
    sample_name = os.path.basename(cov_file).replace(".coverage.tsv", "")
    df = pd.read_csv(cov_file, sep="\t")
    df = df.rename(columns={"Gene": "Gene_name"})

    required_cols = {"Gene_name", "Gene_length", "Num_bases_covered", "Prop_bases_covered",
                      "Mean_depth", "Mean_depth_per_kbp"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"{cov_file} is missing required column(s): {missing_cols}")

    # Task 1: apply horizontal coverage filter (threshold is user-configurable)
    df = df[df["Prop_bases_covered"] >= breadth_threshold].copy()

    if df.empty:
        print(f"[WARN] {sample_name}: no genes passed the breadth threshold ({breadth_threshold}) -- skipping")
        return sample_name, df.assign(Sample=sample_name, Coverage_per_cell=np.nan)

    # Task 3: mean SCG coverage -> Number_of_genomes_sequenced proxy.
    # Per-COG sums use the length-normalised Mean_depth_per_kbp column;
    # COGs with no genes passing the filter contribute NaN and are
    # excluded from the mean (np.nanmean skips them) -- i.e. an
    # observed-only denominator, matching the workflow's normalization rule.
    scg_coverages = defaultdict(float)
    for cog_id, group in scg_map.groupby("COG"):
        matched = df[df["Gene_name"].isin(group["Gene_name"])]
        scg_coverages[cog_id] = matched["Mean_depth_per_kbp"].sum() if not matched.empty else np.nan

    scg_vals = np.array(list(scg_coverages.values()), dtype=float)

    if np.all(np.isnan(scg_vals)):
        print(f"[WARN] {sample_name}: no single-copy marker genes passed the breadth threshold -- "
              f"Number_of_genomes_sequenced cannot be estimated (Coverage_per_cell will be NaN)")

    # Mean SCG coverage used for normalization (Number_of_genomes_sequenced proxy)
    scg_mean_cov = np.nanmean(scg_vals)

    # Task 4: normalize every retained gene's length-normalised coverage
    # by the genome-equivalent estimate
    df["Coverage_per_cell"] = df["Mean_depth_per_kbp"] / scg_mean_cov
    df["Sample"] = sample_name

    # Include extra info for downstream use
    return sample_name, df[["Sample", "Gene_name", "Gene_length", "Num_bases_covered",
                             "Prop_bases_covered", "Mean_depth", "Mean_depth_per_kbp", "Coverage_per_cell"]]

#####
# Execute normalization in parallel across coverage tables
#####
def process_dir(cov_dir, scg_map, breadth_threshold):
    coverage_files = glob.glob(os.path.join(cov_dir, "*.coverage.tsv"))
    results = {}
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_coverage_file, f, scg_map, breadth_threshold): f for f in coverage_files}
        for future in as_completed(futures):
            sample, df = future.result()
            results[sample] = df
    return results

#####
# Apply to both metaG and metaT directories
#####
metaG_data = process_dir(args.metag_dir, scg_map, args.breadth_threshold)
metaT_data = process_dir(args.metat_dir, scg_map, args.breadth_threshold)

#####
# Export normalized coverage profiles
#####
metaG_all = pd.concat(metaG_data.values(), ignore_index=True)
metaT_all = pd.concat(metaT_data.values(), ignore_index=True)

os.makedirs(args.output_dir, exist_ok=True)

metaG_outfile = os.path.join(args.output_dir, f"{args.output_prefix}.genes.reps.MG_cov_normalised.tsv")
metaT_outfile = os.path.join(args.output_dir, f"{args.output_prefix}.genes.reps.MT_cov_normalised.tsv")

metaG_all.to_csv(metaG_outfile, sep="\t", index=False)
metaT_all.to_csv(metaT_outfile, sep="\t", index=False)

print(f"[INFO] Wrote MetaG normalized coverage: {metaG_outfile} ({len(metaG_all)} rows)")
print(f"[INFO] Wrote MetaT normalized coverage: {metaT_outfile} ({len(metaT_all)} rows)")

#####
# For paired samples, compute expression
#####
mapping_df = pd.read_csv(args.name_matching, sep="\t")

expr_records = []
used_pairs = 0
for _, row in mapping_df.iterrows():
    mg_name = row["MetaG_sample"]
    mt_name = row["MetaT_sample"]

    if mg_name in metaG_data and mt_name in metaT_data:
        mg_df = metaG_data[mg_name]
        mt_df = metaT_data[mt_name]

        merged = pd.merge(
            mg_df[["Gene_name", "Coverage_per_cell"]],
            mt_df[["Gene_name", "Coverage_per_cell"]],
            on="Gene_name",
            suffixes=("_MG", "_MT")
        )

        merged["MetaG_sample"] = mg_name
        merged["MetaT_sample"] = mt_name
        merged["log2_expr"] = np.nan
        mask = (merged["Coverage_per_cell_MG"] > 0) & (merged["Coverage_per_cell_MT"] > 0)
        merged.loc[mask, "log2_expr"] = np.log2(
            merged.loc[mask, "Coverage_per_cell_MT"] / merged.loc[mask, "Coverage_per_cell_MG"]
        )

        expr_records.append(merged)
        used_pairs += 1
    else:
        print(f"[WARN] Skipping pair {mg_name} ↔ {mt_name}, missing coverage data")

if expr_records:
    expr_df = pd.concat(expr_records, ignore_index=True)
    expr_outfile = os.path.join(args.output_dir, f"{args.output_prefix}.genes.reps.expression_profile.tsv")
    expr_df.to_csv(expr_outfile, sep="\t", index=False)
    print(f"[INFO] Wrote expression profile: {expr_outfile} ({len(expr_df)} rows)")

print(f"[INFO] Processed {len(metaG_data)} MetaG samples")
print(f"[INFO] Processed {len(metaT_data)} MetaT samples")
print(f"[INFO] Mapped pairs used for expression: {used_pairs}")
