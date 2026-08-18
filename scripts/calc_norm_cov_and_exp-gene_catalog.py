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
parser.add_argument("-t", "--metat_dir", default=None,
                    help="Directory containing metatranscriptomic coverage files (*.tsv). If omitted, MT/expression is skipped.")
parser.add_argument("-n", "--name_matching", default=None,
                    help="File containing metagenome to metatranscriptome name matching. Required only if MT/expression is computed.")
parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
parser.add_argument("-p", "--output_prefix", required=True, help="Output file prefix")
parser.add_argument("-s", "--scgs_list", required=True, help="Path to SCG list mapping COG to Ref_name")
parser.add_argument("--threads", type=int, default=4, help="Number of threads")
parser.add_argument("--epsilon", type=float, default=1e-6, help="Small value to avoid log(0)")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

#####
# Load single-copy marker gene list
#####
scg_map = pd.read_csv(args.scgs_list, sep="\t")

#####
# Process coverage files: normalize by SCG mean
#####
def process_coverage_file(cov_file, scg_map):
    sample_name = os.path.basename(cov_file).replace(".coverage.tsv", "")
    df = pd.read_csv(cov_file, sep="\t")

    # Require the column we're going to normalize
    if "Median_coverage_per_kbp" not in df.columns:
        raise ValueError(
            f"{cov_file} is missing required column 'Median_coverage_per_kbp'. "
            f"Columns are: {list(df.columns)}"
        )

    # Compute median SCG coverage
    scg_coverages = defaultdict(float)
    for cog_id, group in scg_map.groupby("COG"):
        matched = df[df["Gene_name"].isin(group["Gene_name"])]
        scg_coverages[cog_id] = matched["Median_coverage_per_kbp"].sum() if not matched.empty else np.nan

    scg_median_cov = np.nanmedian(list(scg_coverages.values()))
    
    if np.isnan(scg_median_cov) or scg_median_cov == 0:
        raise ValueError(
            f"SCG median coverage is {scg_median_cov} for sample {sample_name}. "
            f"This usually means SCGs were not found or have zero coverage."
        )
    
    # Normalise all genes by the SCG coverage
    df["Mean_depth_per_genome"] = df["Median_coverage_per_kbp"] / scg_median_cov
    df["Sample"] = sample_name


    keep_cols = [
        "Sample",
        "Gene_name",
        "Gene_length",
        "Horizontal_coverage",
        "Median_coverage_per_kbp",
        "Mean_depth_per_genome",
    ]
    # keep only columns that exist (in case some runs have fewer)
    keep_cols = [c for c in keep_cols if c in df.columns]

    return sample_name, df[keep_cols], scg_median_cov

#####
# Execute normalisation in parallel across coverage tables
#####
def process_dir(cov_dir, scg_map):
    coverage_files = glob.glob(os.path.join(cov_dir, "*.coverage.tsv"))
    results = {}
    cells = {}

    if len(coverage_files) == 0:
        return results, cells

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_coverage_file, f, scg_map): f for f in coverage_files}
        for future in as_completed(futures):
            sample, df, n_cells = future.result()
            results[sample] = df
            cells[sample] = n_cells

    return results, cells


#####
# Apply to metaG and, if available, metaT directories
#####
metaG_data, metaG_cells = process_dir(args.metag_dir, scg_map)
if len(metaG_data) == 0:
    raise RuntimeError(f"[ERROR] No MetaG coverage files found in: {args.metag_dir}")


do_metat = False
metaT_data, metaT_cells = {}, {}

if args.metat_dir:
    metaT_data, metaT_cells = process_dir(args.metat_dir, scg_map)
    if len(metaT_data) > 0:
        do_metat = True
    else:
        print(f"[WARN] --metat_dir provided but no *.coverage.tsv files found in: {args.metat_dir}. Skipping MT/expression.")

#####
# Export tables containing the number of cells sequenced
#####
mg_cells_df = pd.DataFrame({
    "Sample": list(metaG_cells.keys()),
    "Number_of_cells_sequenced": list(metaG_cells.values())
}).sort_values("Sample")

mg_cells_out = os.path.join(
    args.output_dir,
    f"{args.output_prefix}.MG_number_of_cells_sequenced.tsv"
)
mg_cells_df.to_csv(mg_cells_out, sep="\t", index=False)
print(f"[INFO] Wrote MetaG Number_of_cells_sequenced: {mg_cells_out}")

if do_metat:
    mt_cells_df = pd.DataFrame({
        "Sample": list(metaT_cells.keys()),
        "Number_of_cells_sequenced": list(metaT_cells.values())
    }).sort_values("Sample")

    mt_cells_out = os.path.join(
        args.output_dir,
        f"{args.output_prefix}.MT_number_of_cells_sequenced.tsv"
    )
    mt_cells_df.to_csv(mt_cells_out, sep="\t", index=False)
    print(f"[INFO] Wrote MetaT Number_of_cells_sequenced: {mt_cells_out}")

#####
# Export normalised coverage profiles
#####
os.makedirs(args.output_dir, exist_ok=True)

metaG_all = pd.concat(metaG_data.values(), ignore_index=True)
metaG_outfile = os.path.join(args.output_dir, f"{args.output_prefix}.genes.reps.MG_cov_normalised.tsv")
metaG_all.to_csv(metaG_outfile, sep="\t", index=False)
print(f"[INFO] Wrote MetaG normalized coverage: {metaG_outfile} ({len(metaG_all)} rows)")

if do_metat:
    metaT_all = pd.concat(metaT_data.values(), ignore_index=True)
    metaT_outfile = os.path.join(args.output_dir, f"{args.output_prefix}.genes.reps.MT_cov_normalised.tsv")
    metaT_all.to_csv(metaT_outfile, sep="\t", index=False)
    print(f"[INFO] Wrote MetaT normalized coverage: {metaT_outfile} ({len(metaT_all)} rows)")

#####
# For paired samples, compute expression
#####
used_pairs = 0
expr_records = []

if do_metat:
    if not args.name_matching:
        print("[WARN] MetaT data present but --name_matching not provided. Skipping expression calculation.")
    else:
        mapping_df = pd.read_csv(args.name_matching, sep="\t")

        for _, row in mapping_df.iterrows():
            mg_name = row["MetaG_sample"]
            mt_name = row["MetaT_sample"]

            if mg_name in metaG_data and mt_name in metaT_data:
                mg_df = metaG_data[mg_name]
                mt_df = metaT_data[mt_name]

                merged = pd.merge(
                    mg_df[["Gene_name", "Mean_depth_per_genome"]],
                    mt_df[["Gene_name", "Mean_depth_per_genome"]],
                    on="Gene_name",
                    suffixes=("_MG", "_MT")
                )

                merged["MetaG_sample"] = mg_name
                merged["MetaT_sample"] = mt_name
                merged["log2_expr"] = np.log2(
                    (merged["Mean_depth_per_genome_MT"] + args.epsilon) /
                    (merged["Mean_depth_per_genome_MG"] + args.epsilon)
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
        else:
            print("[WARN] No expression records generated (no valid MG↔MT pairs found).")

print(f"[INFO] Processed {len(metaG_data)} MetaG samples")
print(f"[INFO] Processed {len(metaT_data)} MetaT samples" if do_metat else "[INFO] MetaT processing skipped")
print(f"[INFO] Mapped pairs used for expression: {used_pairs}" if do_metat and args.name_matching else "[INFO] Expression skipped")
