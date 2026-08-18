import pandas as pd
import argparse
import numpy as np
import os

#####
# Command-line arguments
######
parser = argparse.ArgumentParser(
    description="Compute community-level functional coverage (MG always; MT optional) and expression (optional). "
                "Optionally compute MAG-level functional coverage/expression if --mapping_file is provided "
                "(reps duplicated across all linked MAGs)."
)
parser.add_argument("-g", "--metag_cov", required=True,
                    help="Normalised MetaG coverage of gene catalogue representatives (must include Gene_name, Sample, Mean_depth_per_genome).")
parser.add_argument("-t", "--metat_cov", default=None,
                    help="Normalised MetaT coverage of gene catalogue representatives. If omitted, MT/expression is skipped.")
parser.add_argument("-n", "--sample_pairs", default=None,
                    help="File mapping MetaG_sample to MetaT_sample. Required only if expression is computed.")
parser.add_argument("-a", "--annotations", required=True, help="Annotations of gene cluster representatives.")
parser.add_argument("-m", "--mapping_file", default=None,
                    help="Optional mapping of representative -> member -> MAG. If omitted, MAG-level outputs are skipped.")

parser.add_argument("-s", "--scgs_file", default=None,
                    help="Mapping of representative -> single copy gene COGs (optional; loaded for compatibility).")

parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
parser.add_argument("-p", "--output_prefix", required=True, help="Output file prefix")
parser.add_argument("--epsilon", type=float, default=1e-6, help="Small value to avoid log(0) / division-by-zero")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

#####
# Load inputs
#####
print("[INFO] Loading input tables...")
metaG_reps = pd.read_csv(args.metag_cov, sep="\t")

do_metat = False
metaT_reps = None
sample_pairs = None

if args.metat_cov:
    metaT_reps = pd.read_csv(args.metat_cov, sep="\t")
    do_metat = True

if do_metat and args.sample_pairs:
    sample_pairs = pd.read_csv(args.sample_pairs, sep="\t")
elif do_metat and not args.sample_pairs:
    print("[WARN] MetaT coverage provided but --sample_pairs not provided. Expression outputs will be skipped.")

do_mag = False
map_df = None
if args.mapping_file:
    map_df = pd.read_csv(args.mapping_file, sep="\t")
    do_mag = True

# scgs_file is not used in this script, but keep compatibility if users pass it
if args.scgs_file:
    _ = pd.read_csv(args.scgs_file, sep="\t")

# Column checks
required_cov_cols = {"Gene_name", "Sample", "Mean_depth_per_genome"}
if not required_cov_cols.issubset(metaG_reps.columns):
    raise ValueError(f"[ERROR] MetaG coverage file missing columns {required_cov_cols}. Found: {set(metaG_reps.columns)}")

if do_metat and (not required_cov_cols.issubset(metaT_reps.columns)):
    raise ValueError(f"[ERROR] MetaT coverage file missing columns {required_cov_cols}. Found: {set(metaT_reps.columns)}")

if do_metat and (sample_pairs is not None):
    if not {"MetaG_sample", "MetaT_sample"}.issubset(sample_pairs.columns):
        raise ValueError(
            "[ERROR] sample_pairs file must contain columns: MetaG_sample, MetaT_sample. "
            f"Found: {list(sample_pairs.columns)}"
        )

if do_mag:
    if not {"Gene_cluster_representative", "MAG_name"}.issubset(map_df.columns):
        raise ValueError(
            "[ERROR] mapping_file must contain columns: Gene_cluster_representative, MAG_name. "
            f"Found: {list(map_df.columns)}"
        )

###############################################################
# FUNCTIONAL ANALYSIS
###############################################################

print("[INFO] Starting functional profiling...")

# Load annotation table
annotation_df = pd.read_csv(args.annotations, sep="\t")

if not {"Gene_name", "Annotation"}.issubset(annotation_df.columns):
    raise ValueError(
        "[ERROR] annotations file must contain columns: Gene_name, Annotation. "
        f"Found: {list(annotation_df.columns)}"
    )

annotation_df = annotation_df[["Gene_name", "Annotation"]].dropna()
annotation_df = annotation_df[annotation_df["Annotation"] != "Unknown"]

############ COMMUNITY-LEVEL FUNCTION COVERAGE ############
# Use reps only (no MAG duplication) for community coverage.

metaG_comm = metaG_reps.merge(annotation_df, on="Gene_name", how="inner")
func_MG = (
    metaG_comm.groupby(["Sample", "Annotation"], as_index=False)
              .agg(MG_coverage_per_cell=("Mean_depth_per_genome", "sum"))
)

func_MG.to_csv(
    os.path.join(args.output_dir, f"{args.output_prefix}.genes.functions.MG_cov_normalised.tsv"),
    sep="\t", index=False
)

func_MT = None
if do_metat:
    metaT_comm = metaT_reps.merge(annotation_df, on="Gene_name", how="inner")
    func_MT = (
        metaT_comm.groupby(["Sample", "Annotation"], as_index=False)
                  .agg(MT_coverage_per_cell=("Mean_depth_per_genome", "sum"))
    )
    func_MT.to_csv(
        os.path.join(args.output_dir, f"{args.output_prefix}.genes.functions.MT_cov_normalised.tsv"),
        sep="\t", index=False
    )

############ COMMUNITY-LEVEL EXPRESSION (paired samples) ############
if do_metat and (sample_pairs is not None):
    func_MG_key = func_MG.merge(
        sample_pairs,
        left_on="Sample",
        right_on="MetaG_sample",
        how="inner"
    )

    func_expr = func_MG_key.merge(
        func_MT.rename(columns={"Sample": "MetaT_sample"}),
        on=["Annotation", "MetaT_sample"],
        how="inner"
    )

    func_expr["log2_expr"] = np.log2(
        (func_expr["MT_coverage_per_cell"] + args.epsilon) /
        (func_expr["MG_coverage_per_cell"] + args.epsilon)
    )

    func_expr = func_expr[[
        "Annotation", "MetaG_sample", "MetaT_sample",
        "MG_coverage_per_cell", "MT_coverage_per_cell", "log2_expr"
    ]]

    func_expr.to_csv(
        os.path.join(args.output_dir, f"{args.output_prefix}.genes.functions.expression_profile.tsv"),
        sep="\t", index=False
    )
else:
    if do_metat and (sample_pairs is None):
        print("[INFO] Skipping community-level expression (missing --sample_pairs).")
    elif not do_metat:
        print("[INFO] Skipping community-level MT + expression (no MetaT coverage provided).")

###############################################################
# OPTIONAL: MAG-LEVEL FUNCTION COVERAGE + EXPRESSION
###############################################################
if do_mag:
    print("[INFO] MAG mapping provided: computing MAG-level functional profiles...")

    rep_mag_pairs = (
        map_df[["Gene_cluster_representative", "MAG_name"]]
        .dropna()
        .drop_duplicates()
        .rename(columns={"Gene_cluster_representative": "Gene_name"})
    )

    # Duplicate reps across all linked MAGs
    metaG_func = (
        metaG_reps.merge(annotation_df, on="Gene_name", how="inner")
                  .merge(rep_mag_pairs, on="Gene_name", how="inner")
    )

    mag_func_MG = (
        metaG_func.groupby(["MAG_name", "Sample", "Annotation"], as_index=False)
                  .agg(MG_coverage_per_cell=("Mean_depth_per_genome", "sum"))
    )

    # Proportions (Option 1; same methodology as before)
    mg_totals = (
        mag_func_MG.groupby(["Sample", "Annotation"], as_index=False)
                   .agg(total_func_MG=("MG_coverage_per_cell", "sum"))
    )
    mag_func_MG = mag_func_MG.merge(mg_totals, on=["Sample", "Annotation"], how="left")
    mag_func_MG["Proportion_of_sample_coverage_MG"] = (
        mag_func_MG["MG_coverage_per_cell"] / (mag_func_MG["total_func_MG"] + args.epsilon)
    )

    mag_func_MG[[
        "MAG_name", "Sample", "Annotation", "MG_coverage_per_cell",
        "Proportion_of_sample_coverage_MG"
    ]].to_csv(
        os.path.join(args.output_dir, f"{args.output_prefix}.mags.functions.MG_cov_normalised.tsv"),
        sep="\t", index=False
    )

    mag_func_MT = None
    if do_metat:
        metaT_func = (
            metaT_reps.merge(annotation_df, on="Gene_name", how="inner")
                      .merge(rep_mag_pairs, on="Gene_name", how="inner")
        )

        mag_func_MT = (
            metaT_func.groupby(["MAG_name", "Sample", "Annotation"], as_index=False)
                      .agg(MT_coverage_per_cell=("Mean_depth_per_genome", "sum"))
        )

        mt_totals = (
            mag_func_MT.groupby(["Sample", "Annotation"], as_index=False)
                       .agg(total_func_MT=("MT_coverage_per_cell", "sum"))
        )
        mag_func_MT = mag_func_MT.merge(mt_totals, on=["Sample", "Annotation"], how="left")
        mag_func_MT["Proportion_of_sample_coverage_MT"] = (
            mag_func_MT["MT_coverage_per_cell"] / (mag_func_MT["total_func_MT"] + args.epsilon)
        )

        mag_func_MT[[
            "MAG_name", "Sample", "Annotation", "MT_coverage_per_cell",
            "Proportion_of_sample_coverage_MT"
        ]].to_csv(
            os.path.join(args.output_dir, f"{args.output_prefix}.mags.functions.MT_cov_normalised.tsv"),
            sep="\t", index=False
        )

    # MAG FUNCTION EXPRESSION (paired samples)
    if do_metat and (sample_pairs is not None):
        mag_expr = (
            mag_func_MG.merge(sample_pairs, left_on="Sample", right_on="MetaG_sample", how="inner")
                       .merge(
                            mag_func_MT.rename(columns={"Sample": "MetaT_sample"}),
                            on=["MAG_name", "Annotation", "MetaT_sample"],
                            how="inner"
                        )
        )

        mag_expr["log2_expr"] = np.log2(
            (mag_expr["Proportion_of_sample_coverage_MT"] + args.epsilon) /
            (mag_expr["Proportion_of_sample_coverage_MG"] + args.epsilon)
        )

        mag_expr = mag_expr[[
            "MAG_name", "Annotation", "MetaG_sample", "MetaT_sample",
            "MG_coverage_per_cell", "Proportion_of_sample_coverage_MG",
            "MT_coverage_per_cell", "Proportion_of_sample_coverage_MT",
            "log2_expr"
        ]]

        mag_expr.to_csv(
            os.path.join(args.output_dir, f"{args.output_prefix}.mags.functions.expression_profile.tsv"),
            sep="\t", index=False
        )
    else:
        if do_metat and (sample_pairs is None):
            print("[INFO] Skipping MAG-level expression (missing --sample_pairs).")
        elif not do_metat:
            print("[INFO] Skipping MAG-level MT + expression (no MetaT coverage provided).")
else:
    print("[INFO] No --mapping_file provided: skipping MAG-level functional profiling.")

print("[INFO] Functional profiling complete.")
