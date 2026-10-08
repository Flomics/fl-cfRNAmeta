"""PCA and tSNE of one dataset, coloured by every metadata variable in turn.

This is the per-dataset counterpart of the all-batches dimred: same colour maps,
same orderings, same transformations, one dataset at a time. Colours, marker
shapes and the dataset ordering come from src/dataset_mappings.json; the palettes
for every other variable come from src/dimred_plots.py.

    python src/dimred_perdataset.py \
        --matrix   results/cfRNA_meta_zhu/expression_matrix/gene_tpm_norm.tsv \
        --sampleinfo results/cfRNA_meta_zhu/sampleinfo.tsv \
        --label    zhu \
        --out-dir  dimred_plots_perdataset

Writes <out-dir>/<label>/{pca,tsne}/2d-{pca,tsne}_<variable>.{png,pdf}.

The matrix is a gene x sample table with gene_id and gene_name columns followed
by one numeric column per sample, as produced by snakeDA. gene_tpm_norm is the
matrix the figures use: it is normalised per sample, so a column subset of the
all-batches matrix is exactly the TPM of that dataset on its own. That does not
hold for TMM, whose factors are fitted across the whole sample set.
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dimred_plots import (  # noqa: E402
    discrete_color_map_dict,
    hue_order_map_dict,
    get_discrete_color_map,
    get_continous_color_map,
    new_scatter_dimred_plot,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MAPPINGS = os.path.join(REPO_ROOT, "src", "dataset_mappings.json")

RANDOM_STATE = 31415
ALPHA = 0.7

# tSNE is meaningless, and sklearn errors, below this many samples
MIN_SAMPLES_TSNE = 5
MAX_PERPLEXITY = 20

# Which components to plot against which, as 0-based pairs
TSNE_COMPONENT_PAIRS = [(0, 1)]
PCA_COMPONENT_PAIRS = [(0, 1), (1, 2), (3, 4)]

DISCRETE_PLOT_VARS = [
    # Identity
    "dataset_batch_label",
    "collection_center",
    # SRA
    "assay_type",
    "instrument",
    "organism",
    # Pre-analytical
    "plasma_tubes_short_name",
    "centrifugation_step_1",
    "centrifugation_step_2",
    "biomaterial",
    "nucleic_acid_type",
    "rna_extraction_kit_short_name",
    "dnase",
    "library_prep_kit_short_name",
    "library_selection",
    "cdna_library_type",
    "read_length",
    # Clinical
    "sex",
    "phenotype",
]

CONTINUOUS_PLOT_VARS = [
    # QC
    "exonic_percentage",
    "genes_contributing_to_80%_of_reads",
    "percentage_of_spliced_reads",
    "percentage_of_uniquely_mapped_reads",
    # Biotype
    "misc_rna_pct",
    "protein_coding_pct",
    "spike_in_pct",
    "mt_rna_pct",
    # Donor
    "age",
    "avgspotlen",
    # Cell-type deconvolution
    "platelet",
    "erythrocyte-erythroid_progenitor",
    "erythrocyte",
    "neutrophil",
]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--matrix", required=True,
                   help="gene x sample expression matrix (gene_id, gene_name, samples...)")
    p.add_argument("--sampleinfo", required=True,
                   help="sample metadata, one row per matrix sample column")
    p.add_argument("--label", default=None,
                   help="name of the subdirectory to write into (default: the matrix's dataset dir)")
    p.add_argument("--out-dir", default="dimred_plots_perdataset",
                   help="parent plots directory (default: %(default)s)")
    p.add_argument("--mappings", default=DEFAULT_MAPPINGS,
                   help="dataset labels, palette, markers and order (default: src/dataset_mappings.json)")
    p.add_argument("--marker-size", type=int, default=8,
                   help="scatter point area in points squared (default: %(default)s)")
    p.add_argument("--ext", nargs="+", default=[".png", ".pdf"],
                   help="output formats (default: .png .pdf)")
    return p.parse_args(argv)


def load_inputs(matrix_file, sampleinfo_file):
    """Read the matrix and the rows of the sampleinfo that describe its samples.

    The sampleinfo covers every sample of the study while the matrix is filtered, so
    it is subset and reordered to the matrix columns rather than required to match.
    """
    sample_df = pd.read_csv(sampleinfo_file, sep="\t", low_memory=False)
    counts_df = pd.read_csv(matrix_file, sep="\t", index_col="gene_id")

    if "sample_name" not in sample_df.columns:
        raise SystemExit("sampleinfo has no sample_name column to match the matrix on")

    sample_names = counts_df.select_dtypes(include=[np.number]).columns.tolist()
    n_rows = sample_df.shape[0]

    sample_df = sample_df.drop_duplicates("sample_name").set_index("sample_name")
    unknown = [ss for ss in sample_names if ss not in sample_df.index]
    if unknown:
        raise SystemExit(
            f"{len(unknown)} matrix columns are not in the sampleinfo, "
            f"e.g. {', '.join(unknown[:5])}")
    sample_df = sample_df.loc[sample_names].reset_index()

    print(f"  matrix     : {counts_df.shape[0]} genes x {len(sample_names)} samples")
    print(f"  sampleinfo : {sample_df.shape[1]} columns"
          + (f", {n_rows - len(sample_names)} rows not in the matrix dropped"
             if n_rows > len(sample_names) else ""))
    return counts_df, sample_df, sample_names


def add_dataset_mappings(sample_df, mappings_file):
    """Register the dataset colours, markers and ordering from dataset_mappings.json.

    The JSON is keyed by dataset_batch; the plots use the human-readable
    dataset_batch_label, so each mapping is re-keyed by label as well.
    """
    with open(mappings_file) as fh:
        mappings = json.load(fh)

    if "dataset_batch" not in sample_df.columns:
        raise SystemExit("sampleinfo has no dataset_batch column, which the plots colour by")

    labels = mappings["datasetsLabels"]
    sample_df["dataset_batch_label"] = sample_df["dataset_batch"].map(labels)

    discrete_color_map_dict["dataset_batch"] = mappings["datasetsPalette"]
    discrete_color_map_dict["dataset_batch_label"] = {
        labels[kk]: vv for kk, vv in mappings["datasetsPalette"].items()}

    # datasetsMarkers covers fewer datasets than datasetsLabels: the lab-level merged
    # names have no marker. An unmapped marker becomes NaN, and the scatter groups by
    # marker, which would drop those samples from the plot without saying so.
    markers = mappings["datasetsMarkers"]
    unmarked = [kk for kk in labels if kk not in markers]
    if unmarked:
        print(f"  no marker defined, using 'o': {', '.join(sorted(unmarked))}")
    marker_map_dict = {
        "dataset_batch": {kk: markers.get(kk, "o") for kk in labels},
        "dataset_batch_label": {labels[kk]: markers.get(kk, "o") for kk in labels},
    }

    hue_order_map_dict["dataset_batch"] = mappings["datasetVisualOrder"]
    hue_order_map_dict["dataset_batch_label"] = [
        labels[kk] for kk in mappings["datasetVisualOrder"]]

    return marker_map_dict


def build_color_maps(sample_df):
    """One colour map per metadata variable, skipping those the dataset lacks."""
    plot_vars = ({kk: "discrete" for kk in DISCRETE_PLOT_VARS}
                 | {kk: "continuous" for kk in CONTINUOUS_PLOT_VARS})

    missing = [kk for kk in plot_vars if kk not in sample_df.columns]
    if missing:
        plot_vars = {kk: vv for kk, vv in plot_vars.items() if kk not in missing}
        print(f"  not in sampleinfo, skipped: {', '.join(missing)}")

    color_map_dict = {}
    for var, kind in plot_vars.items():
        if kind == "discrete":
            color_map_dict[var] = get_discrete_color_map(
                sample_df, var, color_map_dict=discrete_color_map_dict)
        else:
            color_map_dict[var] = get_continous_color_map(sample_df, var)
    return color_map_dict


def log_transformed(counts_df, sample_names):
    """log(1 + x) on the sample columns, leaving gene_name alone."""
    out = counts_df.copy()
    out.loc[:, sample_names] = np.log(1 + out.loc[:, sample_names])
    return out


def usable_pairs(pairs, n_samples, n_features):
    """Drop component pairs the dataset is too small to produce.

    PCA yields at most min(n_samples, n_features) components, so the small batches
    cannot reach PC5 and sklearn raises rather than returning what it can.
    """
    n_max = min(n_samples, n_features)
    keep = [(a, b) for a, b in pairs if max(a, b) < n_max]
    dropped = [(a, b) for a, b in pairs if (a, b) not in keep]
    if dropped:
        print(f"  only {n_max} components available, dropping pairs: "
              + ", ".join(f"{a + 1}/{b + 1}" for a, b in dropped))
    return keep


def run_dimred(method, plot_df, sample_names, color_map_dict, marker_map_dict,
               out_dir, pairs, params_config, args):
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n{method}: {len(sample_names)} samples, {plot_df.shape[0]} features -> {out_dir}")
    pc1 = [a for a, _ in pairs]
    pc2 = [b for _, b in pairs]
    new_scatter_dimred_plot(
        df=plot_df,
        dimred_method=method,
        params_config=params_config,
        pc1=pc1,
        pc2=pc2,
        sample_ids=sample_names,
        annotate_plot=False,
        full_matrix=False,
        color_by=color_map_dict.keys(),
        color_map_dict=color_map_dict,
        marker_map_dict=marker_map_dict,
        hue_order_dict=hue_order_map_dict,
        plot_ext=args.ext,
        marker_size=args.marker_size,
        out_dir=out_dir,
        alpha=ALPHA,
    )


def main(argv=None):
    args = parse_args(argv)
    np.random.seed(RANDOM_STATE)

    label = args.label or os.path.basename(
        os.path.dirname(os.path.dirname(os.path.abspath(args.matrix))))
    print(f"=== {label} ===")

    counts_df, sample_df, sample_names = load_inputs(args.matrix, args.sampleinfo)
    marker_map_dict = add_dataset_mappings(sample_df, args.mappings)
    color_map_dict = build_color_maps(sample_df)

    plot_df = log_transformed(counts_df, sample_names)
    dataset_dir = os.path.join(args.out_dir, label)

    n_samples, n_features = len(sample_names), plot_df.shape[0]

    # tSNE: perplexity must stay below the sample count
    perplexity = min(n_samples - 1, MAX_PERPLEXITY)
    if perplexity < MIN_SAMPLES_TSNE:
        print(f"\ntsne: only {n_samples} samples, skipped (PCA only)")
    else:
        run_dimred("tsne", plot_df, sample_names, color_map_dict, marker_map_dict,
                   os.path.join(dataset_dir, "tsne"),
                   usable_pairs(TSNE_COMPONENT_PAIRS, n_samples, n_features),
                   {"perplexity": perplexity}, args)

    pca_pairs = usable_pairs(PCA_COMPONENT_PAIRS, n_samples, n_features)
    if pca_pairs:
        run_dimred("pca", plot_df, sample_names, color_map_dict, marker_map_dict,
                   os.path.join(dataset_dir, "pca"), pca_pairs, None, args)
    else:
        print(f"\npca: only {n_samples} samples, skipped")

    print(f"\ndone: {dataset_dir}/{{pca,tsne}}/")


if __name__ == "__main__":
    main()
