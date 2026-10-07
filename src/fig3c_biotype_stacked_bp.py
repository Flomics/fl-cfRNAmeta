"""Figure 3C: mean RNA biotype representation per dataset batch.

Usage:
    python src/fig3c_biotype_stacked_bp.py <gene_tpm_norm.tsv> <sampleinfo.tsv> <gene_biotypes.tsv> <out_dir>

    python src/fig3c_biotype_stacked_bp.py tables/gene_tpm_norm.tsv tables/sampleinfo_all-batches.tsv \
        tables/gencode_v39_gene_biotypes.tsv figures/fig3c

Each bar is the mean over the batch's samples of the per-sample biotype proportions of TPM,
spike-ins excluded. Ported from the gene_biotype_analysis notebook in fl-snakeDA.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SRC_DIR)
from barplot_stacked_plot import create_stacked_barplot

if len(sys.argv) != 5:
    sys.exit(__doc__)
count_data_file, sample_info_file, gene_biotypes_file, plots_dir = sys.argv[1:5]

# The panel must be in Arial, which cannot be redistributed here. Fail rather than let
# matplotlib fall back to another font with only a log warning.
try:
    font_manager.findfont("Arial", fallback_to_default=False)
except ValueError:
    sys.exit("Arial is not installed: the figure must be drawn in Arial. Install it and "
             "clear the matplotlib font cache (matplotlib.get_cachedir()).")

os.makedirs(plots_dir, exist_ok=True)

biotype_col = "biotype_cfRNAmeta"
group_var = "dataset_batch_label"
boxes_var = "broad_protocol_category_short"
bracket_var = "dataset_short_name"

with open(os.path.join(SRC_DIR, "dataset_mappings.json")) as fh:
    dataset_mappings = json.load(fh)
with open(os.path.join(SRC_DIR, "gene_biotype_mappings.json")) as fh:
    biotype_mappings = json.load(fh)

# The published panel predates the split of Mt_tRNA into its own "MT tRNA" category (5b7fa12).
biotype_map = dict(biotype_mappings["cfrnameta_biotype_map"], Mt_tRNA="Other small RNAs")

broad_protocol_category = {
    batch: group for group, batches in dataset_mappings["datasetGroups"].items() for batch in batches
}
broad_protocol_category_short = {"custom": "Custom", "cfDNA": "cfDNA", "exome_based": "EB", "wro": "WRO", "wrr": "WRR"}
broad_protocol_category_color = {
    "Custom": "#D9BBAE",
    "cfDNA": "#AECAD9",
    "EB": "#BBE0BB",
    "WRO": "#D0AED9",
    "WRR": "#D9D6AE",
}

counts_df = pd.read_csv(count_data_file, sep="\t")
sample_names = counts_df.select_dtypes(include=[np.number]).columns.tolist()
counts_df = counts_df[~counts_df["gene_id"].str.startswith(("ERCC", "SIRV"))]

# The sampleinfo covers more samples than the filtered matrix.
sample_df = pd.read_csv(sample_info_file, sep="\t", low_memory=False).set_index("sample_name")
missing = sorted(set(sample_names) - set(sample_df.index))
if missing:
    sys.exit(f"{len(missing)} matrix samples are missing from the sampleinfo, e.g. {missing[:3]}")
sample_df = sample_df.loc[sample_names].copy()
sample_df[group_var] = sample_df["dataset_batch"].map(dataset_mappings["datasetsLabels"])
sample_df[boxes_var] = sample_df["dataset_batch"].map(broad_protocol_category).map(broad_protocol_category_short)
unmapped = sample_df.loc[sample_df[[group_var, boxes_var]].isna().any(axis=1), "dataset_batch"].unique()
if len(unmapped):
    sys.exit(f"dataset_batch values without a label or protocol category: {list(unmapped)}")

gene_df = pd.read_csv(gene_biotypes_file, sep="\t")
gene_df[biotype_col] = gene_df["custom_biotype"].map(biotype_map)
if gene_df[biotype_col].isna().any():
    sys.exit(f"custom_biotype values without a mapping: {sorted(gene_df.loc[gene_df[biotype_col].isna(), 'custom_biotype'].unique())}")

# Colours follow the published panel: tab20 spread over the categories, most genes first.
biotypes = gene_df[biotype_col].value_counts().index.tolist()
colors_dict = {bb: plt.cm.tab20(vv) for bb, vv in zip(biotypes, np.linspace(0, 1, len(biotypes)))}

biotype_tpm = (
    pd.concat(
        [gene_df.set_index("gene_id")[[biotype_col]], counts_df.set_index("gene_id")[sample_names]],
        join="inner",
        axis=1,
    )
    .groupby(biotype_col)
    .sum()
    .T
)
print(f"{len(counts_df)} genes x {len(sample_names)} samples")

df = pd.concat([biotype_tpm, sample_df[[group_var, boxes_var, bracket_var]]], axis=1)
ordered_cols = [cc for cc in biotype_mappings[biotype_col] if cc in df.columns]
present = set(df[group_var])
ordered_rows = [
    dataset_mappings["datasetsLabels"][bb]
    for bb in dataset_mappings["datasetVisualOrder"]
    if dataset_mappings["datasetsLabels"][bb] in present
]

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Arial"
plt.rcParams["mathtext.it"] = "Arial:italic"
plt.rcParams["mathtext.bf"] = "Arial:bold"
font_size = 5

out_prefix = os.path.join(plots_dir, "fig3c_biotype_stacked_bp")
_, df_mean, _, _ = create_stacked_barplot(
    dataset=df,
    meta_column=group_var,
    value_columns=ordered_cols,
    meta_order=ordered_rows,
    color_map=colors_dict,
    bar_width=0.8,
    fig_width=5,
    fig_height=2.5,
    ax_width=1.01 * 2.99,
    aspect=2,
    dpi=600,
    output=out_prefix,
    scaling="mean",
    normalize_data=True,
    add_category_border=True,
    category_border_width=0.1,
    boxes_column=boxes_var,
    boxes_color_map=broad_protocol_category_color,
    boxes_y_position=-0.04,
    boxes_height=0.03,
    boxes_width=0.95,
    boxes_borderwidth=0.0,
    boxes_legend_pos="bottom",
    boxes_legend_title="Broad Protocol Category (BPC)",
    boxes_legend_fontsize=font_size,
    boxes_legend_y_pos=-0.65,
    group_by_column=bracket_var,
    group_spacing=0.0,
    group_label_y_offset=-0.065,
    group_bracket_linewidth=0.5,
    group_bracket_vertical_line_length=0.01,
    group_bracket_horizontal_line_length=0.4,
    group_position="middle",
    show_group_label=False,
    show_xlabel=False,
    xlabel_fontsize=font_size,
    ylabel_fontsize=font_size,
    hide_bottom_tick=True,
    hide_left_tick=True,
    x_tick_label_fontsize=font_size,
    y_tick_label_fontsize=font_size,
    x_tick_label_rotation=45,
    x_ticks_label_pad=-0.035,
    show_title=False,
    legend_title="Biotype",
    legend_x_pos=1.0,
    legend_y_pos=0.5,
    legend_title_fontsize=font_size,
    legend_fontsize=font_size,
    hide_top_spine=True,
    hide_right_spine=True,
    hide_bottom_spine=True,
    hide_left_spine=True,
    y_upper_pad=0.05,
)

# Proportions per bar, as plotted
df_mean.div(df_mean.sum(axis=1), axis=0).loc[ordered_rows].to_csv(f"{out_prefix}_proportions.tsv", sep="\t")
