#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cramer's V between phenotype and each pre-analytical variable, per dataset (Figure 1C).

Usage:
    python src/metadata_correlations.py <per_sample_metadata.tsv> <out_dir>

    per_sample_metadata.tsv   tables/cfRNA-meta_per_sample_metadata.tsv
    out_dir                   writes metadata_corrs.{png,pdf,svg} and metadata_corrs.tsv

Created on Mon Jun 30 12:01:55 2025
@author: erusu
"""
import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
import seaborn as sns

if len(sys.argv) != 3:
    sys.exit(__doc__)
FILE, OUTPATH = sys.argv[1:]
os.makedirs(OUTPATH, exist_ok=True)

# The panel must be in Arial, which cannot be redistributed here. Fail rather than let
# matplotlib fall back to another font with only a log warning.
try:
    font_manager.findfont("Arial", fallback_to_default=False)
except ValueError:
    sys.exit("Arial is not installed: the figure must be drawn in Arial. Install it and "
             "clear the matplotlib font cache (matplotlib.get_cachedir()).")

df = pd.read_csv(FILE, sep ="\t")

# modify phenotype column
df['simple_phenotype'] = pd.Series(np.nan, index=df.index, dtype=object)

df.loc[df['phenotype'] == "healthy", 'simple_phenotype'] = "healthy"
df.loc[df['phenotype'].isna(), 'simple_phenotype'] = "missing"
df.loc[df['phenotype'] == "", 'simple_phenotype'] = "missing"
df.loc[df['phenotype'].str.contains("Acute Myeloid Leukemia", na=False), 'simple_phenotype'] = "cancer"
df.loc[df['phenotype'].str.contains("Alzheimers disease", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Chronic hepatitis B", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Chronic kidney failure EPO-treated", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Cirrhosis", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Colorectal cancer", na=False), 'simple_phenotype'] = "cancer"
df.loc[df['phenotype'].str.contains("Diffuse large B-cell lymphoma", na=False), 'simple_phenotype'] = "cancer"
df.loc[df['phenotype'].str.contains("Diverticulitis", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Esophagus cancer", na=False), 'simple_phenotype'] = "cancer"
df.loc[df['phenotype'].str.contains("G-CSF-treated healthy donors", na=False), 'simple_phenotype'] = "healthy"
df.loc[df['phenotype'].str.contains("Healthy", na=False), 'simple_phenotype'] = "healthy"
df.loc[df['phenotype'].str.contains("Healthy pregnant woman", na=False), 'simple_phenotype'] = "healthy"
df.loc[df['phenotype'].str.contains("Healthy pregnant woman who delivered preterm", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Liver cancer", na=False), 'simple_phenotype'] = "cancer"
df.loc[df['phenotype'].str.contains("Lung cancer", na=False), 'simple_phenotype'] = "cancer"
df.loc[df['phenotype'].str.contains("Multiple myeloma", na=False), 'simple_phenotype'] = "cancer"
df.loc[df['phenotype'].str.contains("Nonalcoholic fatty liver disease", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Nonalcoholic steatohepatitis", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Pancreatic cancer", na=False), 'simple_phenotype'] = "cancer"
df.loc[df['phenotype'].str.contains("Pre-cancerous condition: cirrhosis", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Pre-cancerous condition: MGUS", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Pre-eclampsia", na=False), 'simple_phenotype'] = "non-cancer disease"
df.loc[df['phenotype'].str.contains("Primary mediastinal B-cell lymphoma", na=False), 'simple_phenotype'] = "cancer"
df.loc[df['phenotype'].str.contains("Stomach cancer", na=False), 'simple_phenotype'] = "cancer"



MAX_MISSING = 0.5

def cramers_v_corrected(x, y):
    # a few samples without the variable should not hide it for the whole dataset:
    # they are left out, and only a variable missing in most samples is Unspecified
    if y.isnull().mean() > MAX_MISSING:
        return "Unspecified"
    known = y.notna() & (x != "missing")
    confusion_matrix = pd.crosstab(x[known], y[known])
    # V is undefined if the samples left have a single value or a single phenotype
    if min(confusion_matrix.shape) < 2:
        return "Uniform"

    chi2, _, _, _ = chi2_contingency(confusion_matrix, correction=False)
    n = confusion_matrix.to_numpy().sum()
    r, k = confusion_matrix.shape
    return np.sqrt((chi2 / n) / min(k - 1, r - 1))



df['new_dataset_short_name'] = df['dataset_short_name']
df.loc[( df['dataset_short_name'] == "ibarra") & (df['biomaterial'] == "buffy coat"), "new_dataset_short_name"] = "Ibarra (buffy coat)"
df.loc[( df['dataset_short_name'] == "ibarra") & (df['biomaterial'] == "serum"), "new_dataset_short_name"] = "Ibarra (serum)"
df.loc[( df['dataset_short_name'] == "ibarra") & (df['biomaterial'] == "plasma"), "new_dataset_short_name"] = "Ibarra (plasma)"

# datasets without second centrifugation, should not be marked as missing info in this plot
df.loc[df['dataset_short_name'] == "block", "centrifugation_step_2"] = "None"
df.loc[df['dataset_short_name'] == "chalasani", "centrifugation_step_2"] = "None"
df.loc[df['dataset_short_name'] == "decruyenaere", "centrifugation_step_2"] = "None"
df.loc[df['new_dataset_short_name'] == "Ibarra (buffy coat)", "centrifugation_step_2"] = "None"
df.loc[(df['dataset_short_name'] == "moufarrej") & (df['centrifugation_step_2'].isna()), "centrifugation_step_2"] = "None"
df.loc[df['dataset_short_name'] == "toden", "centrifugation_step_2"] = "None"


# datasets without DNA treatments should not be marked as missing info in this plot
df.loc[df['dataset_short_name'] == "block", "dnase"] = "None"
df.loc[df['dataset_short_name'] == "reggiardo", "dnase"] = "None"
df.loc[df['dataset_short_name'] == "toden", "dnase"] = "None"


# all sun samples come from the same hospital
df.loc[df['dataset_short_name'] == "sun", "collection_center"] = "Zhongnan Hospital"
df.loc[df['dataset_short_name'] == "tao", "collection_center"] = "Peking University First Hospital"

# Toden has one sample with collection center "Unknown"
df.loc[(df['dataset_short_name'] == "toden") & (df["collection_center"] == "Unknown"), "collection_center"] = np.nan


# block & moufarrej blood collection tube not known
df.loc[(df['dataset_short_name'] == "block"), "plasma_tubes_short_name"] = np.nan
df.loc[(df['dataset_short_name'] == "moufarrej") & (df["plasma_tubes_short_name"] =="Unspecified"), "plasma_tubes_short_name"] = np.nan
df.loc[(df['dataset_short_name'] == "toden"), "plasma_tubes_short_name"] = np.nan
df.loc[(df['dataset_short_name'] == "zhu"), "plasma_tubes_short_name"] = np.nan

# set centrifugation NAs 
df.loc[(df['dataset_short_name'] == "chen"), "centrifugation_step_1"] = np.nan
df.loc[(df['dataset_short_name'] == "chen"), "centrifugation_step_2"] = np.nan
df.loc[(df['dataset_short_name'] == "reggiardo"), "centrifugation_step_1"] = np.nan
df.loc[(df['dataset_short_name'] == "reggiardo"), "centrifugation_step_2"] = np.nan
df.loc[(df['dataset_short_name'] == "tao"), "centrifugation_step_1"] = np.nan
df.loc[(df['dataset_short_name'] == "tao"), "centrifugation_step_2"] = np.nan
df.loc[(df['dataset_short_name'] == "zhu"), "centrifugation_step_1"] = np.nan
df.loc[(df['dataset_short_name'] == "zhu"), "centrifugation_step_2"] = np.nan

# set library prep kit NAs
df.loc[(df['dataset_short_name'] == "chalasani"), "library_prep_kit_short_name"] = np.nan
df.loc[(df['dataset_short_name'] == "ibarra"), "library_prep_kit_short_name"] = np.nan

COL = "new_dataset_short_name"

datasets = df[COL].unique()

tech_vars =["collection_center", "plasma_tubes_short_name", "centrifugation_step_1", 
                 "centrifugation_step_2", "biomaterial", "nucleic_acid_type",
                 "rna_extraction_kit_short_name", "dnase", "library_prep_kit_short_name",
                 "library_selection", "cdna_library_type", "read_length"]

# datasets excluded from the plot because all phenotypes are the same
excluded_datasets = ["flomics_1", "flomics_2", "giraldez", "ngo", "wang", 
                     "rozowsky", "wei"]

#excluded_datasets = ['rozowsky','giraldez_standard', 'wei',"ngo",
#                     'giraldez_phospho-rna-seq', 'flomics_2','wang','flomics_1']

results = pd.DataFrame(index=tech_vars, columns=datasets)
# cells computed after leaving out samples without the variable, marked * in the plot
partial = pd.DataFrame(False, index=tech_vars, columns=datasets)

for dset in datasets:
    if dset not in excluded_datasets: 
        sub_df = df[df[COL] == dset]
        for var in tech_vars:
            try:
                score = cramers_v_corrected(sub_df['simple_phenotype'], sub_df[var])
            except:
                score = np.nan
            results.loc[var, dset] = score
            partial.loc[var, dset] = (isinstance(score, float)
                                      and sub_df[var].isnull().any())
        
# reorder columns accordingly
dataset_order = ["sun", "chalasani", "Ibarra (buffy coat)", 
                 "Ibarra (plasma)", "Ibarra (serum)", "toden", "reggiardo", "block",
                 "chen", "decruyenaere",
                 "moufarrej", "roskams", "tao",  "zhu"]

#dataset_order = ['block_150bp', 'block_300bp','chalasani','chen','decruyenaere', 
#                 'ibarra_buffy_coat', 'ibarra_plasma_cancer', 'ibarra_plasma_non_cancer',
#                 'ibarra_serum', 'moufarrej_site_1', 'moufarrej_site_2', 
#                 'reggiardo_bioivt', 'reggiardo_dls','roskams_pilot', 'roskams_validation',
#                  'sun_2', 'tao', 'toden', 'zhu']


results = results[dataset_order]
partial = partial[dataset_order]
results.to_csv(os.path.join(OUTPATH, "metadata_corrs.tsv"), sep="\t")
print(results)
# drop rows that are all the same: 
sel_results = results.drop(['biomaterial', "nucleic_acid_type", 
                            "rna_extraction_kit_short_name", "dnase", 
                            "library_prep_kit_short_name",
                            "plasma_tubes_short_name",
                            "library_selection", "cdna_library_type"], axis=0)
partial = partial.loc[sel_results.index]
partial.index, partial.columns = range(len(partial.index)), range(len(partial.columns))


# set errors to ignore in case we remove columns 
sel_results = sel_results.rename(index={'collection_center': "Collection center", 
                                        'plasma_tubes_short_name': "Blood collection tube",
                                        'centrifugation_step_1': "Centrifugation, step 1",
                                        'centrifugation_step_2': "Centrifugation, step 2", 
                                        'biomaterial': "Biomaterial", 
                                        'nucleic_acid_type': "Nucleic acid type",
                                        'rna_extraction_kit_short_name': "RNA extraction kit", 
                                        'dnase': "DNAse treatment", 
                                        'library_prep_kit_short_name': "Library prep kit",
                                        'library_selection': "Library selection", 
                                        'cdna_library_type': "cDNA library type", 
                                        'read_length': "Read length"}, errors = "ignore")


sel_results.columns = sel_results.columns.map(str.title)
sel_results = sel_results.rename(columns={'Ibarra (Buffy Coat)': "Ibarra (buffy coat)", 
                                        'Ibarra (Plasma)': "Ibarra (plasma)",
                                        'Ibarra (Serum)': "Ibarra (serum)", 
                                         'Roskams': 'Roskams-Hieter'})
cmap = sns.cubehelix_palette(start=0.4, rot=0, dark=0.4, light=0.91, as_cmap=True)

# Build color matrix
def get_color(val):
    if val == "Uniform":
        return "#B2DF8A"
    elif val == "Unspecified":
        return "#333333"
    else:
        try:
            # Normalize number to colormap
            return cmap(float(val))
        except:
            return "white"  # fallback for unexpected content

color_matrix = sel_results.map(get_color)


fonts = 5
fonts_cells = 4

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["svg.fonttype"] = "none"

# Layout in inches, so that cell size and margins do not depend on the number of
# datasets and the legend and colour bar never overlap the heatmap or each other
n_rows, n_cols = sel_results.shape
cell = 0.22
left, bottom, top, right = 1.0, 0.85, 0.05, 0.9
hm_w, hm_h = n_cols * cell, n_rows * cell
fig_w, fig_h = left + hm_w + right, bottom + hm_h + top
fig = plt.figure(figsize=(fig_w, fig_h))
ax = fig.add_axes([left / fig_w, bottom / fig_h, hm_w / fig_w, hm_h / fig_h])

for i, row in enumerate(sel_results.index):
    for j, col in enumerate(sel_results.columns):
        val = sel_results.loc[row, col]
        ax.add_patch(plt.Rectangle((j, i), 1, 1, color=color_matrix.loc[row, col]))
        if not isinstance(val, float):
            continue
        if partial.loc[i, j]:
            label = "{}*".format(round(val, 1))
        elif row == "Centrifugation, step 2" and col == "Ibarra (plasma)":
            label = "{}**".format(round(val, 1))
        else:
            label = "{}".format(round(val, 1))
        ax.text(j + 0.5, i + 0.5, label, ha='center', va='center', fontsize=fonts_cells,
                color="white" if val > 0.5 else "black")

ax.set_xlim(0, n_cols)
ax.set_ylim(0, n_rows)
ax.set_xticks(np.arange(n_cols) + 0.5)
ax.set_xticklabels(sel_results.columns, fontsize=fonts, rotation=45, ha='right',
                   rotation_mode='anchor')
ax.set_yticks(np.arange(n_rows) + 0.5)
ax.set_yticklabels(sel_results.index, fontsize=fonts)
ax.invert_yaxis()
ax.tick_params(axis='both', which='both', length=0, pad=2)

for spine in ax.spines.values():
    spine.set_visible(False)
for i in range(n_rows + 1):
    ax.hlines(i, 0, n_cols, color='white', linewidth=1)
for j in range(n_cols + 1):
    ax.vlines(j, 0, n_rows, color='white', linewidth=1)

# legend at the top right of the heatmap, colour bar below it
legend_elements = [
    Patch(facecolor='#B2DF8A', edgecolor='black', linewidth=0.3, label='Uniform'),
    Patch(facecolor='#333333', edgecolor='black', linewidth=0.3, label='Unspecified'),
]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1),
          fontsize=fonts, frameon=False, borderaxespad=0., handlelength=1, handleheight=1)

cbar_h = 0.45
cbar_ax = fig.add_axes([(left + hm_w + 0.1) / fig_w, bottom / fig_h, 0.06 / fig_w, cbar_h / fig_h])
sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=0, vmax=1))
cbar = fig.colorbar(sm, cax=cbar_ax, ticks=[0, 0.5, 1])
cbar.outline.set_linewidth(0.3)
cbar.ax.tick_params(labelsize=fonts, width=0.3, length=1.5, pad=1)
cbar.set_label("Cramér's V", fontsize=fonts, labelpad=2)

plt.savefig(os.path.join(OUTPATH, "metadata_corrs.png"), dpi=600, bbox_inches='tight', pad_inches=0.02)
plt.savefig(os.path.join(OUTPATH, "metadata_corrs.pdf"), dpi=600, bbox_inches='tight', pad_inches=0.02)
plt.savefig(os.path.join(OUTPATH, "metadata_corrs.svg"), dpi=600, bbox_inches='tight', pad_inches=0.02)

