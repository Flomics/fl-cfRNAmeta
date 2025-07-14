#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 12:01:55 2025

@author: erusu
"""
import os
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
import seaborn as sns

FILE = "../tables/cfRNA-meta_per_sample_metadata.tsv"

df = pd.read_csv(FILE, sep ="\t")

OUTPATH = "/mnt/efs/home/erusu/workplace/2025-05-cfRNAmeta/5.phenotype_corr"

# modify phenotype column
df['simple_phenotype'] = np.nan

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



def cramers_v_corrected(x, y):
    confusion_matrix = pd.crosstab(x, y, dropna = False)
    
    # if there is NaN in the columns, return Unspecified
    if y.isnull().all(): 
        return "Unspecified"
    elif np.nan in list(confusion_matrix.columns): 
        return "Unspecified"
    else: 
        if confusion_matrix.shape[1] < 2:
            return "Uniform"  
        
        chi2, _, _, _ = chi2_contingency(confusion_matrix, correction=False)
        n = confusion_matrix.to_numpy().sum()
        r, k = confusion_matrix.shape
        result = np.sqrt((chi2 / n) / min(k - 1, r - 1))

        return result



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

# Toden --> remove only 1 not annotated sample to be able to calculate it
df = df[~((df['dataset_short_name'] == "toden") & (df["collection_center"] == "Unknown"))]
#df.loc[(df['dataset_short_name'] == "toden") & (df["collection_center"] == "Unknown"), "collection_center"] = np.nan


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

for dset in datasets:
    if dset not in excluded_datasets: 
        sub_df = df[df[COL] == dset]
        for var in tech_vars:
            try:
                score = cramers_v_corrected(sub_df['simple_phenotype'], sub_df[var])
            except:
                score = np.nan
            results.loc[var, dset] = score
        
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
print(results)
# drop rows that are all the same: 
sel_results = results.drop(['biomaterial', "nucleic_acid_type", 
                            "rna_extraction_kit_short_name", "dnase", 
                            "library_prep_kit_short_name",
                            "plasma_tubes_short_name",
                            "library_selection", "cdna_library_type"], axis=0)


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
# Build color matrix
def get_color(val):
    if val == "Uniform":
        return "#B2DF8A"
    elif val == "Unspecified":
        return "#333333"
    else:
        try:
            # Normalize number to colormap
            norm = Normalize(vmin=0, vmax=1)
            cmap = sns.cubehelix_palette(start = 0.4, rot = 0, dark = 0.4, 
                                         light = 0.91, as_cmap = True)
            return cmap(norm(float(val)))
        except:
            return "white"  # fallback for unexpected content

color_matrix = sel_results.applymap(get_color)


fonts = 8 # original
fonts = 5 
fonts_cells = 4
# Plot using matplotlib
#fig, ax = plt.subplots(figsize=(6, 9))
fig, ax = plt.subplots(figsize=(2.5, 3.5))  # Add this cell_size line above or hardcode height

for i, row in enumerate(sel_results.index):
    for j, col in enumerate(sel_results.columns):
        val = sel_results.loc[row, col]
        color = color_matrix.loc[row, col]
        ax.add_patch(plt.Rectangle((j, i), 1, 1, color=color))
        tcol = "black"
        if isinstance(val, float): 
            if val > 0.5: 
                tcol = "white"
        # Display values
        if row == "Collection center" and col == "Toden": 
            ax.text(j + 0.5, i + 0.5, "{}*".format(round(val, 1)) if isinstance(val, float) else "",
                    ha='center', va='center', fontsize=fonts_cells, color = tcol) 
        elif row == "Centrifugation, step 2" and col == "Ibarra (plasma)": 
            ax.text(j + 0.5, i + 0.5, "{}**".format(round(val, 1)) if isinstance(val, float) else "",
                    ha='center', va='center', fontsize=fonts_cells, color = tcol)             
        else:
            ax.text(j + 0.5, i + 0.5, "{}".format(round(val, 1)) if isinstance(val, float) else "",
                    ha='center', va='center', fontsize=fonts_cells, color = tcol)

ax.set_xlim(0, len(sel_results.columns))
ax.set_ylim(0, len(sel_results.index))
ax.set_xticks(np.arange(len(sel_results.columns)) + 0.5)
ax.set_xticklabels(sel_results.columns, fontsize=fonts)

ax.set_yticks(np.arange(len(sel_results.index)) + 0.5)
ax.set_yticklabels(sel_results.index, fontsize=fonts)
ax.invert_yaxis()
ax.set_aspect("equal")

ax.tick_params(axis='both', which='both', length=0)

plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Add colorbar for numeric values

# Define manual colorbar position: [left, bottom, width, height]
cbar_ax = fig.add_axes([1.1, 0.38, 0.02, 0.1])  # adjust these values to shift/resize
sm = ScalarMappable(cmap=sns.cubehelix_palette(start=0.4, rot=0, dark=0.4, light=0.91, 
                                               as_cmap=True), 
                    norm=Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)

# Add manual label
cbar_ax.text(-1.7, 0.5, "Cramer's V score", rotation=90, va='center', fontsize=fonts,
             ha='center', transform=cbar_ax.transAxes)





# Define manual legend entries
legend_elements = [
    Patch(facecolor='#B2DF8A', edgecolor='black', label='Uniform'),
    Patch(facecolor='#333333', edgecolor='black', label='Unspecified'),
]

# Add to the plot
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), 
          fontsize=fonts, borderaxespad=0.)

# remove spines
for spine in ax.spines.values():
    spine.set_visible(False)

# Add grid lines
for i in range(len(sel_results.index) + 1):
    ax.hlines(i, 0, len(sel_results.columns), color='white', linewidth=1)

for j in range(len(sel_results.columns) + 1):
    ax.vlines(j, 0, len(sel_results.index), color='white', linewidth=1)


plt.savefig(os.path.join(OUTPATH, "metadata_corrs2.png"), dpi=600, bbox_inches='tight')
plt.savefig(os.path.join(OUTPATH, "metadata_corrs2.pdf"), dpi=600, bbox_inches='tight')
plt.savefig(os.path.join(OUTPATH, "metadata_corrs2.svg"), dpi=600, bbox_inches='tight')

