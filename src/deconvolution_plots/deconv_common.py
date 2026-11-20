"""Shared loading code for the deconvolution figure notebooks (Fig 4A/B, Suppl Figs 5/6).

Same steps as the "POST DECONVOLUTION ANALYSIS" cells of src/deconvolution_results.ipynb:
read the nu-SVR fractions of each dataset, join them with the per-sample metadata, add the
plotting batches. Paths are resolved from the repository root, so the notebooks run from
any working directory.

Inputs per dataset, as written by the nu-SVR script (62 cell types + r + rmse): a fractions
CSV and a sample list. Their paths are given as templates with a {dataset} placeholder.
"""
import io
import json
import os
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"

# gDNA-free datasets, keyed as in the deconvolution output folders
DATASETS = ("Flomics_2", "chalasani_merged", "chen", "ibarra", "moufarrej", "reggiardo",
            "roskams", "tao", "toden", "zhu", "decru")

# metadata dataset_short_name -> key above
DATASET_CONV = {"flomics_2": "Flomics_2", "chalasani": "chalasani_merged", "decruyenaere": "decru"}

PLASMA_TUBE_DICT = {"EDTA": "EDTA", "EDTA/Streck": "EDTA/Streck",
                    "BD Vacutainer clotting tubes": "Vacutainer", "K2EDTA": "EDTA",
                    "Unspecified": "Unspecified"}

CELL_TYPE_COLORS = {  # cell_type_colors_dark in the original notebook
    "erythrocyte": "#DC143C", "platelet": "#B8860B", "monocyte": "#FFA500", "b cell": "#4169E1",
    "duct epithelial cell": "#808000", "neutrophil": "#ADD8E6", "myeloid progenitor": "#DAA520",
    "goblet cell": "#6A5ACD", "nk cell": "#9370DB", "basophil": "#4B0082",
    "endothelial cell": "#008080", "macrophage": "#B22222", "innate lymphoid cell": "#8A2BE2",
    "kidney epithelial cell": "#32CD32", "fibroblast": "#8B4513", "mesothelial cell": "#2E8B57",
}

BPC_SHORT = {"Custom": "Custom", "cfDNA": "cfDNA", "Exome-based (EB)": "EB",
             "Whole RNA-Seq (oligo-dT pr.) (WRO)": "WRO", "Whole RNA-Seq (random pr.) (WRR)": "WRR"}
BPC_COLORS = {"Custom": "#D9BBAE", "cfDNA": "#AECAD9", "EB": "#BBE0BB", "WRO": "#D0AED9", "WRR": "#D9D6AE"}

# x-axis order of Fig 4A and Suppl Fig 5 (labels from src/dataset_mappings.json)
BATCH_ORDER = ["Chalasani", "Ibarra (buffy coat)", "Ibarra (plasma C)", "Ibarra (plasma NC)",
               "Ibarra (serum)", "Toden", "Reggiardo (BioIVT)", "Reggiardo (DLS)",
               "Block (2x75bp)", "Block (2x150bp)", "Chen (Peking University)", "Chen (Others)",
               "Decruyenaere", "Flomics (Cordoba)", "Flomics (Others)", "Moufarrej (Site 1)",
               "Moufarrej (Site 2)", "Roskams (pilot)", "Roskams (validation)", "Tao",
               "Zhu (Cancer)", "Zhu (healthy)"]


def set_fonts(family="Arial"):
    """Arial as in the paper; Liberation Sans (metric-compatible) where Arial is missing."""
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [family, "Liberation Sans", "DejaVu Sans"]


def git_show(rev, path):
    """Contents of <path> at commit <rev> of this repository."""
    return subprocess.check_output(["git", "-C", str(REPO), "show", f"{rev}:{path}"], text=True)


def load_metadata(git_commit=None):
    """Per-sample metadata, indexed by run. git_commit: read the table as it was at that commit (default: current file)."""
    path = "tables/cfRNA-meta_per_sample_metadata.tsv"
    src = io.StringIO(git_show(git_commit, path)) if git_commit else REPO / path
    meta = pd.read_csv(src, sep="\t", low_memory=False)
    meta["dataset_short_name"] = meta["dataset_short_name"].replace(DATASET_CONV)
    return meta.set_index("run")


def load_labels(git_commit=None):
    """datasetsLabels from src/dataset_mappings.json (as it was at git_commit, if given)."""
    path = "src/dataset_mappings.json"
    txt = git_show(git_commit, path) if git_commit else (REPO / path).read_text()
    return json.loads(txt)["datasetsLabels"]


def load_deconvolution(metadata, fractions_file, sample_list_file, datasets=DATASETS):
    """Fractions joined with metadata, one DataFrame per dataset key.

    fractions_file, sample_list_file: paths with a {dataset} placeholder, filled with each key
    of `datasets`. Returns (df_dict, cell_types). Datasets without a fractions file are skipped.
    """
    df_dict, cell_types = {}, None
    for key in datasets:
        frac_path = Path(fractions_file.format(dataset=key))
        if not frac_path.exists():
            print(f"{key}: NOT found ({frac_path})")
            continue
        frac = pd.read_csv(frac_path, index_col=0)
        samples = Path(sample_list_file.format(dataset=key)).read_text().splitlines()
        frac.columns = frac.columns.str.split("-").str[0]  # "<run>-NUSVR--nu=..-C=.." -> run
        assert frac.columns.is_unique, f"{key}: duplicated sample columns"
        assert set(frac.columns) == set(samples), f"{key}: samples differ from sample_list.txt"
        frac = frac.T
        frac.columns = frac.columns.str.split("/", n=1).str[0]  # "erythrocyte/erythroid progenitor" -> "erythrocyte"
        frac = frac.drop(columns=["r", "rmse"])  # nu-SVR fit statistics, not cell types
        cell_types = frac.columns.tolist()

        df = frac.join(metadata, how="inner")
        df["plasma_tubes_short_name"] = df["plasma_tubes"].replace(PLASMA_TUBE_DICT)
        if key == "Flomics_2":
            df["collection_center"] = df["collection_subcenter"]
        print(f"{key}: {len(df)} samples")
        df_dict[key] = df
    return df_dict, cell_types


# Order in which the MCTs are passed to the plots. The original notebook used list(set(...)),
# whose order changes between Python sessions and decides which leaf the dendrograms draw on
# top. This order reproduces the row order of all published heatmaps (Fig 4B, Suppl Fig 6).
MCT_ORDER = ["erythrocyte", "platelet", "monocyte", "b cell", "neutrophil", "macrophage"]


def select_mcts(df_dict, cell_types, min_abundance=0.1, min_fraction=0.1):
    """Major cell types: >= min_abundance in >= min_fraction of the samples of at least one dataset.

    Returned in MCT_ORDER; any cell type not listed there is appended, sorted.
    """
    mcts = set()
    for df in df_dict.values():
        mcts |= {c for c in cell_types if (df[c] >= min_abundance).mean() >= min_fraction}
    return [c for c in MCT_ORDER if c in mcts] + sorted(mcts - set(MCT_ORDER))


def add_plot_batches(df, labels):
    """dataset_batch -> display label; dataset_batch_custom splits Chen, Zhu and Flomics."""
    df = df.copy()
    df["dataset_batch"] = df["dataset_batch"].replace(labels)
    custom = df["dataset_batch"].copy()
    chen, zhu, flo = (df["dataset_batch"] == x for x in ("Chen", "Zhu", "Flomics"))
    peking = df["collection_center"] == "Peking University First Affiliated Hospital"
    custom[chen & peking] = "Chen (Peking University)"
    custom[chen & ~peking] = "Chen (Others)"
    custom[zhu & (df["phenotype"] == "Healthy")] = "Zhu (healthy)"
    custom[zhu & (df["phenotype"] == "Liver cancer")] = "Zhu (Cancer)"
    cordoba = df["collection_subcenter"] == "Cordoba"
    custom[flo & cordoba] = "Flomics (Cordoba)"
    custom[flo & ~cordoba] = "Flomics (Others)"
    df["dataset_batch_custom"] = custom
    df["broad_protocol_category_short"] = df["broad_protocol_category"].replace(BPC_SHORT)
    return df


# Annotation colours of the metadata tracks in the heatmaps.
HEATMAP_COLORS = {  # heatmap_color_dict of src/deconvolution_results.ipynb, verbatim
    "dataset_batch": {
        'Flomics': '#144d6b',
        'Block (2x75bp)': '#b3b3b3',
        'Block (2x150bp)': '#7d7a7a',
        'Decruyenaere': '#009E73',
        'Zhu': '#ffd633',
        'Chen': '#997a00',
        'Ngo': '#fa8072',
        'Roskams-Hieter (pilot)': '#944dff',
        'Roskams-Hieter (validation)': '#5b2e9e',
        'Moufarrej (Site 1)': '#CC79A7',
        'Moufarrej (Site 2)': '#CC79A7',
        'Sun': '#8a3d00',
        'Tao': '#0072B2',
        'Toden': '#800099',
        'Ibarra (buffy coat)': '#800000',
        'Ibarra (plasma, C)': '#A52A2A',
        'Ibarra (plasma, NC)': '#A52A2A',
        'Ibarra (serum)': '#B22522',
        'Chalasani': '#800040',
        'ENCODE\n(bulk tissue RNA-Seq)': '#006600',
        'Wei (cfDNA)': '#B32400',
        'Giráldez (phospho-RNA-seq)': '#B1CC71',
        'Giráldez (standard)': '#9cc43b',
        'Reggiardo (BioIVT)': '#F1085C',
        'Reggiardo (DLS)': '#F1085C',
        'Wang': '#FE8F42'
    },
    "collection_center_short_name":{
        'BBSSPA': '#A6CEE3',
        'BioIVT': '#1F78B4',
        'Discovery Life Sciences': '#B2DF8A',
        'GAPPS': '#33A02C',
        'GEMs (Indiana)': '#FB9A99',
        'IDIBAPS': '#E31A1C',
        'LPCH': '#FDBF6F',
        'NCC': '#FF7F00',
        'PKUFH': '#CAB2D6',
        'PUMC': '#6A3D9A',
        'SMMU': '#FFFF99',
        'SMMU (EHBH)': '#B15928',
        'San Diego Blood Bank': '#F0027F',
        'Scripps': '#666666',
        'Scripps (BMT)': '#BF5B17',
        'Scripps (CC)': '#6A3D9A',
        'Southwest Hospital': '#FCCDE5',
        'UCSD': '#BC80BD',
        'UKY': '#CCEBC5',
        'Unknown':"#333333",
        'WashU': '#FFED6F',
        'NA':'#CCCCCC'
    },
    'phenotype': {
        'Acute Myeloid Leukemia': '#8DD3C7',
        'Alzheimers disease': '#FFFFB3',
        'Chronic kidney failure EPO-treated': '#BEBADA',
        'Colorectal cancer': '#FFFFB3',
        'Diffuse large B-cell lymphoma': '#BEBADA',
        'Diverticulitis': '#FDB462',
        'Esophagus cancer': '#FB8072',
        'G-CSF-treated healthy donors': '#FCCDE5',
        'Healthy': '#1F78B4',
        'Liver cancer': '#80B1D3',
        'Lung cancer': '#FDB462',
        'Multiple myeloma': '#B3DE69',
        'Nonalcoholic fatty liver disease': '#C51B7D',
        'Nonalcoholic steatohepatitis': '#E99A6E',
        'Pancreatic cancer': '#FCCDE5',
        'Pre-cancerous condition: MGUS': '#98DF8A',
        'Pre-cancerous condition: cirrhosis': '#FF9896',
        'Pre-eclampsia': '#C49C94',
        'Primary mediastinal B-cell lymphoma': '#D9D9D9',
        'Stomach cancer': '#BC80BD',
        'healthy': '#1F78B4',
        'NA':'#CCCCCC'
    },
    "stage": {
        '0': '#F7FBFF',
        '1': '#DEEBF7',
        '1A': '#C6DBEF',
        '1B': '#9ECAE1',
        '2A': '#6BAED6',
        '2B': '#4292C6',
        '2C': '#2171B5',
        '3': '#08519C',
        '3A': '#084594',
        '3B': '#05317B',
        '3C': '#032152',
        '4': '#021029',
        'A': '#FEE0D2',
        'B': '#FCBBA1',
        'C': '#FC9272',
        'I-B': '#F8826F',
        'II-B': '#FB6A4A',
        'III': '#EF3B2C',
        'IV': '#CB181D',
        'No stage (healthy donor)': '#D9F0A3',
        'NA':'#CCCCCC'
    },
    'disease_stage': {
        'F0': '#FCFBFD',
        'F1': '#EFEDF5',
        'F2': '#DADAEB',
        'F3': '#BCBDDC',
        'F4': '#9E9AC8',
        'NA':'#CCCCCC'
    },
    'sex': {
        'female': '#FBB4AE',
        'male': '#B3CDE3',
        'NA':'#CCCCCC'
    },
    "biomaterial": {
        "buffy coat": "#E41A1C",
        "plasma": "#377EB8",
        "serum": "#4DAF4A",
        "tissue": "#984EA3",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    },
    "plasma_tubes_short_name": {
        "EDTA": "#E6A0C4",
        "EDTA/Streck": "#C6CDF7",
        "PAXgene blood ccfDNA": "#D8A499",
        "Vacutainer": "#7294D4",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    },
    "libraryselection": {
        "'Hybrid Selection'": "#A6CEE3",
        "cDNA": "#1F78B4",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    },
    "rna_extraction_kit_short_name": {
        "Maxwell (plasma and serum)": "#9E0142",
        "MOF": "#D53E4F",
        "Norgen": "#F46D43",
        "Norgen (slurry)": "#FDAE61",
        "Qiagen ExoRNeasy": "#FEE08B",
        "Qiagen miRNeasy": "#FFFFBF",
        "Qiagen miRNeasy Plasma": "#E6F598",
        "Qiagen QIAamp": "#ABDDA4",
        "Qiagen RNeasy Plasma": "#66C2A5",
        "Qiagen RNeasy Tissue": "#3288BD",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    },
    "library_prep_kit_short_name": {
        "Custom_Rozowsky": "#A6CEE3",
        "Custom_Wang": "#1F78B4",
        "Illumina TruSeq small RNA PNK-treated": "#B2DF8A",
        "Illumina TruSeq small RNA": "#33A02C",
        "SMART-Seq": "#FB9A99",
        "SMARTer Pico v1": "#E31A1C",
        "SMARTer Pico v2": "#FDBF6F",
        "SMARTer Pico v3": "#FF7F00",
        "Swift 2S kit ": "#CAB2D6",
        "Truseq Nano DNA HT": "#6A3D9A",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    },
    "dnase": {
        "Baseline-ZERO": "#7FC97F",
        "DNase I": "#FDC086",
        "DNase I + HL-dsDNase":"#BEAED4",
        "HL-dsDNase": "#FFFF99",
        "Turbo DNAse": "#386CB0",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    },
    "cdna_library_type": {
        "Reverse": "#8DD3C7",
        "Unstranded": "#BEBADA",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    },
    "read_length": {
        "1x50": "#EDF8E9",
        "1x75": "#BAE4B3",
        "2x75": "#74C476",
        "2x100": "#31A354",
        "2x150": "#006D2C",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    },
    "instrument": {
        'Illumina HiSeq 4000': '#A6CEE3',
        'Illumina NextSeq 500': '#1F78B4',
        'Illumina NovaSeq 6000': '#B2DF8A',
        'Ilumina HiSeq X Ten': '#33A02C',
        'NA': '#CCCCCC',
        'Unspecified': '#333333',
        'placeholder': '#FFF0F5',
        'None': '#B3B3B3'
    },
    "centrifugation_step_1": {
        "1000g": "#FFF5F0",
        "1500g": "#FEE0D2",
        "1600g": "#FCBBA1",
        "1900g": "#FC9272",
        "2000g": "#FB6A4A",
        "2500g": "#EF3B2C",
        "3000g": "#CB181D",
        "3400g": "#A50F15",
        "12000g": "#67000D",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    },
    "centrifugation_step_2": {
        "1940g": "#FFF5EB",
        "2500g": "#FEE6CE",
        "3000g": "#FDD0A2",
        "6000g": "#FDAE6B",
        "12000g": "#FD8D3C",
        "13000g": "#F16913",
        "15000g": "#D94801",
        "16000g": "#8C2D04",
        "NA": "#CCCCCC",
        "Unspecified": "#333333",
        "placeholder": "#FFF0F5",
        "None": "#B3B3B3"
    }
}

# Colours as used for Suppl Fig 6 (July 2025), reconstructed from the published panels:
# same as HEATMAP_COLORS except for the differences below. Values absent from this mapping
# get automatic colours (tab10, Set3, Pastel1, Dark2, Accent, in turn), as in the published panels.
HEATMAP_COLORS_20250702 = {k: dict(v) for k, v in HEATMAP_COLORS.items()
                           if k not in ("phenotype", "sex", "instrument", "stage", "libraryselection", "disease_stage")}
HEATMAP_COLORS_20250702["dataset_batch"].update({  # labels of src/dataset_mappings.json @ 7b4078a
    "Flomics 2": "#144d6b", "Ibarra (plasma, cancer)": "#A52A2A",
    "Ibarra (plasma, non-cancer)": "#A52A2A"})
HEATMAP_COLORS_20250702["rna_extraction_kit_short_name"].pop("Maxwell (plasma and serum)")
HEATMAP_COLORS_20250702["nucleic_acid_type"] = {"Total RNA": "#A1D76A"}
