"""Colour maps, orderings and the PCA/tSNE scatter used by dimred_perdataset.py.

Extracted from Flomics' internal bioinfo_utils so this repository runs on its own:
natural_sort_key from scripts/utils.py, the palettes and colour-map builders from
scripts/cfrna.py, and the dimred/plotting functions from scripts/ml.py. Only the
closure the dimred script actually reaches is kept; utils.py in particular pulls in
HTSeq at import time, which nothing here needs.

Two upstream patches are carried over: TSNE takes max_iter rather than the removed
n_iter, and the pdb.set_trace() calls are gone.
"""

import os
import re
import time

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA, FastICA
from sklearn.manifold import TSNE, MDS, SpectralEmbedding

random_state = 31415
np.random.seed(random_state)


# --- sorting -----------------------------------------------------------------

def natural_sort_key(ss):
    if not pd.isna(ss):
        return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', ss)]
    else:
        return [None]


# --- palettes ----------------------------------------------------------------

palette_version = "deep"
palette_version = "colorblind"
seaborn_palette = sns.color_palette(palette_version)
# Using: sns.color_palette("deep")
darkblue_col  = seaborn_palette[0]  # Dark blue   [0]
orange_col    = seaborn_palette[1]  # Orange      [1]
green_col     = seaborn_palette[2]  # Green       [2]
red_col       = seaborn_palette[3]  # Red         [3]
violet_col    = seaborn_palette[4]  # Violet      [4]
brown_col     = seaborn_palette[5]  # Brown       [5]
pink_col      = seaborn_palette[6]  # Pink        [6]
gray_col      = seaborn_palette[7]  # Gray        [7]
yellow_col    = seaborn_palette[8]  # Yellow      [8]
lightblue_col = seaborn_palette[9]  # Light blue  [9]

# Fix common labels
healthy_col = green_col
cancer_col  = red_col
ncd_col     = yellow_col
# Cancer-type
colorectal_col = darkblue_col
lung_col       = gray_col    
breast_col     = pink_col     
prostate_col   = lightblue_col
pancreatic_col = violet_col
# Non-cancer disease
ncd_colorectal_col = orange_col 
ncd_lung_col       = brown_col  
ncd_breast_col     = red_col    
ncd_prostate_col   = yellow_col 
ncd_pancreatic_col = 'gold'


def set_fontsize_screen():
    sns.set_style('whitegrid')
    mpl.rcParams['font.size'] = 20


# --- per-variable orderings and colour maps ----------------------------------

hue_order_map_dict = {
    'label':      [0, 1],
    'label_pred': [0, 1],
    'status':['healthy', 'cancer', 'non-cancer disease', 'nan'],
    'status_subtype':[
        'healthy', 'non-cancer disease',
        #'cancer',
        'colorectal', 'lung', 'breast', 'prostate', 'pancreatic', 'nan'],
    'phenotype_subtype':[
        "healthy", 
        "colorectal cancer", "lung cancer", "breast cancer", "pancreatic cancer", "prostate cancer",
        "NCD colorectal", "NCD lung", "NCD breast", "NCD pancreatic", "NCD prostate",
        'nan'
        ],
    'cancer_type':['colorectal', 'lung', 'breast', 'prostate', 'pancreatic', 'nan'],
    'stage_cancer_simple':['I', 'II', 'III', 'IV', 'healthy', '0', 'nan'],
    #'cancer_stage_simple':['I', 'II', 'III', 'IV', 'healthy'],
    'cancer_stage': ['healthy', 'early-stage', 'late-stage', 'nan'],
    'sex': ['male', 'female', 'Hombre', 'Mujer', 'nan'],
    'sequencing_batch':[
        # INSERM-CRC23 (Montpellier)
        'FL-SEQB-19',
        # LiquiDx Data
        'FL-SEQB-20','FL-SEQB-21','FL-SEQB-24','FL-SEQB-25','FL-SEQB-26','FL-SEQB-27','FL-SEQB-28',
        'FL-SEQB-29','FL-SEQB-31','FL-SEQB-32','FL-SEQB-35','FL-SEQB-36','FL-SEQB-37','FL-SEQB-38',
        'FL-SEQB-39','FL-SEQB-40','FL-SEQB-41','FL-SEQB-42','FL-SEQB-44','FL-SEQB-45','FL-SEQB-46',
        'FL-SEQB-47','FL-SEQB-50','FL-SEQB-53', 'FL-SEQB-54', 'FL-SEQB-56',
        'FL-SEQB-67', # Last batch!
        # Protocol Optimization
        'FL-SEQB-72',
        # LiquiDx Data - Other features
        'liquidx_results_00.part', 'liquidx_results_01.part', 'liquidx_results_02.part', 'liquidx_results_03.part', 
        'liquidx_results_04.part', 'liquidx_results_05.part', 'liquidx_results_06.part', 'liquidx_results_07.part',
        'liquidx_results_08.part', 'liquidx_results_09.part', 'liquidx_results_10.part', 'liquidx_results_11.part',
        'liquidx_results_12.part', 'liquidx_results_13.part', 'liquidx_results_14.part', 'liquidx_results_15.part',
        'liquidx_results_16.part', 'liquidx_results_17.part', 'liquidx_results_18.part', 'liquidx_results_19.part',
        # Non-LiquiDx Data (e.g. cf-meta)
        'Flomics_1', 'flomics_1', 'Flomics_2', 'flomics_2', 
        'block', 'chalasani', 'chen', 'decruyenaere', 'decru', 'giraldez', 'ibarra', 'moufarrej', 
        'ngo', 'reggiardo', 'roskams', 'rozowsky', 'sun', 'tao',
        'taowei', 'wei',
        'toden',
        'wang', 'wang_read_2',
        'zhu',
        # CRS
        "FL-SEQB-59",
        ],
    'collection_center': [
        # LiquidX ------------------------------------
        # Montpellier
        'INSERM-CRC23', 
        # IDIBAPS
        'Hospital Clínic de Barcelona',
        # Biobanco del Sistema Sanitario Público de Andalucía (BBSPA)
        "BBSSPA",
        # cfRNA-meta ---------------------------------
        # Lu Lab
        'Peking University First Affiliated Hospital',
        'Chinese Academy of Medical Sciences and Peking Union Medical College',
        'Southwest Hospital',
        'Eastern Hepatobiliary Surgery Hospital, Second Military Medical University',
        'National Center for Liver Cancer',
        'Second Military Medical University',
        'nan'
        ],
    'collection_subcenter': [
        # LiquidX ------------------------------------
        # Montpellier
        'INSERM-CRC23', 
        # IDIBAPS
        'Hospital Clínic de Barcelona',
        # Biobanco del Sistema Sanitario Público de Andalucía (BBSPA)
        'Cordoba', 'Granada', 'Granada Hospital', 'Sevilla', 'Jaen', 'Jaen CTTC', 'Malaga', 'Cadiz',
        # cfRNA-meta ---------------------------------
        # Lu Lab
        'Peking University First Affiliated Hospital',
        'Chinese Academy of Medical Sciences and Peking Union Medical College',
        'Southwest Hospital',
        'Eastern Hepatobiliary Surgery Hospital, Second Military Medical University',
        'National Center for Liver Cancer',
        'Second Military Medical University'
        #'nan'
        ],
    'inserm_crc23_report': [
        "BBSSPA single-spin",
        "BBSSPA double-spin",
        "IDIBAPS old (>4600 days)",
        "IDIBAPS young (<4600 days)",
        "INSERM-CRC23: Healthy",
        "INSERM-CRC23: Cancer",
        "nan",
    ],
    'inserm_crc23_report_anon': [
        "Biobank A single-spin",
        "Biobank A double-spin",
        "Biobank B old (>4600 days)",
        "Biobank B young (<4600 days)",
        "INSERM-CRC23: Healthy",
        "INSERM-CRC23: Cancer",
        "nan",
    ],
    'visual_inspection_color': [
        'Yellow',
        'Orange',
        'Red',
        'Unusual color',
        'nan'
    ],    
    'visual_inspection_cloudiness': [
        'Transparent',
        'Cloudy',
        'Turbid',
        'nan'
    ],
    #'resequenced_sample': ['True', 'False'], # Careful with type, boolean or string, need to be cosnsistent (!)
    #'resequence_version': ['old', 'new'],
    ## BBSSPA Table
    'centrifugation_step': [
        'Single_SpinTA', 'Single_Spin4C', 
        'Double_Spin', 'Double_SpinTA', 'Double_Spin4C',
        'Spin4C_m24h_Spin4C', 'Spin4C_mFrz_Spin4C', 
        'SpinNA_m24h_Spin4C', 'SpinNA_m24h_SpinNA', # v1
        'SpinTA_m24h_Spin4C',                       # v2
        '?',  # NA
        'nan' # Non-BBSPA
    ],
    'raw_processing_time': [
        '<1h', '<2h',
        "30'-2.5h", '1-2h', '1-3h', '1-4h',
        '?',  # NA
        'nan' # Non-BBSPA
    ],
    'anticoagulants': [
        'EDTA K2', 'EDTA K3',
        '?',  # NA
        'nan' # Non-BBSPA
    ],
}


discrete_color_map_dict = {
    'label': {
        0:healthy_col, 1:cancer_col
        },
    'label_pred': {
        0:healthy_col, 1:cancer_col
        },
    'status':{
        'healthy':healthy_col, 'cancer':cancer_col, 'non-cancer disease':ncd_col
        },
    'status_subtype':{
        'healthy':healthy_col, 'non-cancer disease':ncd_col,
        #'cancer':cancer_col, 
        "colorectal":colorectal_col, "lung":lung_col, "breast":breast_col, "pancreatic":pancreatic_col, "prostate":prostate_col,
        },
#     # more colorblind-friendly
#    'status_subtype':{
#         'healthy':seaborn_palette[2], 'non-cancer disease':seaborn_palette[8],
#         "colorectal":seaborn_palette[0], "lung":seaborn_palette[7], "breast":seaborn_palette[4], "pancreatic":seaborn_palette[5], "prostate":seaborn_palette[9],
#         #kk:seaborn_palette[ii] for ii, kk in enumerate(hue_order_map_dict['status_subtype'])
#         },
    'phenotype_subtype':{
        "healthy": green_col, 
        "colorectal cancer":colorectal_col, "lung cancer":lung_col, "breast cancer":breast_col, "pancreatic cancer":pancreatic_col, "prostate cancer":prostate_col,
        "NCD colorectal":ncd_colorectal_col, "NCD lung":ncd_lung_col, "NCD breast":ncd_breast_col, "NCD pancreatic":ncd_pancreatic_col, "NCD prostate":ncd_prostate_col
    },
    # https://www.fredhutch.org/en/news/center-news/2015/12/cancer-awareness-colors-cascade.html
    'cancer_type':  {
        "colorectal":colorectal_col, "lung":lung_col, "breast":breast_col, "pancreatic":pancreatic_col, "prostate":prostate_col,
        },
    'cancer_stage': {
        'healthy':healthy_col,
        'early-stage':yellow_col, 'late-stage':red_col
        },
    'stage_cancer_simple': {
        'healthy':healthy_col,
        '0':healthy_col,
        'I':'lightseagreen', 'II':'gold', 'III':'sandybrown', 'IV':red_col
        },
    # 'cancer_stage_simple': {
    #     'healthy':healthy_col,
    #     '0':healthy_col,
    #     'I':'lightseagreen', 'II':'gold', 'III':'sandybrown', 'IV':red_col
    #     },
    #'stage_cancer_1': {'0', 'I', 'IIa', 'IIIb', 'IV'},
    #'tnm_cancer_1': {'TisN0M0', 'T2N0M0', 'T3N0M0', 'T1N0M0', 'T3N2aM0', 'T3N1M0','T3N2bM1'}
    'sex': {
        'male':green_col, 'female':red_col,
        'Hombre': green_col, 'Mujer':red_col,
        },
    #'sequencing_batch': {'FL-SEQB-20':green_col, 'FL-SEQB-21':red_col, 'FL-SEQB-24':yellow_col},
    'dataset_name': {
        #'Flomics_1', 'Flomics_2', 
        "block":"#b3b3b3", "decruyenaere": "#009E73" , "zhu":"#ffd633", "chen":"#997a00", "ngo": "#fa8072", "roskams":"#944dff" ,
        "moufarrej":"#CC79A7", "sun":"#D55E00", "tao":"#0072B2", "toden":"#800099",  "ibarra":"#800000", "chalasani":"#800040",
        "rozowsky":"#006600", "taowei":"#B32400", "giraldez":"#B1CC71", "reggiardo":"#F1085C", "wang":"#FE8F42"
        }, 
    #'experiment_id': {'FL-EXP-139':green_col, 'FL-EXP-140':red_col,},
    #'library_prep_batch': {1.0:green_col, 2.0:red_col, 3.0:yellow_col, 4.0:'yellow'},
    # 'collection_center': {
    #     # LiquidX ------------------------------------
    #     # Montpellier
    #     'INSERM-CRC23':darkblue_col, 
    #     # IDIBAPS
    #     'Hospital Clínic de Barcelona':red_col,
    #     'HCB':red_col,
    #     # Biobanco del Sistema Sanitario Público de Andalucía
    #     "BBSSPA": green_col,
    #     #'LiquiDx':red_col
    #     # cfRNA-meta ---------------------------------
    #     # Lu Lab
    #     'Peking University First Affiliated Hospital': yellow_col,
    #     'Chinese Academy of Medical Sciences and Peking Union Medical College': yellow_col,
    #     'Southwest Hospital': yellow_col,
    #     'Eastern Hepatobiliary Surgery Hospital, Second Military Medical University': yellow_col,
    #     'National Center for Liver Cancer': yellow_col,
    #     'Second Military Medical University': yellow_col,
    #     },
    'collection_subcenter': {
        # LiquidX ------------------------------------
        # Montpellier
        'INSERM-CRC23':darkblue_col, 
        # IDIBAPS
        'Hospital Clínic de Barcelona':red_col,
        #'Cordoba':"brown", 'Granada':green_col, 'Sevilla':'gold', 'Jaen':'purple', 'Malaga':'magenta', 'Cadiz':'lightskyblue'
        'Cordoba':brown_col, 'Granada':green_col, 'Granada Hospital':green_col, 'Sevilla':yellow_col, 'Jaen':violet_col, 'Jaen CTTC':violet_col, 'Malaga':violet_col, 'Cadiz':lightblue_col,
        # cfRNA-meta ---------------------------------
        # Lu Lab
        'Peking University First Affiliated Hospital': yellow_col,
        'Chinese Academy of Medical Sciences and Peking Union Medical College': yellow_col,
        'Southwest Hospital': yellow_col,
        'Eastern Hepatobiliary Surgery Hospital, Second Military Medical University': yellow_col,
        'National Center for Liver Cancer': yellow_col,
        'Second Military Medical University': yellow_col,
        },
    'inserm_crc23_report': {
        "BBSSPA single-spin": violet_col,
        "BBSSPA double-spin": lightblue_col,
        "IDIBAPS old (>4600 days)": gray_col,
        "IDIBAPS young (<4600 days)": yellow_col,
        "INSERM-CRC23: Healthy": green_col,
        "INSERM-CRC23: Cancer": red_col,
    },
    'inserm_crc23_report_anon': {
        "Biobank A single-spin": yellow_col,
        "Biobank A double-spin": lightblue_col,
        "Biobank B old (>4600 days)": violet_col,
        "Biobank B young (<4600 days)": darkblue_col,
        "INSERM-CRC23: Healthy": green_col,
        "INSERM-CRC23: Cancer": red_col,
    },
    'user_1': {'LGA':green_col, 'LSA':red_col, 'SGG':yellow_col, 'BVN':pink_col},
    'user_2': {'LGA':green_col, 'LSA':red_col, 'SGG':yellow_col, 'BVN':pink_col},
    'dna_rna_isolation_lab_technician': {'LGA':green_col, 'LSA':red_col, 'SGG':yellow_col, 'BVN':pink_col},
    'library_prep_lab_technician':      {'LGA':green_col, 'LSA':red_col, 'SGG':yellow_col, 'BVN':pink_col},
    'visual_inspection_color': {
        'Yellow': yellow_col,
        'Orange': orange_col,
        'Red': red_col,
        'Unusual color':green_col,
    },    
    'visual_inspection_cloudiness': {
        'Transparent':lightblue_col,
        'Cloudy':gray_col,
        'Turbid':brown_col,
    },
    'hemolysis_comments': {
        'Yellow fluorescent': 'yellow',
        'Orange':            (0.9948327566320646, 0.874555940023068, 0.7530334486735871, 1.0),     # 'Orange_1
        'Orange ++':         (0.9921568627450981, 0.726797385620915, 0.49150326797385624, 1.0),    # 'Orange_2
        'Orange +++':        (0.9137254901960784, 0.3686274509803921, 0.050980392156862744, 1.0),  # 'Orange_3
        'Orange-Red':        (0.4980392156862745, 0.15294117647058825, 0.01568627450980392, 1.0),
        'Orange and cloudy': (0.9948327566320646, 0.874555940023068, 0.7530334486735871, 1.0),     # 'Orange_0 ("cloudy")
        'Red':               (0.9835755478662053, 0.4127950788158401, 0.28835063437139563, 1.0),   # 'Red_1'
        'Red ++':            (0.9344867358708189, 0.2286812764321415, 0.17139561707035755, 1.0),   # 'Red_2'
        'Red +++':           (0.7925720876585928, 0.09328719723183392, 0.11298731257208766, 1.0),  # 'Red_3'
        'Red++++':           'red',
        'Cloudy':            (0.8501191849288735, 0.8501191849288735, 0.8501191849288735, 1.0),     # 'Cloud_0', 
        'Cloudy +':          (0.7393771626297578, 0.7393771626297578, 0.7393771626297578, 1.0),     # 'Cloud_1'  
        'Cloudy ++':         (0.586082276047674, 0.586082276047674, 0.586082276047674, 1.0),        # 'Cloud_2'
        'Very cloudy':       (0.44844290657439445, 0.44844290657439445, 0.44844290657439445, 1.0),  # 'Cloud_3'
        },
    'resequenced_sample': {
        False:healthy_col, True:cancer_col
        },
    'resequence_version': {
        'new':healthy_col, 'old':cancer_col
        },
    ## BBSSPA Table
    'centrifugation_step': {
        'Single_SpinTA': yellow_col, 'Single_Spin4C': orange_col, 
        'Double_Spin': darkblue_col, 'Double_SpinTA': lightblue_col, 'Double_Spin4C':pink_col,
        'Spin4C_m24h_Spin4C': red_col, 'Spin4C_mFrz_Spin4C':green_col, 
        'SpinNA_m24h_Spin4C': gray_col, 'SpinNA_m24h_SpinNA':brown_col, # v1
        'SpinTA_m24h_Spin4C': gray_col,                                  # v2
        '?': 'k', # NA
    },
    'raw_processing_time': {
        '<1h': darkblue_col, '<2h': lightblue_col,
        "30'-2.5h":pink_col, '1-2h': green_col, '1-3h': orange_col, '1-4h': yellow_col  ,
        '?': 'k', # NA
    },
    'anticoagulants': {
        'EDTA K2': healthy_col, 'EDTA K3': cancer_col,
        '?': 'k', # NA
    },
}


# --- colour-map builders -----------------------------------------------------

def get_discrete_color_map(df, col_name, color_map_dict=None):

    # check if color scheme was manually defined
    if col_name not in color_map_dict.keys():
        color_map_dict=None
    
    col_map_df = df.set_index('sample_name')[[col_name]]
    # define data-driven color scheme
    if isinstance(color_map_dict, type(None)):
        label_map = (
            df[col_name]
            .drop_duplicates()
            .sort_values()
            .reset_index(drop=True).reset_index() # Keep new index
            .set_index(col_name)['index']         # Series
            .to_dict()
        )
        # cm.get_cmap was removed in matplotlib 3.9; plt.get_cmap takes the same lut
        cmap = plt.get_cmap('Dark2', len(label_map))
        col_map_df['color'] = col_map_df[col_name].map(lambda xx: cmap(label_map[xx]) if not pd.isnull(xx) else 'k')
        #cmap_dict = df.set_index('sample_name')[col_name].map(lambda xx: cm_viridis(label_map[xx])).to_dict()
    else:
        # try:
        #     assert all([(xx in color_map_dict[col_name].keys()) or (pd.isnull(xx)) for xx in col_map_df[col_name].unique()])
        # except AssertionError:
        #     print(f'{col_name=}')
        try:
            col_map_df['color'] = col_map_df[col_name].map(lambda xx: color_map_dict[col_name][xx] if not pd.isnull(xx) else 'k')
        except KeyError:
            print(f'\nError in: {col_name=}')
            missing_keys = [xx for xx in col_map_df[col_name].unique() if xx not in color_map_dict[col_name].keys()]
            print(f'\tMissing values: {missing_keys}')
            raise KeyError(f'{col_name}: no colour defined for {missing_keys}') from None
        #cmap_dict = df.set_index('sample_name')[col_name].map(lambda xx: color_map_dict[col_name][xx]).to_dict() 
    
    return col_map_df


def get_continous_color_map(df, col_name, discretize=False, n_clusters=3, debug=False):
    
    # only contains NA's
    empty_col = df[col_name].isnull().all()

    if not discretize or empty_col:
        cmap = plt.get_cmap('viridis')
        #cmap = mpl.colormaps['viridis'] 

        # Need to adjust for "outliers" that lead to lack of contrast
        #var_values = df[col_name].unique()
        var_values = df[col_name].tolist()
        mean_val = np.nanmean(var_values)
        std_val  = np.nanstd(var_values)
        # min_val/max_val or 'outlier' (4 * sd)
        n_sd = 3
        min_val = max(np.nanmin(var_values), mean_val - n_sd*std_val)
        max_val = min(np.nanmax(var_values), mean_val + n_sd*std_val)
        norm = plt.Normalize(
            min_val,
            max_val
        )

        cmap.set_under('magenta') # bottom arrow
        cmap.set_over('red')      # top arrow

        col_map_df = df.set_index('sample_name')[[col_name]]
        col_map_df['color'] = (
            col_map_df[col_name].map(lambda xx: cmap(norm(xx)) if not pd.isnull(xx) else 'k')
        )

        #return {kk: cmap(norm(vv)) for kk, vv in df.set_index('sample_name')[col_name].items()}

    else:
        # remove NA's
        kk_df = df[~df[col_name].isnull()]
        
        sample_names = kk_df.loc[:, 'sample_name'].tolist()
        xx = kk_df[col_name].to_numpy()
        # Reshape the sample for clustering (KMeans expects a 2D array)
        xx_reshaped = xx.reshape(-1, 1)
        
        # Apply K-means clustering to classify the data into 3 groups
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state).fit(xx_reshaped)
        class_labels = kmeans.labels_
        
        # Create a DataFrame to hold the sample and their cluster labels
        kk_df = pd.DataFrame({
            'sample_name': sample_names,
            'value': xx,
            'cluster_id': class_labels
        })
        
        cmap_dict = {'Low':green_col, 'Medium':yellow_col, 'High':red_col}
        # Map numerical cluster labels to meaningful names
        cluster_means = kk_df.groupby('cluster_id')['value'].mean().sort_values().reset_index()
        if (cluster_means.shape[0] == 3):
            cluster_names = ['Low', 'Medium', 'High']            
        elif (cluster_means.shape[0] == 2):
            cluster_names = ['Low', 'High']
        else:
            cluster_names = ['Medium']
        #display(cluster_means)
        cluster_means['cluster_name'] = cluster_names 
        sorted_clusters_map = cluster_means.set_index('cluster_id')['cluster_name'].to_dict()
        kk_df['cluster_name'] = kk_df['cluster_id'].map(sorted_clusters_map)


        # DEBUG  -----------------------------------------------------------------------------
        if debug:
            #print(kk_df)
            # Optional: Display the distribution of classes
            print(kk_df['cluster_name'].value_counts())
            
            # Plot the results
            pp_scatter = plt.scatter(
                kk_df.index, kk_df['value'], 
                c=kk_df['cluster_name'].map({kk:ii for ii, kk in enumerate(cluster_names)}), 
                cmap='viridis',
                label=kk_df['cluster_name']
            )
    
            plt.xlabel('Sample index')
            plt.ylabel('Value')
            plt.title('Sample Classification')
            # Cluster means not the boundaries (!)
            for ii in range(cluster_means.shape[0] - 1):
                cluster_boundary = cluster_means.value[ii:(ii+2)].mean()
                print(f'{cluster_boundary:.2f}')
                plt.axhline(y=cluster_boundary, color='r', linestyle='--')
            # Add labels to each point
            for ii, rr in kk_df.iterrows():
                plt.text(ii, rr['value'], rr.sample_name, fontsize=5, ha='right')
            # Add a legend
            handles, labels = pp_scatter.legend_elements()
            legend_labels = cluster_names
            plt.legend(handles, legend_labels, title=col_name)
            # Put the legend out of the figure (old)
            # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            # # sns.move_legend(ax_n_raw_reads_bp, "upper left", bbox_to_anchor=(1, 1)
            plt.show()
        # ------------------------------------------------------------------------------------
        col_map_df = (
            pd.merge(
                # add samples with NA's
                df[['sample_name']].reset_index(drop=True),
                kk_df,
                left_on='sample_name',
                right_on='sample_name',
                how='left'
                )
                .set_index('sample_name')[['cluster_name']]
                .rename(columns={'cluster_name': col_name})
        )
        col_map_df['color'] = col_map_df[col_name].map(lambda xx: cmap_dict[xx] if not pd.isnull(xx) else 'k')

        #return {kk: cmap_dict[vv] for kk, vv in df.set_index('sample_name')['cluster_name'].items()}

    return col_map_df


# --- dimensionality reduction and plotting -----------------------------------

def calculate_dimred_embedding(df, dimred_method, params_config=None, pc1=0, pc2=1, sample_ids=None, random_state=1):
    """
    Calculates the low-dimensional embedding using the specified dimensionality reduction method.
    """
    
    possible_dimred_methods = ['pca', 'tsne', 'ica', 'mds', 'se']
    if dimred_method not in possible_dimred_methods:
        raise ValueError(
            f"Unknown value for dimensionality reduction method argument: {dimred_method=}. "
            f"Valid arguments: {possible_dimred_methods}."
        )

    # Pre-processing  ----------------------------------------------------------------------

    if isinstance(sample_ids, type(None)):
        # include all numerical columns
        sample_ids = df.select_dtypes(include=[np.number]).columns
        
    # Filter for subset of samples: 'sample_id's should be present as columns in 'df'
    matrix = df[sample_ids].values.T
    
    # get number of components to compute
    # => from single-dimred or multiple-dimred (i.e. scatter-plot)
    n_components = max(pc1, pc2) + 1 if isinstance(pc1, int) else max(pc1 + pc2) + 1

    # Compute dimensionality reduction: 'dimred_method' -------------------------------------

    if dimred_method == 'pca':
        # Principal component analysis (PCA) 
        dimred_obj = PCA(
            n_components=n_components,
            random_state=random_state,
            whiten=True
        )
        X = dimred_obj.fit_transform(matrix)
        axis_name = 'PCA'

    elif dimred_method == 'tsne':
        # t-distributed stochastic neighbor embedding (t-SNE) 
        if isinstance(params_config, type(None)):
            params_config = {'perplexity': 5}
            print(f"\nWarning! No perplexity parameter was given. Default: perplexity={params_config['perplexity']}\n")
        dimred_obj = TSNE(
            n_components=n_components,
            random_state=random_state,
            max_iter=5000,
            init = "pca",
            early_exaggeration=6,
            **params_config
        )
        X = dimred_obj.fit_transform(matrix)
        axis_name = 't-SNE'

    elif dimred_method == 'ica':
        # Independent component analysis (ICA)
        dimred_obj = FastICA(
            n_components=n_components,
            random_state=random_state,
            whiten='unit-variance'
        )
        X = dimred_obj.fit_transform(matrix)
        axis_name = 'ICA'

    elif dimred_method == 'se':
        # Spectral Embedding (SE)
        dimred_obj = SpectralEmbedding(
            n_components=n_components,
            random_state=random_state
        )
        X = dimred_obj.fit_transform(matrix)
        axis_name = 'SE'

    else:  
        # Multidimensional scaling (MDS)
        dimred_obj = MDS(
            n_components=n_components,
            random_state=random_state
        )
        X = dimred_obj.fit_transform(matrix)
        axis_name = 'MDS'

    return X, dimred_obj, axis_name, sample_ids


def prepare_color_scheme(sample_ids, sample_names, color_map=None, marker_map=None, color_by=None):

    n_samples = len(sample_ids)

    # init Dataframe: sample-level
    sample_label_df = pd.DataFrame({
        'sample_id':sample_ids,
        'sample_name':sample_names
    })

    # if no 'color_map' is given - use 'sample_name's as 'label' to color samples 
    if isinstance(color_map, type(None)):

        if not isinstance(color_by, type(None)):
            print(f"\nIgnoring {color_by=} - Missing `color_map` argument!") 
        color_by = 'sample_name'
        
        # map between 'sample_id' and 'label' (i.e. 'sample_name's)
        sample_label_df['label'] = sample_label_df['sample_name']

        # init Dataframe: label-level
        label_color_df = pd.DataFrame({
            #  'sample_name's (i.e. 'label's) might not be unique (e.g. when ploting 'genes' instead of 'samples')
            'label': sorted(sample_label_df['label'].unique(), key=natural_sort_key)
            })
        n_labels = label_color_df.shape[0]

        # default 'color_map' for 'label's: viridis
        print(f"\nUsing default `color_map` (viridis) with {color_by=}.\n")
        cmap = mpl.colormaps['viridis'].resampled(n_labels)
        
        # map between category 'label' (i.e. 'sample_name's) and 'color'
        label_color_df['color'] = [cmap(ii) for ii in range(n_labels)]

    else:

        if isinstance(color_by, type(None)):
            color_by = color_map.columns[0]
            print(f"\nDeterminining which {color_by=} to use from `color_map`, default to first column.\n")

        # Color map -----------------------------------------------------------

        # map between 'sample_id' and 'label' (i.e. sample metadata values) from color_map
        sample_label_df = sample_label_df.merge(
            color_map[color_by].reset_index().rename(columns={'sample_name':'sample_id', color_by:'label'}),
            on=['sample_id'],
            how='left'
        )

        # init Dataframe: label-level
        label_color_df = (
            color_map
            .rename(columns={color_by:'label'})
            [['label', 'color']]
            .drop_duplicates()
            .copy()
            # sort for consistent alphabetical ordering of labels within each color (for 'marker' assignment)
            .sort_values(by=['color', 'label'])
        )

    # Marker map -----------------------------------------------------------
    
    # Only necessary for 'discrete' color schemes
    if (color_by == 'sample_name') or (np.issubdtype(label_color_df.dtypes['label'], np.number)):
        # assign marker index within each color group
        # => for consistency, in this case 'color's are unique and thus 'marker's are the default: 'o'
        label_color_df['marker_idx'] = 0
        # map index to marker (default: 'o')
        label_color_df['marker'] = 'o'

    else:
        if isinstance(marker_map, type(None)):
            # assign marker index within each color group
            label_color_df['marker_idx'] = label_color_df.groupby('color').cumcount()
            # map index to marker
            marker_types = ['o', '^', 's', 'D', 'P', 'X', '*', 'v']
            label_color_df['marker'] = label_color_df['marker_idx'].map(lambda x: marker_types[x])
        else:
            # Manual marker assignment
            label_color_df['marker'] = label_color_df['label'].map(marker_map)


    # add color scheme to sample-level DataFrame:
    sample_label_df = sample_label_df.merge(
        label_color_df, 
        on=['label'],
        how='left'
    )
    assert sample_label_df.shape[0] == n_samples

    return sample_label_df, label_color_df, color_by


def plot_dimred_embedding(X, dimred_obj, pc1=0, pc2=1, sample_ids=None, sample_names=None, annotate_plot=False, full_matrix=False, color_by=None, color_map=None, marker_map=None, hue_order=None, marker_size=60, show_cbar=False, axis_name='', title='', outpath=None, **kwargs):
    """
    Visualize the low-dimensional embedding.
    """
    
    # Set font: Arial
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['mathtext.bf'] = 'Arial:bold'
    # Export .svg
    plt.rcParams["svg.fonttype"] = "none"

    # Default
    layout_config = {
        'subplot_size': 10, # inches
        # Font size ----
        # Figure
        'axis_label': 20,
        'axis_label_pad': 20,
        'axis_ticks': 15,
        'axis_ticks_width': 1,
        'axis_ticks_length': 4,
        'axis_ticks_pad':    4,
        'subtitle':   20,
        'spine_linewidth': 3,
        # Legend
        'legend_fontsize': 10,
        'legend_title': 10,
        # Color-bar
        'cbar_ticks': 5,
        'cbar_ticks_width': 1,
        'cbar_ticks_length': 4,
        'cbar_ticks_pad':    4,
        'cbar_label': 10,
        'cbar_linewidth': 1,
        # Sample
        'annotation': 10,
    }

    font_size = 5
    # cfRNA-meta
    layout_config = {
        'subplot_size': 1.55, # inches
        # Font size ----
        # Figure
        'axis_label': font_size,
        'axis_label_pad': 2,
        'axis_ticks': font_size,
        'axis_ticks_width': 0.25,
        'axis_ticks_length': 1,
        'axis_ticks_pad':    1.,
        'subtitle':   font_size,
        'spine_linewidth': 0.5,
        # Legend
        'legend_fontsize': font_size,
        'legend_title': font_size,
        # Color-bar
        'cbar_ticks': font_size,
        'cbar_ticks_width': 0.25,
        'cbar_ticks_length': 1,
        'cbar_ticks_pad':    1.,
        'cbar_label': font_size,
        'cbar_linewidth': 0.25,
        # Sample
        'annotation': font_size,
    }

    # use 'sample_id's to annotate samples if 'sample_name's not given
    if isinstance(sample_names, type(None)):
        if isinstance(sample_ids, type(None)):
            sample_ids = [f"sample_{ii}" for ii in range(X.shape[0])]
        sample_names = sample_ids
    # check length consistent
    assert len(sample_ids) == len(sample_names)

    # Prepare Data: Coloring scheme ---------------------------------------------------------

    # Both at sample-level and label-level (for legend generation)
    sample_label_df, label_color_df, color_by = prepare_color_scheme(
        sample_ids, sample_names, 
        color_map=color_map, 
        marker_map=marker_map,
        color_by=color_by
    )
    # Ensure correct order (same as in X, based on sample_id's)
    sample_label_df = sample_label_df.set_index('sample_id').loc[sample_ids]

    # number of elements in legend
    n_labels = label_color_df.shape[0]
    
    # Plot variable: categorical or numerical
    numeric_values = np.issubdtype(label_color_df.dtypes['label'], np.number)
    unique_markers = label_color_df['marker'].unique()

    # Plot: 2D Embedding ---------------------------------------------------------

    # Default: use legend (instead of color bar)
    # => except if (numerical) and (too many elements in legend)
    if (numeric_values) and (n_labels > 20):
        try:
            # check if elements are numerical
            # => True if string is a number
            show_cbar = True
            print(f'Number of elements in legend is too big: {n_labels=} (legend will be hidden).')
        except:
            pass

    # Figure layout ----------------------------------

    # if single-dimred
    if isinstance(pc1, int):
        pc1, pc2 = [pc1], [pc2]
    assert len(pc1) == len(pc2)
    n_plots = len(pc1)

    if full_matrix:
        n_cols = max(pc1) + 1
        n_rows = max(pc2) + 1
    else:
        n_cols = min(n_plots, 3)
        n_rows = (n_plots - 1) // 3 + 1

    # before padding
    fig_size_base = (n_cols * layout_config['subplot_size'], n_rows * layout_config['subplot_size'])

    x_padding = layout_config['subplot_size']
    y_padding = 0
    # add some padding: adjust size of plot with legend/cbar
    fig_size = (fig_size_base[0] + x_padding, fig_size_base[1] + y_padding)

    # Init figure
    if (n_rows == 1) and (n_cols == 1):
        single_plot=True
        # fig = plt.figure(figsize=fig_size)
        # gs  = gridspec.GridSpec(
        #     nrows=1,
        #     ncols=2, 
        #     width_ratios=[fig_size_base[0], x_padding], 
        #     #wspace=0.05,
        #     wspace=0,
        #     hspace=0,
        # )
        # axs = [fig.add_subplot(gg) for gg in gs]
        # Use axes of absolute size in a large enough figure canvas
        aspect     = 1
        # Axis
        ax_width   = layout_config['subplot_size']
        ax_height  = ax_width / aspect
        # emmbed
        #x_pad = 0.2
        fig_width  = ax_width  * 2
        fig_height = ax_height * 1.2
        fig = plt.figure(figsize=(fig_width, fig_height))
        # Figure coordinates: [0, 1]
        axs = [
            fig.add_axes([
                (ax_width * ii) / fig_width,         # left
                0.1,                    # bottom
                ax_width / fig_width,  # width
                ax_height / fig_height  # height
            ]) for ii in range(2)
            ]
    else:
        single_plot=False
        fig, axs = plt.subplots(
            nrows=n_rows,
            ncols=n_cols, 
            figsize=fig_size,
            squeeze=False
        )
        # Flatten for looping
        axs = axs.flatten()

    # --------
    # Dim-red: PC_{ii}  vs PC_{jj}
    # --------
    
    # loop over axis
    for ii, ax in enumerate(axs):

        #ax.set_aspect('equal')  # Set equal aspect ratio
        
        # Plot layout ------------------------------------------------

        # Avoid plotting diagonal and empty figures
        if (ii >= len(pc1)) or (pc1[ii] == pc2[ii]):
            ax.axis('off')
            continue
        # Upper diagonal only
        elif (full_matrix) and (pc1[ii] >= pc2[ii]):
            ax.axis('off')
            continue

        ax.spines['top'].set_linewidth(layout_config['spine_linewidth'])
        ax.spines['right'].set_linewidth(layout_config['spine_linewidth'])
        ax.spines['bottom'].set_linewidth(layout_config['spine_linewidth'])
        ax.spines['left'].set_linewidth(layout_config['spine_linewidth'])

        # Axis-Labels
        if isinstance(dimred_obj, PCA):
            xlabel = f"Principal Component {pc1[ii] + 1}"
            ylabel = f"Principal Component {pc2[ii] + 1}"
        else:
            xlabel = f"{axis_name} - Component {pc1[ii] + 1}" if axis_name else f"Component {pc1[ii] + 1}"
            ylabel = f"{axis_name} - Component {pc2[ii] + 1}" if axis_name else f"Component {pc2[ii] + 1}"
        if hasattr(dimred_obj, 'explained_variance_ratio_'):
            xlabel += f" ({100 * dimred_obj.explained_variance_ratio_[pc1[ii]]:.2f}%)"
            ylabel += f" ({100 * dimred_obj.explained_variance_ratio_[pc2[ii]]:.2f}%)"

        ax.set_ylabel(ylabel, fontsize=layout_config['axis_label'], labelpad=layout_config['axis_label_pad'])
        ax.set_xlabel(xlabel, fontsize=layout_config['axis_label'], labelpad=layout_config['axis_label_pad'])

        ax.tick_params(
            labelsize=layout_config['axis_ticks'], 
            length=layout_config['axis_ticks_length'], width=layout_config['axis_ticks_width'], pad=layout_config['axis_ticks_pad']
        )
        
        # plt.xticks(fontsize=layout_config['axis_ticks'], rotation=45)
        # plt.yticks(fontsize=layout_config['axis_ticks'])

        if title:
            #plt.suptitle(title, fontsize=layout_config['subtitle'], y=1)
            plt.suptitle(title, fontsize=layout_config['subtitle'], y=0.92)
            
        # Plot the points and annotate ------------------------------

        # x = X[:, pc1[ii]]
        # y = X[:, pc2[ii]]
        
        sample_label_df['x'] = X[:, pc1[ii]]
        sample_label_df['y'] = X[:, pc2[ii]]

        # Annotate points with 'sample_names'
        # => Loop over samples: very slow!
        if annotate_plot:
            #for idx, ss in enumerate(sample_label_df['sample_name']):
            for rr in sample_label_df.itertuples(index=False):
                ax.annotate(
                    #sample_names[idx], (x[idx], y[idx]),
                    rr.sample_name, (rr['x'], rr['y']),
                    fontsize=layout_config['annotation'],
                    textcoords="offset points",
                    xytext=(10, 10), ha='center', style='italic'
                )

        if len(unique_markers) == 1:
            # Plot all samples
            ax.scatter(
                sample_label_df['x'].values, sample_label_df['y'].values,
                s=marker_size,
                c=sample_label_df['color'].values,
                linewidths=marker_size*0.05,
                label=sample_label_df['label'].values,
                marker=unique_markers[0],
                **kwargs
            )
        else:
            # Plot each (subgroup) of samples with identical parameters (i.e. 'marker') individually
            # matplotlib.pyplot.scatter() does not support passing an array of marker styles
            # => it expects a single marker style string (like 'o', 's', etc.) for the whole call.
            for mm, group in sample_label_df.groupby('marker'):
                ax.scatter(
                    group['x'].values, group['y'].values,
                    s=marker_size,
                    linewidths=marker_size*0.05,
                    c=group['color'].values,
                    label=group['label'].values,
                    marker=mm,
                    **kwargs
                )

    # Legend or Colorbar -----------------------------------------------------

    # avoid legend redundancy when multiple plots
    if n_plots > 1:
        if full_matrix:
            # use first since diagonal is empty
            ax_legend = axs[0]
        else:
            # use second plot instead of last, since last will be empty
            ax_legend = axs[1]
    else:
        ax_legend = axs[-1]
        
    if (not show_cbar):
        # # Create a legend with unique labels
        # unique_labels = {str(ll):cc for ll, cc in label_color_map.items()}
        # if not isinstance(hue_order, type(None)):
        #     try:
        #         assert all([ii in hue_order for ii in unique_labels.keys()])
        #     except AssertionError:
        #         pass
        #     assert all([ii in hue_order for ii in unique_labels.keys()])
        #     # re-order labels
        #     unique_labels = {ii:unique_labels[ii] for ii in hue_order if ii in unique_labels.keys()}
        # handles = [
        #     plt.Line2D(
        #         [0], [0], marker='o', color='w', markerfacecolor=cc, markersize=10, label=ll
        #         ) for ll, cc in unique_labels.items()
        # ]

        # Create a legend with unique labels
        unique_labels_df = label_color_df.astype({'label':str}).set_index('label')

        if not isinstance(hue_order, type(None)):
            try:
                assert all([ii in hue_order for ii in unique_labels_df.index])
            except AssertionError:
                pass
            #assert all([ii in hue_order for ii in unique_labels_df.index])
            # re-order labels
            unique_labels_df = unique_labels_df.loc[[ii for ii in hue_order if ii in unique_labels_df.index]]

        handles = []
        for rr in unique_labels_df.itertuples(index=True):
            handles.append(
                plt.Line2D(
                    [0], [0], marker=rr.marker, color='w', markerfacecolor=rr.color, markersize=marker_size, label=rr.Index
                )
            )

        # if too many labels - move legend outside plot
        #if len(handles) < 10 or full_matrix:
        if full_matrix:
            ax_legend.legend(
                handles=handles,
                title=color_by,
                loc='best',
                fontsize=layout_config['legend_fontsize'],
                title_fontsize=layout_config['legend_title'],
                )
        else:
            # TODO: Currently added it into a separate subplot (!)              
            ax_legend.legend(
                handles=handles,
                title=color_by,
                loc='center left',
                #bbox_to_anchor=(1.0, 0.5), borderaxespad=0.5,
                fontsize=layout_config['legend_fontsize'],
                title_fontsize=layout_config['legend_title'],
            )

    elif show_cbar:

        #ax.get_legend().remove()

        # Default: 'viridis'
        cmap = mpl.colormaps['viridis']  

        # Need to adjust for "outliers" that lead to lack of contrast
        # => TODO: Ideally the min, max values should be the same than when defining the color scheme (!)
        label_values = sample_label_df['label'].values
        
        # outlier definition
        mean_val = np.nanmean(label_values)
        std_val  = np.nanstd(label_values)
        # min_val/max_val or 'outlier' (4 * sd)
        n_sd = 3
        min_val = max(np.nanmin(label_values), mean_val - n_sd*std_val)
        max_val = min(np.nanmax(label_values), mean_val + n_sd*std_val)
        norm = plt.Normalize(
            min_val,
            max_val
        )
        # define outlier color
        cmap.set_under('magenta') # bottom arrow
        cmap.set_over('red')      # top arrow

        # add colorbar
        cbar = fig.colorbar(
            cm.ScalarMappable(norm=norm, cmap=cmap),
            ax=ax_legend,
            location='left' if single_plot else 'right',
            shrink=0.7,
            extend='both'
        )
        #cbar.set_label(f'{color_by}', fontsize=layout_config['cbar_label'], loc='center', transform_rotates_text=True)
        #cbar.ax.set_xlabel(f'{color_by}', fontsize=layout_config['cbar_label'], loc='center')
        cbar.ax.set_ylabel(
            f'{color_by}', fontsize=layout_config['cbar_label'],
            loc='center', 
            #rotation=270, labelpad=0.1,
            transform_rotates_text=True
        )
        cbar.ax.yaxis.set_label_position('right')
        cbar.ax.tick_params(
            labelsize=layout_config['cbar_ticks'],
            length=layout_config['cbar_ticks_length'], width=layout_config['cbar_ticks_width'], pad=layout_config['cbar_ticks_pad']
        )
        cbar.outline.set_linewidth(layout_config['cbar_linewidth'],)

    # Store figure to file --------------------------------

    if not isinstance(outpath, type(None)):
        if isinstance(outpath, str):
            outpath = [outpath]
        for pp_outpath in outpath:
            plt.savefig(
                pp_outpath,
                bbox_inches='tight',
                dpi=600
            )
    plt.close()

    return dimred_obj, fig


def new_scatter_dimred_plot(df, dimred_method, params_config=None, pc1=0, pc2=1, sample_ids=None, sample_names=None, annotate_plot=False, full_matrix=False, random_state=1, color_by=None, color_map_dict=None, marker_map_dict=None, hue_order_dict=None, marker_size=60, show_cbar=False, show_title=True, plot_sx='', plot_ext='.pdf', out_dir=None, **kwargs):    
    # Embedding
    start_time = time.time()
    X, dimred_obj, axis_name, sample_ids = calculate_dimred_embedding(
        df=df,
        dimred_method=dimred_method,
        params_config=params_config,
        pc1=pc1,
        pc2=pc2,
        sample_ids=sample_ids,
        random_state=random_state 
    )
    print(f'\nCalculate dimred embedding: {time.time() - start_time:.2f} s.')

    if isinstance(hue_order_dict, type(None)):
        hue_order_dict = {}    
    if isinstance(marker_map_dict, type(None)):
        marker_map_dict = {}
    
    # Visualization - loop over coloring vars
    for kk in color_by:

        # Plotting extensions
        if not isinstance(out_dir, type(None)):
            if isinstance(plot_ext, str):
                plot_ext = [plot_ext]
            dimred_pp = []
            for pp_ext in plot_ext:
                dimred_pp.append(os.path.join(out_dir, f'2d-{dimred_method}_{kk}{plot_sx}{pp_ext}'))
            print(dimred_pp)
        else:
            dimred_pp = None

        start_time = time.time()
        plot_dimred_embedding(
            X, dimred_obj,
            pc1=pc1 if isinstance(pc1, int) else list(pc1),
            pc2=pc2 if isinstance(pc2, int) else list(pc2),
            sample_ids=sample_ids,
            sample_names=sample_names,
            annotate_plot=annotate_plot,
            full_matrix=full_matrix,
            color_by=kk,
            # DataFrame with 'names' as index, 'labels' column (e.g. 'status'), and corresponding 'colors'
            color_map=color_map_dict[kk],
            # Dictionary
            marker_map=marker_map_dict[kk] if kk in marker_map_dict.keys() else None,
            hue_order=hue_order_dict[kk] if kk in hue_order_dict.keys() else None,
            marker_size=marker_size, # Default: 60
            show_cbar=show_cbar,
            axis_name=axis_name,
            title=f'{axis_name} - {kk} ' if show_title else '',
            outpath=dimred_pp,
            **kwargs
        )
        print(f'\t- Visualize dimred embedding ({kk}): {time.time() - start_time:.2f} s.')

    # Store components (?)
    print('Done.\n')

    return dimred_obj, X
