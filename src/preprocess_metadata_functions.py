import os
import pandas as pd
import numpy as np
import re
import pdb
from sra_columns_mapping import rename_columns_and_values
pd.set_option('display.max_columns', None)

import GEOparse


def simplify_column_names(cols):
    if type(cols) is not pd.Series:
        cols = pd.Series(cols)
    cols = (cols
            .str.lower()
            .str.replace(r'\s', r'_', regex=True)
            .str.replace(r'[()]', r'', regex=True)
           ).to_list()
    return cols

reserved_vars_samplesheet = [
    'sample_idx_sequencing_batch', 'sample_display_name', 
    'fastq_1', 'fastq_2', 'out_dir',
    'resequenced_sample'
]


dataset_column_list = [
    'dataset_short_name',
    'SRA link',
    'Biomaterial',
    'Nucleic acid type',
    'Library selection',
    'Assay name',
    'Plasma volume',
    'Plasma tubes',
    'RNA extraction kit',
    'RNA extraction kit (short name)',
    'DNAse',
    'Library prep kit',
    'Library prep kit (short name)',
    'cDNA library type',
    ]
dataset_column_list = simplify_column_names(dataset_column_list)


def merge_sample_with_dataset_metadata(df, dataset_metadata, keep_sample_cols=[]):
    # Before merging, verify that the columns we have defined manually
    # are not present in the dataset_metadata table.
    df.columns = simplify_column_names(df.columns)
    cols_common = (set(df.columns).intersection(set(dataset_column_list))
                   - {'dataset_short_name'}
                   - set(keep_sample_cols))
    if len(cols_common) > 0:
        raise RuntimeError(f"The following columns are present in the sample-level metadata dataframe and also in the dataset-level dataframe: {cols_common}")
    df = df.merge(
        dataset_metadata[list(set(dataset_column_list) - set(keep_sample_cols))],on='dataset_short_name', how='left')
    return df


def preprocess_chen(dataset_metadata):
    print("### Dataset: chen")
    csv_path = "../sra_metadata/chen_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "chen" 
    df["dataset_short_name"] = "chen"
    df["dataset_batch"] = "chen"
    df["read_length"] = "2x150"
    df["centrifugation_step_1"] = "Unspecified"
    df["centrifugation_step_2"] = "Unspecified"

    # Exclude the two E. coli samples and the brain tissue sample
    df = df[~(df['sample_name'].isin([
        'ET_L2',
        'EH_L2',
        'NC_L2',
    ]))]

    # Parse the GEO series matrix file, which contains the
    # the GEO/GSM ids and the collection center (hospital).
    with open("../sra_metadata/chen_GSE174302_series_matrix.txt") as f:
        GEO_matrix = f.readlines()

    line = [line for line in GEO_matrix if re.search(r'!Series_sample_id', line)][0]
    s = re.search(r'^!Series_sample_id\t"(.+)"\s*', line).group(1).strip()
    GSM_ids = s.split()

    line = [line for line in GEO_matrix if re.search(r'!Sample_description', line)][1]
    s = re.search(r'^!Sample_description\t"(.+)"\s*', line).group(1).strip()
    sample_descriptions = s.split('"\t"')
    collection_centers = [re.search(r'Provided by (.+?)\.?;', s).group(1) for s in sample_descriptions]

    if len(collection_centers) != len(GSM_ids):
        raise RuntimeError("Parsing the series matrix file, found different number of sample ids and collection_centers.")
    gsm_id_collection_center = pd.DataFrame(np.transpose([GSM_ids, collection_centers]),
                                            columns=['GSM_id', 'collection_center'])
    df = (df
        .merge(gsm_id_collection_center, left_on='geo_accession_exp',
                right_on='GSM_id', how='left')
    )
    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/chen_metadata_preprocessed.csv", index=False)
    return df

def preprocess_zhu(dataset_metadata):
    print("### Dataset: zhu")
    csv_path = "../sra_metadata/zhu_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "zhu" 
    df["dataset_short_name"] = "zhu"
    df["dataset_batch"] = "zhu"
    df["read_length"] = "2x150"
    df["centrifugation_step_1"] = "Unspecified"
    df["centrifugation_step_2"] = "Unspecified"

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/zhu_metadata_preprocessed.csv", index=False)
    return df

def preprocess_roskams(dataset_metadata):
    print("### Dataset: roskams")
    csv_path = "../sra_metadata/roskams_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)    

    df["sequencing_batch"] = "roskams" 
    df["dataset_short_name"] = "roskams"
    df["dataset_batch"] = np.where(df["cohort"] == "pilot", "roskams_pilot", "roskams_validation")
    df["read_length"] = np.where(df["dataset_batch"] == "roskams_pilot", "2x100", "2x150")
    df["centrifugation_step_1"] = "1000g"
    df["centrifugation_step_2"] = "15000g" 

    # Phenotype is described in the source_name column
    df['source_name'].unique()
    df['phenotype'] = df['source_name'].replace({
        'Human non-cancer donor plasma':'Healthy',
        'Human multiple myeloma plasma':'Multiple myeloma',
        'Human MGUS plasma':'Pre-cancerous condition: MGUS',
        'Human liver cancer plasma':'Liver cancer',
        'Human liver cirrhosis plasma':'Pre-cancerous condition: cirrhosis',
    })

    # The supplementary table 2 contains sample-level information about the RNA
    # extraction and library preparation batches
    supp_table = pd.read_csv("../sra_metadata/roskams_supp_table_2.tsv", sep='\t')
    supp_table['RNA extraction batch'] = (supp_table['RNA Extraction']
                                        .str.extract(r'batch\s*(\d+)'))
    supp_table['Library preparation batch'] = (supp_table['Library Preparation']
                                            .str.extract(r'batch\s*(\d+)'))
    # Parse the GEO series matrix file, which contains the mapping between
    # the GEO/GSM ids and the sample names in the study (PP02, etc)
    with open("../sra_metadata/roskams_GSE182824_series_matrix.txt") as f:
        GEO_matrix = f.readlines()

    line = [line for line in GEO_matrix if re.search(r'!Series_sample_id', line)][0]
    s = re.search(r'^!Series_sample_id\t"(.+)"\s*', line).group(1).strip()
    GSM_ids = s.split()

    line = [line for line in GEO_matrix if re.search(r'!Sample_title', line)][0]
    s = re.search(r'^!Sample_title\t"(.+)"\s*', line).group(1).strip()
    sample_names = s.split('"\t"')
    sample_names
    if len(sample_names) != len(GSM_ids):
        raise RuntimeError("Parsing the series matrix file, found different number of sample ids and sample names.")
    sample_id_mapping = pd.DataFrame(np.transpose([GSM_ids, sample_names]),
                                    columns=['GSM_id', 'Sample_id'])
    # Merge the metadata table with the supp table using the sample id mapping
    df = (df
        .merge(sample_id_mapping, left_on='geo_accession_exp',
            right_on='GSM_id', how='outer')
        .merge(supp_table[['SeqID', 'Library preparation batch', 'RNA extraction batch']]
                .rename(columns={'SeqID':'Sample_id'}), on='Sample_id', how='outer')
    )

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/roskams_metadata_preprocessed.csv", index=False)
    return df

def preprocess_ngo(dataset_metadata):
    print("### Dataset: ngo")
    csv_path = "../sra_metadata/ngo_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "ngo" 
    df["dataset_short_name"] = "ngo"
    df["dataset_batch"] = "ngo"
    df["read_length"] = "2x75"
    df["centrifugation_step_1"] = "Unspecified"
    df["centrifugation_step_2"] = "Unspecified"


    df['phenotype'] = df['isolate'].replace({
        "Pregnant adult female who delivered full-term":"Healthy pregnant women",
        "Pregnant adult female who delivered spontaneously preterm":"Healthy pregnant women who delivered preterm"
    })
    
    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/ngo_metadata_preprocessed.csv", index=False)
    return df

def preprocess_ibarra(dataset_metadata):
    print("### Dataset: ibarra")
    csv_path = "../sra_metadata/ibarra_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "ibarra" 
    df["dataset_short_name"] = "ibarra"
    df["biomaterial"] = df["tissue"].str.lower()
    df["read_length"] = "2x75"

    # Guess phenotype
    def parse_phenotype(s):
        m = re.search(r'MM', s)
        if m:
            return pd.Series({"phenotype": "Multiple myeloma",
                            "collection_center": "Scripps Bone Marrow Transplant Center"})
        m = re.search(r'EPO', s)
        if m:
            return pd.Series({"phenotype": "Chronic kidney failure EPO-treated",
                            "collection_center": "Scripps Clinic Cancer Center"})
        m = re.search(r'AML', s)
        if m:
            return pd.Series({"phenotype": "Acute Myeloid Leukemia",
                            "collection_center": ""})
        m = re.search(r'GCSF', s)
        if m:
            return pd.Series({"phenotype": "G-CSF-treated healthy donors",
                            "collection_center": "Scripps"})
        m = re.search(r'SDBB', s)
        if m:
            return pd.Series({"phenotype": "Healthy",
                            "collection_center": "San Diego Blood Bank"})
        m = re.search(r'Diverticulitis', s, re.I)
        if m:
            return pd.Series({"phenotype": "Diverticulitis",
                            "collection_center": np.nan})
        return pd.Series(dtype='object')

    df = df.join(df['sample_name'].apply(parse_phenotype))

    # Plasma tubes
    df["plasma_tubes"] = df["biomaterial"].apply(
        lambda x: "EDTA" if x == "plasma" else "BD Vacutainer clotting tubes" if x == "serum" else "")

    # Assign dataset_batch including plasma subtype logic
    def assign_batch(row):
        if row["biomaterial"] == "plasma":
            if row["phenotype"] in ["Multiple myeloma", "Acute Myeloid Leukemia"]:
                return "ibarra_plasma_cancer"
            else:
                return "ibarra_plasma_non_cancer"
        elif row["biomaterial"] == "serum":
            return "ibarra_serum"
        elif row["biomaterial"] == "buffy coat":
            return "ibarra_buffy_coat"
        else:
            return ""

    df["dataset_batch"] = df.apply(assign_batch, axis=1)

    # Assign centrifugation steps based on dataset_batch
    def assign_centrifugation_steps(batch):
        if batch == "ibarra_buffy_coat":
            return pd.Series({"centrifugation_step_1": "1900g", "centrifugation_step_2": "NA"})
        elif batch == "ibarra_serum":
            return pd.Series({"centrifugation_step_1": "1900g", "centrifugation_step_2": "16000g"})
        elif batch == "ibarra_plasma_non_cancer":
            return pd.Series({"centrifugation_step_1": "1900g", "centrifugation_step_2": "16000g"})
        elif batch == "ibarra_plasma_cancer":
            return pd.Series({"centrifugation_step_1": "1900g", "centrifugation_step_2": "6000g"})
        else:
            return pd.Series({"centrifugation_step_1": "NA", "centrifugation_step_2": "NA"})

    df = df.join(df["dataset_batch"].apply(assign_centrifugation_steps))

    df = merge_sample_with_dataset_metadata(
        df, dataset_metadata, keep_sample_cols=["biomaterial", "plasma_tubes"])

    df.to_csv("../sra_metadata/ibarra_metadata_preprocessed.csv", index=False)
    return df


def preprocess_toden(dataset_metadata):
    print("### Dataset: toden")
    csv_path = "../sra_metadata/toden_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "toden" 
    df["dataset_short_name"] = "toden"
    df["dataset_batch"] = "toden"
    df["read_length"] = "2x75"
    df["centrifugation_step_1"] = "12000g"
    df["centrifugation_step_2"] = "NA" 


    df_first = df.drop_duplicates(subset="isolate", keep="first").copy()

    df_grouped = df.groupby("isolate")["run"].apply(lambda x: "|".join(sorted(x))).reset_index()
    df_grouped.columns = ["isolate", "merged_runs"]

    df_merged = df_first.merge(df_grouped, on="isolate", how="left")

    df_merged["run"] = "SRRISOLATE_" + df_merged["isolate"].astype(str)

    cols = df_merged.columns.tolist()
    cols.insert(0, cols.pop(cols.index("run")))
    df_merged = df_merged[cols]

    df_merged['phenotype'] = df_merged['ad_status'].replace({
        'AD':"Alzheimers disease",
        'NCI':'Healthy',
        'None':np.nan
    })

    df_merged = merge_sample_with_dataset_metadata(df_merged, dataset_metadata)

    df_merged.to_csv("../sra_metadata/toden_metadata_preprocessed.csv", index=False)
    return df_merged

def preprocess_chalasani(dataset_metadata):
    print("### Dataset: chalasani")
    csv_path = "../sra_metadata/chalasani_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "chalasani" 
    df["dataset_short_name"] = "chalasani"
    df["dataset_batch"] = "chalasani"
    df["read_length"] = "2x75"
    df["centrifugation_step_1"] = "1900g"
    df["centrifugation_step_2"] = "NA" 


    # Create a new column with the cleaned isolate ID
    df["isolate_base"] = df["isolate"].str.replace(r"-[A-Z]\d*$", "", regex=True)

    # Aggregate Bases by isolate
    df_sum = df.groupby("isolate_base", as_index=False)["bases"].sum()

    df_concat_run = df.groupby("isolate_base")["run"].apply(lambda x: "|".join(sorted(x))).reset_index()
    df_concat_run.columns = ["isolate_base", "merged_runs"]

    df_first = df.drop_duplicates(subset="isolate_base", keep="first").copy()
    df_merged = (df_first.drop(columns=["bases", "run"])
                .merge(df_sum, on="isolate_base", how="left")
                .merge(df_concat_run, on="isolate_base", how="left"))

    # Set Run = isolate_base, and prepend an "X" to match the output files of fl-rnaseq.
    df_merged["run"] = "X" + df_merged["isolate_base"].astype("str")

    cols = df_merged.columns.tolist()
    cols.insert(0, cols.pop(cols.index("run")))
    df_merged = df_merged[cols]
    df_merged = df_merged.drop(columns="isolate_base")

    df_merged = merge_sample_with_dataset_metadata(df_merged, dataset_metadata)

    # Hot-fix: Temporary exclude missing (n=16) Chalasani samples
    chalasani_ids_to_exclude = [
        'X2835',
        'X3391',
        'X3392',
        'X3752',
        'X3770',
        'X3817',
        'X3823',
        'X3827',
        'X3833',
        'X4253',
        'X4269',
        'X4275',
        'X4363',
        'X9709',
        'X9737',
        'X9760'
    ]
    df_merged = df_merged[~df_merged['run'].isin(chalasani_ids_to_exclude)]

    df_merged.to_csv("../sra_metadata/chalasani_metadata_preprocessed.csv", index=False)
    return df_merged


def preprocess_block(dataset_metadata):
    print("### Dataset: block")
    csv_path = "../sra_metadata/block_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "block" 
    df["dataset_short_name"] = "block"
    df["dataset_batch"] = np.where(
        abs(df["avgspotlen"] - 150) < abs(df["avgspotlen"] - 300),
        "block_150bp", "block_300bp"
    )
    df["read_length"] = np.where(df["dataset_batch"] == "block_150bp", "2x75", "2x150")
    df["centrifugation_step_1"] = "2000g"
    df["centrifugation_step_2"] = "NA" 

    # Exclude non-plasma samples: Tissue, and Plasma-derived vesicles.
    df = df.rename(columns={'tissue':'biomaterial'})
    n1 = len(df)
    df = df[df['biomaterial'] == 'Plasma']
    n2 = len(df)
    print(f"Exclude non-plasma samples: Tissue, and Plasma-derived vesicles. N = {n1 - n2}")

    # Add information from supplementary table 7.
    # Characteristics of HCC and CCA patients whose specimens were used. CH, NJ: Capital Health Cancer Center, NJ. UPEN, PA: Veterans hospital University of Pennysylvania, PA. Biochemed.
    # Note that there is no available information for the LC patients (9 samples).
    supp_table = pd.read_excel("../sra_metadata/block_supp_table_7.xlsx", skiprows=3, index_col=1)
    supp_table = supp_table.iloc[1:, :].drop('Unnamed: 0', axis=1)
    supp_table = (supp_table[['Source/Place', 'Bleed date']].dropna(how='all')
                  .reset_index()
                  .rename(columns={'Patient ID':'patient_id', 
                                   'Source/Place':'collection_center'}))
    # Merge with the sample-level metadata dataframe
    df['patient_id'] = df['sample_id'].str.extract(r'(.+\d)[^\d]*$')
    df = df.merge(supp_table, on='patient_id', how='left')

    df = merge_sample_with_dataset_metadata(
        df, dataset_metadata, keep_sample_cols=["biomaterial"])

    df.to_csv("../sra_metadata/block_metadata_preprocessed.csv", index=False)

    return df

def preprocess_rozowsky(dataset_metadata):
    print("### Dataset: block")
    # Warning: this metadata does not come from SRA 
    csv_path = "../sra_metadata/rozowsky_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "rozowsky" 
    df["dataset_short_name"] = "rozowsky"
    df["dataset_batch"] = "rozowsky"
    df['phenotype'] = 'Healthy'
    df["read_length"] = "2x100"
    df["centrifugation_step_1"] = "NA"
    df["centrifugation_step_2"] = "NA" 

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/rozowsky_metadata_preprocessed.csv", index=False)

    return df

def preprocess_tao(dataset_metadata):
    print("### Dataset: tao")
    csv_path = "../sra_metadata/tao_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "tao" 
    df["dataset_short_name"] = "tao"
    df["dataset_batch"] = "tao"
    df["read_length"] = "2x150"
    df["centrifugation_step_1"] = "Unspecified"
    df["centrifugation_step_2"] = "Unspecified"    

    # Select only the cfRNA-seq samples, filter out tissue and PBMC, and other assays like MeDIP-Seq, miRNA-Seq
    n1 = len(df)
    df = df[(df['source_name'] == 'plasma') &
            (df['assay_type'] == 'RNA-Seq')]
    n2 = len(df)
    print(f"Exclude tissue and PBMC, and other assays like MeDIP-Seq, miRNA-Seq. N = {n1 - n2}")

    # Parse the GEO series matrix file, which contains the mapping between
    # the GEO/GSM ids and the sample names in the study
    GSM_ids = []
    sample_titles = []
    for matrix_file in ["../sra_metadata/tao_GSE186607_series_matrix.txt"]:
        with open(matrix_file) as f:
            GEO_matrix = f.readlines()

        line = [line for line in GEO_matrix if re.search(r'!Sample_geo_accession', line)][0]
        s = re.search(r'^!Sample_geo_accession\t"(.+)"\s*', line).group(1)
        GSM_ids_1 = s.split()
        GSM_ids_1 = [e.strip('" ') for e in GSM_ids_1]
        GSM_ids += GSM_ids_1
        
        line = [line for line in GEO_matrix if re.search(r'!Sample_title', line)][0]
        s = re.search(r'^!Sample_title\t"(.+)"\s*', line).group(1).strip()
        sample_titles_1 = s.split('"\t"')
        sample_titles += sample_titles_1

    if len(sample_titles) != len(GSM_ids):
        raise RuntimeError("Parsing the series matrix file, found different number of sample ids and sample titles.")
    matrix_metadata = pd.DataFrame(np.transpose([GSM_ids, sample_titles]),
                                columns=['GSM_id', 'sample_title'])

    df = (df.merge(matrix_metadata, left_on='sample_name',
                right_on='GSM_id', how='left'))
    df['phenotype'] = (df['sample_title'].str.extract(r'^(.+)-PKU.*$')
                    .replace({
                        'CRC':'Colorectal cancer',
                        'NC':'Healthy',
                        'STAD':'Stomach cancer'
                    }))

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/tao_metadata_preprocessed.csv", index=False)

    return df

def preprocess_wei(dataset_metadata):
    print("### Dataset: wei")
    csv_path = "../sra_metadata/wei_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "taowei" 
    df["dataset_short_name"] = "wei"
    df["dataset_batch"] = "wei"
    df["read_length"] = "2x150"    
    df["centrifugation_step_1"] = "3000g"
    df["centrifugation_step_2"] = "12000g" 

    # Select only the plasma samples, remove the tissue "tissue" ones
    n1 = len(df)
    df = df[(df['tissue'] == 'plasma')]
    n2 = len(df)
    print(f"Exclude tissue samples. N = {n1 - n2}")

    df['phenotype'] = "Pancreatic cancer"

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/wei_metadata_preprocessed.csv", index=False)

    return df

def preprocess_moufarrej(dataset_metadata):
    print("### Dataset: moufarrej")
    # Note: we need to find a way to separate this dataset into the "sites" and into the cohorts. also decide if we want to merge by replicate
    csv_path = "../sra_metadata/moufarrej_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)
    # Improve compatibility with snakeDA (reserved var: ['sequencing_batch'])
    df = df.rename(columns={'sequencing_batch':'sequencing_batch_other'})

    df["sequencing_batch"] = "moufarrej" 
    df["dataset_short_name"] = "moufarrej"
    #df["dataset_batch"] = "moufarrej" # use 'cohort' column from GEO
    df["read_length"] = "2x75"    

    # Dataset_metadata table says EDTA/Streck ?
    #df["plasma_tubes"] = "EDTA"

    # Cohorts is composed of preeclampsia and healthy control (normotensive)
    #df.loc[df['disease'].isnull(), 'phenotype'] = 'Healthy pregnant women'

    # Parse GEO metadata
    geo_series_meta = "../sra_metadata/moufarrej_geo_series_metadata.csv"
    if not os.path.exists(geo_series_meta):

        # Load series matrix file
        gse = GEOparse.get_GEO(geo="GSE192902", destdir="../", silent=True)

        # loop over sample data
        meta_dict = {}
        for gsm_name, gsm in gse.gsms.items():
            library_name = gsm.metadata['title'][0]
            meta_dict[library_name] = {ss.split(": ")[0]:ss.split(": ")[1] for ss in gsm.metadata['characteristics_ch1']}
            
        # create DataFrame
        geo_df = (
            pd.DataFrame
            .from_dict(meta_dict, orient='index')
            .reset_index()
            .rename(columns={'index':'library_name'})
        )
        geo_df.columns = simplify_column_names(geo_df.columns)
        
        geo_cols = ["library_name", "disease", "sampling_time_group", "cohort"]
        # Store processed GEO metadata to file
        geo_df = geo_df[geo_cols]
        geo_df.to_csv(geo_series_meta, index=False)
        
    # Load processed GEO metadata
    geo_df = pd.read_csv(geo_series_meta).rename(columns={'disease':'phenotype'})
    # add GEO metadata
    df = df.merge(geo_df, on='library_name', how='left')

    # re-map 'cohort' to 'collection_center'
    moufarrej_cohort_map = {
        'Discovery':    'Lucile Packard Children Hospital',
        'Validation 1': 'Lucile Packard Children Hospital',
        'Validation 2': 'Global Alliance to Prevent Prematurity and Stillbirth (GAPPS)',
    }

    df["collection_center"] = df["cohort"].map(moufarrej_cohort_map)

    def assign_site(center):
        if center == "Lucile Packard Children Hospital":
            return "moufarrej_site_1"
        elif center == "Global Alliance to Prevent Prematurity and Stillbirth (GAPPS)":
            return "moufarrej_site_2"
        else:
            return "moufarrej_unknown"

    df["dataset_batch"] = df["collection_center"].apply(assign_site)
    df['collection_center'] = df['cohort'].apply(lambda x: moufarrej_cohort_map[x])


    def assign_centrifugation_steps(batch):
        if batch == "moufarrej_site_1":
            return pd.Series({"centrifugation_step_1": "1600g", "centrifugation_step_2": "13000g"})
        elif batch == "moufarrej_site_2":
            return pd.Series({"centrifugation_step_1": "2500g", "centrifugation_step_2": "NA"})
        else:
            return pd.Series({"centrifugation_step_1": "NA", "centrifugation_step_2": "NA"})

    df = df.join(df["dataset_batch"].apply(assign_centrifugation_steps))

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)
    df.to_csv("../sra_metadata/moufarrej_metadata_preprocessed.csv", index=False)

    return df

def preprocess_wang(dataset_metadata):
    print("### Dataset: wang")
    csv_path = "../sra_metadata/wang_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "wang_read_2" 
    df["dataset_short_name"] = "wang"
    df["dataset_batch"] = "wang"
    df["read_length"] = "2x150"    
    df["centrifugation_step_1"] = "1600g"
    df["centrifugation_step_2"] = "16000g" 


    # There are 3 samples that were processed with a different library prep kit.
    # Discard them.
    n1 = len(df)
    index = df[df['sample_name'].str.startswith("NEB")].index
    df.loc[index, "library_prep_kit"] = "NEBNext Small RNA Library Prep Set"
    df.loc[index, "library_prep_kit_short"] = "NEBNext"
    df = df.loc[df.index.difference(index)]
    n2 = len(df)
    print(f"Exclude other library prep. N = {n1 - n2}")

    def parse_plasma_volume_wang(s):
        m = re.search(r'^Input\(([\d.]+?)\).*$', s, re.I)
        if m:
            return float(m.group(1))
        else:
            return np.nan

    df['plasma_volume'] = df['sample_name'].map(parse_plasma_volume_wang).dropna()
    
    # During technology optimization, different modifications of the library prep
    # protocol were tested. We discard the samples for which irrelevant modifications
    # were applied to the protocol, such as removing the ExoI enzyme or changing
    # the concentration of the RT primers.
    n1 = len(df)
    df = (df.set_index('sample_name').drop([
        'USER(-)-1',
        'ExoI(-)-3',
        'ExoI(-)-2',
        'RT(80)-3',
        'RT(80)-2',
        'RT(80)-1',
        'RT(40)-3',
        'RT(40)-2',
        'RT(40)-1',
        'RT(20)-3',
        'RT(20)-2',
        'RT(20)-1',
        'ExoI(-)-1',
        'RT(10)-3',
        'RT(10)-2',
        'RT(10)-1',
        'RT(2.5)-3',
        'RT(2.5)-2',
        'RT(2.5)-1',
        'RT(1.25)-3',
        'RT(1.25)-2',
        'RT(1.25)-1',
        'USER(-)-3',
        'USER(-)-2',
    ]).reset_index())
    n2 = len(df)
    print(f"Exclude other library prep protocols. N = {n1 - n2}")

    # All samples from SRA are from the technology optimizations experiment,
    # in which all samples are healthy.
    df['phenotype'] = "Healthy"

    df = merge_sample_with_dataset_metadata(
        df, dataset_metadata, keep_sample_cols=["library_prep_kit",
                                                "library_prep_kit_short",
                                                "plasma_volume"])

    df.to_csv("../sra_metadata/wang_metadata_preprocessed.csv", index=False)

    return df

def preprocess_giraldez(dataset_metadata):
    print("### Dataset: giraldez")
    csv_path = "../sra_metadata/giraldez_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "giraldez" 
    df["dataset_short_name"] = "giraldez"
    
    # Filter out the 2 synthetic samples
    n1 = len(df)
    df = df[~df['source_name'].str.contains('Synthetic sRNA equimolar pool')]
    n2 = len(df)
    print(f"Exclude synthetic samples. N = {n1 - n2}")
    
    # Filter out protocol optimization samples
    n1 = len(df)
    df = df[df["sra_study"] == "SRP183468"]
    n2 = len(df)
    print(f"Exclude protocol optimization samples. N = {n1 - n2}")

    # Standard library prep
    index = df[df['treatment'].isin(['none', 'Untreated'])].index
    df.loc[index, "library_prep_kit"] = "Illumina TruSeq small RNA"
    df.loc[index, "library_prep_kit_short"] = "Illumina TruSeq small RNA"
    df.loc[index, "assay_name"] = "RNA-seq"
    df.loc[index, "dataset_batch"] = "giraldez_standard"
    df.loc[index, "centrifugation_step_1"] = "3400g"
    df.loc[index, "centrifugation_step_2"] = "1940g" 

    # phospho-RNA-seq library prep
    index = df[df['treatment'].isin(['T4PNK', 'PNK'])].index
    df.loc[index, "library_prep_kit"] = "polynucleotide kinase (PNK) treated, Illumina TruSeq small RNA"
    df.loc[index, "library_prep_kit_short"] = "PNK-treated Illumina TruSeq small RNA"
    df.loc[index, "assay_name"] = "phospho-RNA-seq"
    df.loc[index, "dataset_batch"] = "giraldez_phospho-rna-seq"
    df.loc[index, "centrifugation_step_1"] = "3400g"
    df.loc[index, "centrifugation_step_2"] = "1940g" 

    df["read_length"] = np.where(df["dataset_batch"] == "giraldez_standard", "1x50", "1x75")

    # All remaining samples are from 5 healthy donors
    df.loc[df['sra_study'] == 'SRP183467', 'phenotype'] = "Healthy"

    df = merge_sample_with_dataset_metadata(
        df, dataset_metadata, keep_sample_cols=["library_prep_kit",
                                                "library_prep_kit_short",
                                                "assay_name"])

    df.to_csv("../sra_metadata/giraldez_metadata_preprocessed.csv", index=False)

    return df

def preprocess_sun(dataset_metadata):
    print("### Dataset: sun")
    csv_path = "../sra_metadata/sun_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "sun" 
    df["dataset_short_name"] = "sun"
    df["biomaterial"] = df["tissue"].apply(lambda x: "blood plasma" if x == "plasma" else "blood serum" if x == "serum" else "")
    df["dataset_batch"] = np.where(df["biomaterial"] == "blood plasma", "sun_1", "sun_2")
    df["read_length"] = "2x150"
    df["centrifugation_step_1"] = "1500g"
    df["centrifugation_step_2"] = "3000g" 

    # Exclude cfDNA samples
    n1 = len(df)
    df = df[df['assay_type'] == 'RNA-Seq']
    n2 = len(df)
    print(f"Exclude cfDNA samples. N = {n1 - n2}")

    # Assign sequencing batches based on the sequencing platform instrument
    # These also correspond to two separate datasets on GEO
    # => 'sequencing_batch_other': Improve compatibility with snakeDA (reserved var: ['sequencing_batch'])
    df.loc[df['instrument'] == 'HiSeq X Ten', 'sequencing_batch_other'] = 's1'
    df.loc[df['instrument'] == 'Illumina NovaSeq 6000', 'sequencing_batch_other'] = 's2'

    # Parse the GEO series matrix file, which contains the mapping between
    # the GEO/GSM ids and the sample names in the study
    GSM_ids = []
    sample_titles = []
    for matrix_file in ["../sra_metadata/sun_GSE210775-GPL20795_series_matrix.txt",
                        "../sra_metadata/sun_GSE210775-GPL24676_series_matrix.txt"]:
        with open(matrix_file) as f:
            GEO_matrix = f.readlines()

        line = [line for line in GEO_matrix if re.search(r'!Sample_geo_accession', line)][0]
        s = re.search(r'^!Sample_geo_accession\t"(.+)"\s*', line).group(1)
        GSM_ids_1 = s.split()
        GSM_ids_1 = [e.strip('" ') for e in GSM_ids_1]
        GSM_ids += GSM_ids_1
        
        line = [line for line in GEO_matrix if re.search(r'!Sample_title', line)][0]
        s = re.search(r'^!Sample_title\t"(.+)"\s*', line).group(1).strip()
        sample_titles_1 = s.split('"\t"')
        sample_titles += sample_titles_1

    if len(sample_titles) != len(GSM_ids):
        raise RuntimeError("Parsing the series matrix file, found different number of sample ids and sample titles.")
    matrix_metadata = pd.DataFrame(np.transpose([GSM_ids, sample_titles]),
                                columns=['GSM_id', 'sample_title'])
    # RNA extraction kit (inferred from "Internal ID" and manuscript)
    matrix_metadata['rna_extraction_kit'] = np.nan
    matrix_metadata['rna_extraction_kit_short_name'] = np.nan
    index = (matrix_metadata[matrix_metadata['sample_title']
             .str.contains('with_MOF|Training_set|Validation_set')].index)
    matrix_metadata.loc[index, 'rna_extraction_kit'] = "MOF method"
    matrix_metadata.loc[index, 'rna_extraction_kit_short_name'] = "MOF"

    df = (df.merge(matrix_metadata, left_on='sample_name',
                   right_on='GSM_id', how='left'))

    # Exclude the 15 protocol optimization samples
    n1 = len(df)
    GSM_ids_to_exclude = """GSM6437012
    GSM6437013
    GSM6437014
    GSM6437015
    GSM6437016
    GSM6437017
    GSM6437018
    GSM6437019
    GSM6437020
    GSM6437021
    GSM6437022
    GSM6437023
    GSM6437024
    GSM6437025
    GSM6437026""".split()
    df = df[~df['sample_name'].isin(GSM_ids_to_exclude)]
    n2 = len(df)
    print(f"Exclude the 15 protocol optimization samples. N = {n1 - n2}")

    df = merge_sample_with_dataset_metadata(
        df, dataset_metadata, keep_sample_cols=["biomaterial", "rna_extraction_kit", "rna_extraction_kit_short_name"])

    df.to_csv("../sra_metadata/sun_metadata_preprocessed.csv", index=False)

    return df

def preprocess_decruyenaere(dataset_metadata):
    print("### Dataset: decruyenaere")
    # Most relevant variables:
    col_names = [
        'sample_alias',
        # 'sample_accession_id',
        # 'biosample_id',
        # 'run_accession_id',
        # 'experiment_accession_id',
        # 'study_accession_id',
        # 'instrument_platform',
        'instrument_model',
        # 'library_layout',
        # 'library_name',
        # 'library_strategy',
        # 'library_source',
        # 'library_selection',
        'sample_title',
        # 'run_file_type',
        # 'design_description',
        'description',
        'biological_sex',
        'subject_id',
        'phenotype',
        'case_or_control',
        # 'cell_line',
        # 'ENA-CHECKLIST',
        'organism_part',
        # 'region',
        #'file_name_1' # needed to extract the 'sample_name'
    ]

    csv_path = "../sra_metadata/decruyenaere_metadata.csv"
    csv_path = "../sra_metadata/decruyenaere_metadata_ext.tsv"
    df = pd.read_csv(csv_path, sep='\t')[col_names]
    df.columns = simplify_column_names(df.columns)

    # Rename columns from the EGA-archive to the corresponding SRA metadata columns
    df = df.rename(columns={
        "sample_alias":"run",
        "instrument_model":"instrument",
        "biosample_id":"biosample",
        })

    #df["run"] = df["file_name_1"].apply(lambda x: '_'.join(x.split('_')[:2]))
    # remove 'file_name_1' col
    #df = df.drop('file_name_1', axis=1)

    df["sequencing_batch"]   = "decruyenaere" 
    df["dataset_short_name"] = "decruyenaere"
    df["dataset_batch"]      = "decruyenaere"
    df["disease"]            = df["phenotype"]
    df["read_length"]        = "2x100"
    df["centrifugation_step_1"] = "1900g"
    df["centrifugation_step_2"] = "NA" 

    # Exclude FFPE samples
    n1 = len(df)
    df = df[df['organism_part'] != 'FFPE homo sapiens']
    n2 = len(df)
    print(f"Exclude FFPE samples. N = {n1 - n2}")

    # All remaining samples are blood plasma, as described in the column "organism_part"
    df['biomaterial'] = df['organism_part'].replace({'blood plasma homo sapiens':'plasma'})
    df = df.drop(columns=['organism_part'])

    df = merge_sample_with_dataset_metadata(df, dataset_metadata,
                                            keep_sample_cols=["biomaterial"])

    df.to_csv("../sra_metadata/decruyenaere_metadata_preprocessed.csv", index=False)

    return df

def preprocess_reggiardo(dataset_metadata):
    print("### Dataset: reggiardo")
    csv_path = "../sra_metadata/reggiardo_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["sequencing_batch"] = "reggiardo" 
    df["dataset_short_name"] = "reggiardo"
    df["read_length"] = "2x150"
    df["centrifugation_step_1"] = "Unspecified"
    df["centrifugation_step_2"] = "Unspecified" 
    
    # Assign collection_center based on whether 'stage' is missing
    df["collection_center"] = df["subject_status"].str.strip().map({
        "normal healthy donor": "Discovery Life Sciences",
        "pancreatic cancer patient": "BioIVT"
    })

    # Assign dataset_batch accordingly
    df["dataset_batch"] = df["collection_center"].map({
        "Discovery Life Sciences": "reggiardo_dls",
        "BioIVT": "reggiardo_bioivt"
    })

    # Select only Illumina samples and exclude ONT samples
    n1 = len(df)
    df = df[(df['platform'] == 'ILLUMINA')]
    n2 = len(df)
    print(f"Exclude ONT samples. N = {n1 - n2}")

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/reggiardo_metadata_preprocessed.csv", index=False)

    return df


def preprocess_flomics_1(dataset_metadata):
    print("### Dataset: flomics_1")
    csv_path = "../sra_metadata/flomics_1_metadata.tsv"
    df = pd.read_csv(csv_path, sep='\t')
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "flomics_1"
    df["dataset_batch"] = "flomics_1"
    df["read_length"] = "2x150"
    df["centrifugation_step_1"] = "placeholder"
    df["centrifugation_step_2"] = "placeholder" 
    
    df["run"] = df["sample_name"]
    df["run"] = df["sample_analysis_run_id"].apply(lambda x: x.split('_')[0])

    # exclude samples
    samples_to_remove = ["Flomics_1_1", "Flomics_1_2"]
    df = df[~df["sample_name"].isin(samples_to_remove)]
    
    reserved_cols = set(df.columns).intersection(reserved_vars_samplesheet)
    if reserved_cols:
        df = df.drop(reserved_cols, axis=1)

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/flomics_1_metadata_preprocessed.csv", index=False)

    return df


def preprocess_flomics_2(dataset_metadata):
    print("### Dataset: flomics_2")
    csv_path = "../sra_metadata/flomics_2_metadata.tsv"
    df = pd.read_csv(csv_path, sep='\t')
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "flomics_2"
    df["dataset_batch"] = "flomics_2"
    df["read_length"] = "2x150"
    df["centrifugation_step_1"] = "1500g"
    df["centrifugation_step_2"] = "2500g" 

    #df["run"] = df["sample_name"]
    df["run"] = df["sample_analysis_run_id"].apply(lambda x: x.split('_')[0])
    
    reserved_cols = set(df.columns).intersection(reserved_vars_samplesheet)
    if reserved_cols:
        df = df.drop(reserved_cols, axis=1)

    df = df.rename(columns={"status_subtype":"phenotype"})

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/flomics_2_metadata_preprocessed.csv", index=False)

    return df


def summarize_metadata_batch_level(sample_metadata):
    sample_metadata2 = sample_metadata.dropna(axis=1, how='all').copy()
    # Count how many unique values each column contains
    sample_metadata2.loc[:, 'n_samples'] = sample_metadata2.groupby('dataset_batch')['run'].transform('count')
    n_values_per_group = sample_metadata2.groupby('dataset_batch').nunique()
    n_values_per_group2 = n_values_per_group.reset_index()
    n_values_per_group2['dataset_batch'] = 1
    cols_with_unique_value = n_values_per_group2.apply(lambda col: all(col <= 1), axis=0)
    cols_with_nonunique_value = n_values_per_group2.apply(lambda col: any(col > 1), axis=0)
    cols_with_nonunique_value = cols_with_nonunique_value[cols_with_nonunique_value].index.to_list()

    # Columns with multiple values in the same batch
    # Some of these columns might be used to define new batches within datasets
    # for col in cols_with_nonunique_value:
    #     print("\nColumn:", col)
    #     batch_list = n_values_per_group[n_values_per_group[col] > 1].index
    #     col_values_in_batch = sample_metadata2[sample_metadata2['dataset_batch'].isin(batch_list)][['dataset_batch', col]].drop_duplicates()
    #     print(col_values_in_batch)
    # interesting_cols = [
    #     'instrument',
    #     'collection_center',
    #     'plasma_tubes',
    #     'plasma_volume',
    #     'library_preparation_batch',
    #     'rna_extraction_batch'
    # ]

    # Keep only columns which contain a unique value within each batch, drop all other columns
    sample_metadata3 = sample_metadata2.loc[:, cols_with_unique_value].copy()
    sample_metadata3 = sample_metadata3.dropna(axis=1, how='all')
    batch_metadata = sample_metadata3.drop_duplicates()
    if batch_metadata['dataset_batch'].nunique() != len(batch_metadata):
        raise ValueError("The dataset_batch columns are duplicated values.")
    # Reorder columns
    first_cols = ['dataset_short_name', 'dataset_batch']
    batch_metadata = pd.concat([batch_metadata[first_cols],
                                batch_metadata[[c for c in batch_metadata.columns if c not in first_cols]]],
                                axis=1)
    return batch_metadata


def main():
    # Read plasma cfRNA-seq workflow comparison
    csv_path = "../sra_metadata/dataset_metadata.tsv"
    dataset_metadata = pd.read_csv(csv_path, sep = "\t")
    # Simplify column names
    dataset_metadata.columns = simplify_column_names(dataset_metadata.columns)

    chen = preprocess_chen(dataset_metadata)
    zhu = preprocess_zhu(dataset_metadata)
    roskams = preprocess_roskams(dataset_metadata)
    ngo = preprocess_ngo(dataset_metadata)
    ibarra = preprocess_ibarra(dataset_metadata)
    toden = preprocess_toden(dataset_metadata)
    chalasani = preprocess_chalasani(dataset_metadata)
    block = preprocess_block(dataset_metadata)
    rozowsky = preprocess_rozowsky(dataset_metadata)
    tao = preprocess_tao(dataset_metadata)
    wei = preprocess_wei(dataset_metadata)
    moufarrej = preprocess_moufarrej(dataset_metadata)
    wang = preprocess_wang(dataset_metadata)
    giraldez = preprocess_giraldez(dataset_metadata)
    sun = preprocess_sun(dataset_metadata)
    decruyenaere = preprocess_decruyenaere(dataset_metadata)
    reggiardo = preprocess_reggiardo(dataset_metadata)
    flomics_1 = preprocess_flomics_1(dataset_metadata)
    flomics_2 = preprocess_flomics_2(dataset_metadata)

    dfs = [chen, zhu, roskams, ngo, ibarra, toden, chalasani, block, rozowsky, tao, wei, moufarrej, wang, giraldez, sun, decruyenaere, reggiardo, flomics_1, flomics_2]

    sample_metadata = pd.concat(dfs, axis=0, join="outer", ignore_index=True)

    sample_metadata = rename_columns_and_values(sample_metadata)

    sample_metadata.to_csv("../tables/cfRNA-meta_per_sample_metadata.tsv", sep="\t", index=False)

    # Summarize metadata at the batch level
    batch_metadata = summarize_metadata_batch_level(sample_metadata)
    batch_metadata.to_csv("../tables/cfRNA-meta_per_batch_metadata.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
