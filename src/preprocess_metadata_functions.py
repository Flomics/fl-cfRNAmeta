import pandas as pd
import numpy as np
import re
import pdb
from sra_columns_mapping import rename_columns_and_values
pd.set_option('display.max_columns', None)



def simplify_column_names(cols):
    if type(cols) is not pd.Series:
        cols = pd.Series(cols)
    cols = (cols
            .str.lower()
            .str.replace(r'\s', r'_', regex=True)
            .str.replace(r'[()]', r'', regex=True)
           ).to_list()
    return cols


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
    csv_path = "../sra_metadata/chen_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "chen"
    df["dataset_batch"] = "chen"
    # Exclude the two E. coli samples and the brain tissue sample
    df = df[~df['sample_name'].isin([
        'ET_L2'
        'EH_L2'
        'NC_L2'
    ])]

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/chen_metadata_preprocessed.csv", index=False)
    return df

def preprocess_zhu(dataset_metadata):
    csv_path = "../sra_metadata/zhu_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "zhu"
    df["dataset_batch"] = "zhu"

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/zhu_metadata_preprocessed.csv", index=False)
    return df

def preprocess_roskams(dataset_metadata):
    csv_path = "../sra_metadata/roskams_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)    

    df["dataset_short_name"] = "roskams"
    df["dataset_batch"] = np.where(df["cohort"] == "pilot", "roskams_1", "roskams_2")

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
    csv_path = "../sra_metadata/ngo_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "ngo"
    df["dataset_batch"] = "ngo"

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/ngo_metadata_preprocessed.csv", index=False)
    return df

def preprocess_ibarra(dataset_metadata):
    csv_path = "../sra_metadata/ibarra_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "ibarra"
    df["dataset_batch"] = "ibarra"
    # We overwrite the values that we merge above from the dataset_metadata table
    df["biomaterial"] = df["tissue"].str.lower()
    df["plasma_tubes"] = df["biomaterial"].apply(lambda x: "EDTA" if x == "plasma" else "BD Vacutainer clotting tubes" if x == "serum" else "")
    # TODO Parse the phenotype from the Sample Name

    df = merge_sample_with_dataset_metadata(
        df, dataset_metadata, keep_sample_cols=["biomaterial", "plasma_tubes"])

    df.to_csv("../sra_metadata/ibarra_metadata_preprocessed.csv", index=False)
    return df

def preprocess_toden(dataset_metadata):
    csv_path = "../sra_metadata/toden_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "toden"
    df["dataset_batch"] = "toden"

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/toden_metadata_preprocessed.csv", index=False)
    return df

def preprocess_chalasani(dataset_metadata):
    csv_path = "../sra_metadata/chalasani_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "chalasani"
    df["dataset_batch"] = "chalasani"
    # TODO parse the donor id from the sample name
    # TODO define sample id for the merged fastq file processed by fl-rnaseq.

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/chalasani_metadata_preprocessed.csv", index=False)
    return df

def preprocess_block(dataset_metadata):
    csv_path = "../sra_metadata/block_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "block"
    df["dataset_batch"] = np.where(
        abs(df["avgspotlen"] - 150) < abs(df["avgspotlen"] - 300),
        "block_1", "block_2"
    )
    # Exclude non-plasma samples: Tissue, and Plasma-derived vesicles.
    df = df.rename(columns={'tissue':'biomaterial'})
    df = df[df['biomaterial'] == 'Plasma']

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
    # Warning: this metadata does not come from SRA 
    csv_path = "../sra_metadata/rozowsky_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "rozowsky"
    df["dataset_batch"] = "rozowsky"
    
    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/rozowsky_metadata_preprocessed.csv", index=False)

    return df

def preprocess_tao(dataset_metadata):
    csv_path = "../sra_metadata/tao_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "tao"
    df["dataset_batch"] = "tao"
    # Select only the cfRNA-seq samples, filter out tissue and PBMC, and other assays like MeDIP-Seq, miRNA-Seq
    df = df[(df['source_name'] == 'plasma') &
            (df['assay_type'] == 'RNA-Seq')]

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/tao_metadata_preprocessed.csv", index=False)

    return df

def preprocess_wei(dataset_metadata):
    csv_path = "../sra_metadata/wei_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "wei"
    df["dataset_batch"] = "wei"
    # Select only the plasma samples, remove the tissue "tissue" ones
    df = df[(df['tissue'] == 'plasma')]

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/wei_metadata_preprocessed.csv", index=False)

    return df

def preprocess_moufarrej(dataset_metadata):
    # Note: we need to find a way to separate this dataset into the "sites" and into the cohorts. also decide if we want to merge by replicate
    csv_path = "../sra_metadata/moufarrej_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "moufarrej"
    df["dataset_batch"] = "moufarrej"

    # Dataset_metadata table says EDTA/Streck ?
    #df["plasma_tubes"] = "EDTA"

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/moufarrej_metadata_preprocessed.csv", index=False)

    return df

def preprocess_wang(dataset_metadata):
    csv_path = "../sra_metadata/wang_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "wang"
    df["dataset_batch"] = "wang"

    # There are 3 samples that were processed with a different library prep kit.
    # Discard them.
    index = df[df['sample_name'].str.startswith("NEB")].index
    df.loc[index, "library_prep_kit"] = "NEBNext Small RNA Library Prep Set"
    df.loc[index, "library_prep_kit_short"] = "NEBNext"
    df = df.loc[df.index.difference(index)]

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

    df = merge_sample_with_dataset_metadata(
        df, dataset_metadata, keep_sample_cols=["library_prep_kit",
                                                "library_prep_kit_short",
                                                "plasma_volume"])

    df.to_csv("../sra_metadata/wang_metadata_preprocessed.csv", index=False)

    return df

def preprocess_giraldez(dataset_metadata):
    csv_path = "../sra_metadata/giraldez_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "giraldez"
    
    # Filter out the 2 synthetic samples
    df = df[~df['source_name'].str.contains('Synthetic sRNA equimolar pool')]
    
    # Filter out protocol optimization samples
    df = df[df["sra_study"] == "SRP183468"]

    # Standard library prep
    index = df[df['treatment'].isin(['none', 'Untreated'])].index
    df.loc[index, "library_prep_kit"] = "Illumina TruSeq small RNA"
    df.loc[index, "library_prep_kit_short"] = "Illumina TruSeq small RNA"
    df.loc[index, "assay_name"] = "RNA-seq"
    df.loc[index, "dataset_batch"] = "giraldez_1"

    # phospho-RNA-seq library prep
    index = df[df['treatment'].isin(['T4PNK', 'PNK'])].index
    df.loc[index, "library_prep_kit"] = "polynucleotide kinase (PNK) treated, Illumina TruSeq small RNA"
    df.loc[index, "library_prep_kit_short"] = "PNK-treated Illumina TruSeq small RNA"
    df.loc[index, "assay_name"] = "phospho-RNA-seq"
    df.loc[index, "dataset_batch"] = "giraldez_2"

    df = merge_sample_with_dataset_metadata(
        df, dataset_metadata, keep_sample_cols=["library_prep_kit",
                                                "library_prep_kit_short",
                                                "assay_name"])

    df.to_csv("../sra_metadata/giraldez_metadata_preprocessed.csv", index=False)

    return df

def preprocess_sun(dataset_metadata):
    csv_path = "../sra_metadata/sun_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "sun"
    df["biomaterial"] = df["tissue"].apply(lambda x: "blood plasma" if x == "plasma" else "blood serum" if x == "serum" else "")
    df["dataset_batch"] = np.where(df["biomaterial"] == "blood plasma", "sun_1", "sun_2")

    df = merge_sample_with_dataset_metadata(
        df, dataset_metadata, keep_sample_cols=["biomaterial"])

    df.to_csv("../sra_metadata/sun_metadata_preprocessed.csv", index=False)

    return df

def preprocess_decruyenaere(dataset_metadata):
    # TODO wait on Pablo to push his branch
    # Most relevant variables:
    col_names = [
        'sample_alias', 'sample_accession_id', 'biosample_id', 'run_accession_id', 'experiment_accession_id', 'study_accession_id', 'instrument_platform', 'instrument_model', 'library_layout', 'library_name', 'library_strategy', 'library_source', 'library_selection', 'run_file_type', 'design_description', 'description',
        'sample_title', 'phenotype', 'case_or_control', 'biological_sex', 'subject_id', 'cell_line', 'ENA-CHECKLIST',
        'organism_part', 'region'
    ]
    
    csv_path = "../sra_metadata/decruyenaere_metadata.csv"
    # TODO include parsing of the original metadata column "sample_attributes"
    csv_path = "../sra_metadata/decruyenaere_metadata_ext.tsv"
    df = pd.read_csv(csv_path, sep='\t')[col_names]
    df.columns = simplify_column_names(df.columns)
    
    # only dataset from EGA-archive
    df = df.rename(columns={"sample_alias":"Run"})

    df["dataset_short_name"] = "decruyenaere"
    df["dataset_batch"] = "decruyenaere"
    df = df.drop('library_selection', axis=1)
    
    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/decruyenaere_metadata_preprocessed.csv", index=False)

    return df

def preprocess_reggiardo(dataset_metadata):
    csv_path = "../sra_metadata/reggiardo_metadata.csv"
    df = pd.read_csv(csv_path)
    df.columns = simplify_column_names(df.columns)

    df["dataset_short_name"] = "reggiardo"
    df["dataset_batch"] = "reggiardo"

    # Select only Illumina samples and exclude ONT samples
    df = df[(df['platform'] == 'ILLUMINA')]

    df = merge_sample_with_dataset_metadata(df, dataset_metadata)

    df.to_csv("../sra_metadata/reggiardo_metadata_preprocessed.csv", index=False)

    return df


if __name__ == "__main__":
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

    dfs = [chen, zhu, roskams, ngo, ibarra, toden, chalasani, block, rozowsky, tao, wei, moufarrej, wang, giraldez, sun, decruyenaere, reggiardo]

    merged_metadata = pd.concat(dfs, axis=0, join="outer", ignore_index=True)

    merged_metadata = rename_columns_and_values(merged_metadata)

    merged_metadata.to_csv("../tables/cfRNA-meta_per_sample_metadata.tsv", sep="\t", index=False)
