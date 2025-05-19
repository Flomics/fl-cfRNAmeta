import pandas as pd
import numpy as np
import re
import pdb


def preprocess_chen():
    csv_path = "../sra_metadata/chen_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "chen"
    df["dataset_batch"] = "chen"
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "EDTA-coated vacutainer"
    #df["rna_extraction_kit"]="Norgen Plasma/Serum Circulating RNA and Exosomal Purification kit"
    #df["rna_extraction_kit_short"]="Norgen"
    #df["dnase"]="DNase I"
    #df["library_prep_kit"]="SMARTer Stranded Pico v2"
    #df["library_prep_kit_short"]="SMARTer Pico v2"
    df.to_csv("../sra_metadata/chen_metadata_preprocessed.csv", index=False)

    return df

def preprocess_zhu():
    csv_path = "../sra_metadata/zhu_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "zhu"
    df["dataset_batch"] = "zhu"
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "Unspecified"
    #df["rna_extraction_kit"]="Norgen Plasma/Serum Circulating RNA and Exosomal Purification kit"
    #df["rna_extraction_kit_short"]="Norgen"
    #df["dnase"]="DNase I"
    #df["library_prep_kit"]="SMARTer Stranded Pico v2"
    #df["library_prep_kit_short"]="SMARTer Pico v2"
    df.to_csv("../sra_metadata/zhu_metadata_preprocessed.csv", index=False)

    return df

def preprocess_roskams():
    csv_path = "../sra_metadata/roskams_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "roskams"
    df["dataset_batch"] = np.where(df["Cohort"] == "pilot", "roskams_1", "roskams_2")
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "EDTA-anticoagulated vacutainers"
    #df["rna_extraction_kit"]="Norgen Plasma/Serum Circulating RNA and Exosomal Purification kit"
    #df["rna_extraction_kit_short"]="Norgen"
    #df["dnase"]="Baseline-ZERO"
    #df["library_prep_kit"]="SMARTer Stranded Pico v2"
    #df["library_prep_kit_short"]="SMARTer Pico v2"

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
        raise RunTimeError("Parsing the series matrix file, found different number of sample ids and sample names.")
    sample_id_mapping = pd.DataFrame(np.transpose([GSM_ids, sample_names]),
                                     columns=['GSM_id', 'Sample_id'])
    # Merge the metadata table with the supp table using the sample id mapping
    df = (df
        .merge(sample_id_mapping, on='GEO_Accession (exp)', how='outer')
        .merge(supp_table[['SeqID', 'Library preparation batch', 'RNA extraction batch']]
                .rename(columns={'SeqID':'Sample_id'}), on='Sample_id', how='outer')
    )
    df.to_csv("../sra_metadata/roskams_metadata_preprocessed.csv", index=False)

    return df

def preprocess_ngo():
    csv_path = "../sra_metadata/ngo_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "ngo"
    df["dataset_batch"] = "ngo"
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "EDTA"
    #df["rna_extraction_kit"]="Norgen Plasma/Serum Circulating RNA and Exosomal Purification kit"
    #df["rna_extraction_kit_short"]="Norgen"
    #df["dnase"]="Baseline-ZERO"
    #df["library_prep_kit"]="SMARTer Stranded Pico v2"
    #df["library_prep_kit_short"]="SMARTer Pico v2"
    df.to_csv("../sra_metadata/ngo_metadata_preprocessed.csv", index=False)

    return df

def preprocess_ibarra():
    csv_path = "../sra_metadata/ibarra_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "ibarra"
    df["dataset_batch"] = "ibarra"
    df["biomaterial"] = df["tissue"].str.lower()
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-exome capture"
    df["plasma_tubes"] = df["biomaterial"].apply(lambda x: "EDTA" if x == "plasma" else "BD Vacutainer clotting tubes" if x == "serum" else "")
    #df["rna_extraction_kit"]="Qiagen QIAamp Circulating Nucleic Acid Kit"
    #df["rna_extraction_kit_short"]="Qiagen QIAamp"
    #df["dnase"]="Turbo DNAse"
    #df["library_prep_kit"]="Unspecified"
    #df["library_prep_kit_short"]="Unspecified"
    df.to_csv("../sra_metadata/ibarra_metadata_preprocessed.csv", index=False)

    return df

def preprocess_toden():
    csv_path = "../sra_metadata/toden_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "toden"
    df["dataset_batch"] = "toden"
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-exome capture"
    #df["plasma_tubes"] = "Unspecified"
    #df["rna_extraction_kit"]="Qiagen QIAamp Circulating Nucleic Acid Kit"
    #df["rna_extraction_kit_short"]="Qiagen QIAamp"
    #df["dnase"]="No"
    #df["library_prep_kit"]="Swift 2S kit"
    #df["library_prep_kit_short"]="Swift 2S kit"
    df.to_csv("../sra_metadata/toden_metadata_preprocessed.csv", index=False)

    return df

def preprocess_chalasani():
    csv_path = "../sra_metadata/chalasani_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "chalasani"
    df["dataset_batch"] = "chalasani"
    #df["biomaterial"] = "blood serum"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-exome capture"
    #df["plasma_tubes"] = "BD Vacutainer clotting tubes"
    #df["rna_extraction_kit"]="Qiagen QIAamp Circulating Nucleic Acid Kit"
    #df["rna_extraction_kit_short"]="Qiagen QIAamp"
    #df["dnase"]="Turbo DNAse"
    #df["library_prep_kit"]="Unspecified"
    #df["library_prep_kit_short"]="Unspecified"
    df.to_csv("../sra_metadata/chalasani_metadata_preprocessed.csv", index=False)

    return df

def preprocess_block():
    csv_path = "../sra_metadata/block_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "block"
    df["dataset_batch"] = np.where(
        abs(df["AvgSpotLen"] - 150) < abs(df["AvgSpotLen"] - 300),
        "block_1", "block_2"
    )
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "Unspecified"
    #df["rna_extraction_kit"]="Qiagen RNeasy Serum/Plasma Kit"
    #df["rna_extraction_kit_short"]="Qiagen RNeasy Plasma"
    #df["dnase"]="No"
    #df["library_prep_kit"]="SMARTer Stranded Pico v2"
    #df["library_prep_kit_short"]="SMARTer Pico v2"
    df.to_csv("../sra_metadata/block_metadata_preprocessed.csv", index=False)

    return df

def preprocess_rozowsky(): #this metadata does not come from SRA 
    csv_path = "../sra_metadata/rozowsky_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "rozowsky"
    df["dataset_batch"] = "rozowsky"
    #df["biomaterial"] = "tissue"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "NA"
    #df["rna_extraction_kit"]="Qiagen RNeasy Serum/Plasma Kit"
    #df["rna_extraction_kit_short"]="Qiagen RNeasy Plasma"
    #df["dnase"]="DNAse I"
    #df["library_prep_kit"]="custom"
    #df["library_prep_kit_short"]="custom"
    df.to_csv("../sra_metadata/rozowsky_metadata_preprocessed.csv", index=False)

    return df

def preprocess_tao():
    csv_path = "../sra_metadata/tao_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "tao"
    df["dataset_batch"] = "tao"
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "EDTA"
    #df["rna_extraction_kit"]="Norgen Plasma/Serum Circulating RNA and Exosomal Purification kit"
    #df["rna_extraction_kit_short"]="Norgen"
    #df["dnase"]="DNAse I"
    #df["library_prep_kit"]="SMARTer Stranded Pico v2"
    #df["library_prep_kit_short"]="SMARTer Pico v2"
    df.to_csv("../sra_metadata/tao_metadata_preprocessed.csv", index=False)

    return df

def preprocess_wei():
    csv_path = "../sra_metadata/wei_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "wei"
    df["dataset_batch"] = "wei"
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "DNA"
    df["library_selection"] = "whole-genome"
    #df["plasma_tubes"] = "EDTA"
    #df["rna_extraction_kit"]="NA"
    #df["rna_extraction_kit_short"]="NA"
    #df["dnase"]="No"
    #df["library_prep_kit"]="Truseq Nano DNA HT"
    #df["library_prep_kit_short"]="Truseq Nano DNA HT"
    df.to_csv("../sra_metadata/wei_metadata_preprocessed.csv", index=False)

    return df

def preprocess_moufarrej(): # we need to find a way to separate this dataset into the "sites" and into the cohorts. also decide if we want to merge by replicate
    csv_path = "../sra_metadata/moufarrej_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "moufarrej"
    df["dataset_batch"] = "moufarrej"
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "EDTA"
    #df["rna_extraction_kit"]="Norgen Plasma/Serum Circulating RNA and Exosomal Purification kit Slurry"
    #df["rna_extraction_kit_short"]="Norgen (slurry)"
    #df["dnase"]="Baseline-ZERO"
    #df["library_prep_kit"]="SMARTer Stranded Pico v2"
    #df["library_prep_kit_short"]="SMARTer Pico v2"
    df.to_csv("../sra_metadata/moufarrej_metadata_preprocessed.csv", index=False)

    return df

def preprocess_wang():
    csv_path = "../sra_metadata/wang_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "wang"
    df["dataset_batch"] = "wang"
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "EDTA"
    #df["rna_extraction_kit"]="Apostle MiniMaxTM High-Efficiency cfRNA Isolation Kit"
    #df["rna_extraction_kit_short"]="Apostle"
    #df["dnase"]="No"
    df["library_prep_kit"] = "SLiPiR-seq"
    df["library_prep_kit_short"] = "SLiPiR-seq"
    # There are 3 samples that were processed with a different library prep kit
    index = df[df['Sample Name'].str.startswith("NEB")].index
    df.loc[index, "library_prep_kit"] = "NEBNext Small RNA Library Prep Set"
    df.loc[index, "library_prep_kit_short"] = "NEBNext"
    df.to_csv("../sra_metadata/wang_metadata_preprocessed.csv", index=False)

    return df

def preprocess_giraldez():
    csv_path = "../sra_metadata/giraldez_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "giraldez"
    df["dataset_batch"] = np.where(
        abs(df["AvgSpotLen"] - 25) < abs(df["AvgSpotLen"] - 50),
        "giraldez_1", "giraldez_2"
    )
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "EDTA"
    #df["rna_extraction_kit"]="Qiagen miRNeasy Mini Kit"
    #df["rna_extraction_kit_short"]="Qiagen miRNeasy"
    #df["dnase"]="No"
    #df["library_prep_kit"]="Illumina TruSeq small RNA"
    #df["library_prep_kit_short"]="Illumina TruSeq small RNA"
    df.to_csv("../sra_metadata/giraldez_metadata_preprocessed.csv", index=False)

    return df

def preprocess_sun():
    csv_path = "../sra_metadata/sun_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "sun"
    df["biomaterial"] = df["tissue"].apply(lambda x: "blood plasma" if x == "plasma" else "blood serum" if x == "serum" else "")
    df["dataset_batch"] = np.where(df["biomaterial"] == "blood plasma", "sun_1", "sun_2")
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "EDTA"
    #df["rna_extraction_kit"]="MOF method"
    #df["rna_extraction_kit_short"]="MOF"
    #df["dnase"]="DNAse I"
    #df["library_prep_kit"]="SMARTer Stranded Total RNA-Seq"
    #df["library_prep_kit_short"]="SMARTer Total"
    df.to_csv("../sra_metadata/sun_metadata_preprocessed.csv", index=False)

    return df

def preprocess_decruyenaere():
    # Most relevant variables:
    col_names = [
        'sample_alias', 'sample_accession_id', 'biosample_id', 'run_accession_id', 'experiment_accession_id', 'study_accession_id', 'instrument_platform', 'instrument_model', 'library_layout', 'library_name', 'library_strategy', 'library_source', 'library_selection', 'run_file_type', 'design_description', 'description',
        'sample_title', 'phenotype', 'case_or_control', 'biological_sex', 'subject_id', 'cell_line', 'ENA-CHECKLIST',
        'organism_part', 'region'
    ]
    
    csv_path = "../sra_metadata/decruyenaere_metadata.csv"
    #df = pd.read_csv(csv_path)
    csv_path = "../sra_metadata/decruyenaere_metadata_ext.tsv"
    df = pd.read_csv(csv_path, sep='\t')[col_names]
    
    # only dataset from EGA-archive
    df = df.rename(columns={"sample_alias":"Run"})

    df["dataset_short_name"] = "decruyenaere"
    df["dataset_batch"] = "decruyenaere"
    df["biomaterial"] = "Blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome" # Overwrites existing column ( with values='cDNA_randomPriming')
    #df["plasma_tubes"] = "PAXgene blood ccfDNA"
    #df["rna_extraction_kit"]="Qiagen miRNeasy Serum/Plasma kit"
    #df["rna_extraction_kit_short"]="Qiagen miRNeasy Plasma"
    #df["dnase"]="HL-dsDNase"
    #df["library_prep_kit"]="SMARTer stranded Pico v3"
    #df["library_prep_kit_short"]="SMARTer Pico v3"
    df.to_csv("../sra_metadata/decruyenaere_metadata_preprocessed.csv", index=False)

    return df

def preprocess_reggiardo():
    csv_path = "../sra_metadata/reggiardo_metadata.csv"
    df = pd.read_csv(csv_path)

    df["dataset_short_name"] = "reggiardo"
    df["dataset_batch"] = "reggiardo"
    #df["biomaterial"] = "blood plasma"
    df["nucleic_acid_type"] = "total RNA"
    df["library_selection"] = "whole-transcriptome"
    #df["plasma_tubes"] = "K2EDTA"
    #df["rna_extraction_kit"]="Qiagen ExoRNeasy kit"
    #df["rna_extraction_kit_short"]="Qiagen ExoRNeasy"
    #df["dnase"]="No"
    #df["library_prep_kit"]="Takara SMART-Seq HT kit+Illumina Nextera XT DNA Prep "
    #df["library_prep_kit_short"]="SMART-Seq"
    df.to_csv("../sra_metadata/reggiardo_metadata_preprocessed.csv", index=False)

    return df


if __name__ == "__main__":
    preprocess_chen()
    preprocess_zhu()
    preprocess_roskams()
    preprocess_ngo()
    preprocess_ibarra()
    preprocess_toden()
    preprocess_chalasani()
    preprocess_block()
    preprocess_rozowsky()
    preprocess_tao()
    preprocess_wei()
    preprocess_moufarrej()
    preprocess_wang()
    preprocess_giraldez()
    preprocess_sun()
    preprocess_decruyenaere()
    preprocess_reggiardo()