#!/usr/bin/env python3


import pandas as pd
from preprocess_metadata_functions import *
from sra_columns_mapping import rename_columns_and_values
pd.set_option('display.max_columns', None)

def main():
    # Read plasma cfRNA-seq workflow comparison
    csv_path = "../sra_metadata/dataset_metadata.tsv"
    dataset_metadata = pd.read_csv(csv_path, sep = "\t")

    chen = preprocess_chen()
    zhu = preprocess_zhu()
    roskams = preprocess_roskams()
    ngo = preprocess_ngo()
    ibarra = preprocess_ibarra()
    toden = preprocess_toden()
    chalasani = preprocess_chalasani()
    block = preprocess_block()
    rozowsky = preprocess_rozowsky()
    tao = preprocess_tao()
    wei = preprocess_wei()
    moufarrej = preprocess_moufarrej()
    wang = preprocess_wang()
    giraldez = preprocess_giraldez()
    sun = preprocess_sun()
    decruyenaere = preprocess_decruyenaere()
    reggiardo = preprocess_reggiardo()

    chen = chen.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    zhu = zhu.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    roskams = roskams.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    ngo = ngo.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    ibarra = ibarra.merge(dataset_metadata[['dataset_short_name', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    toden = toden.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    chalasani = chalasani.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    block = block.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    rozowsky = rozowsky.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    tao = tao.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    wei = wei.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    moufarrej = moufarrej.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    wang = wang.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    giraldez = giraldez.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    sun = sun.merge(dataset_metadata[['dataset_short_name', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    decruyenaere = decruyenaere.merge(dataset_metadata[['dataset_short_name', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')
    reggiardo = reggiardo.merge(dataset_metadata[['dataset_short_name','Biomaterial', 'Plasma tubes', 'DNAse','RNA extraction kit','RNA extraction kit (short name)','Library prep kit','Library prep kit (short name)']], on='dataset_short_name', how='left')

    dfs = [chen, zhu, roskams, ngo, ibarra, toden, chalasani, block, rozowsky, tao, wei, moufarrej, wang, giraldez, sun, decruyenaere, reggiardo]

    merged_metadata = pd.concat(dfs, axis=0, join="outer", ignore_index=True)

    merged_metadata = rename_columns_and_values(merged_metadata)

    merged_metadata.to_csv("../tables/cfRNA-meta_per_sample_metadata.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()