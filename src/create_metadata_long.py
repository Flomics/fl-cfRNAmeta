#!/usr/bin/env python3


import pandas as pd
from preprocess_metadata_functions import *
from sra_columns_mapping import rename_columns_and_values
pd.set_option('display.max_columns', None)

def main():


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