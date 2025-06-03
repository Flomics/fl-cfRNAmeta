import numpy as np

def rename_columns_and_values(df):
    # disease, disease_state, subject_status
    mapping = {
        'Esophagus Cancer patient':'Esophagus cancer',
        'Lung Cancer patient':'Lung cancer',
        'Liver Cancer patient':'Liver cancer',
        'Stomach Cancer patient':'Stomach cancer',
        'Colorectal Cancer patient':'colorectal cancer',
        'Healthy donor':'Control',
        'NASH':'Nonalcoholic steatohepatitis',
        'NAFLD':'Nonalcoholic fatty liver disease',
        'Healthy':'Control',
        'Control':'Control',
        'Hepatoma':'Liver cancer',
        'Cirrhosis':'Cirrhosis',
        'pre-eclampsia':'Pre-eclampsia',
        'severe pre-eclampsia':'Pre-eclampsia',
        'liver cancer patient':'Liver cancer',
        'healthy donor':'Control',
        'normal':'Control',
        'HCC':'Liver cancer',
        'CHB':'Chronic hepatitis B',
        'pancreatic cancer patient':'Pancreatic cancer',
        'normal healthy donor':'Control'
    }
    df['disease'] = df['disease'].fillna(df['disease_state']).fillna(df['subject_status'])
    df = df.drop(['disease_state', 'subject_status'], axis=1)
    df['disease'] = df['disease'].replace(mapping)
    print("df['disease'].unique():\n", df['disease'].unique())

    # Instrument
    mapping = {
        'HiSeq X Ten':'Ilumina HiSeq X Ten',
        'Illumina HiSeq 4000':'Illumina HiSeq 4000',
        'NextSeq 500':'Illumina NextSeq 500',
        'Illumina HiSeq 2500':'Illumina HiSeq 2500',
        'Illumina HiSeq 2000':'Illumina HiSeq 2000',
        'Illumina NovaSeq 6000':'Illumina NovaSeq 6000',
        'MinION':'MinION'
    }
    df['instrument'] = df['instrument'].replace(mapping)
    print("df['instrument'].unique():\n", df['instrument'].unique())

    # sex
    mapping = {
        'missing':np.nan,
        'pooled male and female':np.nan
    }
    df['sex'] = df['sex'].fillna(df['gender'])
    df['sex'] = df['sex'].replace(mapping)
    df = df.drop(['gender'], axis=1)
    print("df['sex'].unique():\n", df['sex'].unique())

    # biomaterial, sample_type, source_name
    # The first column name in the mapping dic will be considered as the main column
    # The second, third, ... columns will be copied over in the first column.
    mapping = {
        'biomaterial': {
            'serum':'serum',
            'plasma':'plasma',
            'buffy coat':'buffy coat',
            'blood serum':'serum',
            'blood plasma':'plasma',
            'Blood plasma':'plasma',
            },
        'sample_type': {
            'cell culture':'cell culture',
            'plasma':'plasma',
            'serum':'serum',
            'rna':'tissue',
            'tissue':'tissue'},
        'source_name': {
            'plasma':'plasma',
            "liver cancer patient's plasma":'plasma',
            "healthy donor's plasma":'plasma',
            'human non-cancer donor plasma':'plasma',
            'human multiple myeloma plasma':'plasma',
            'human mgus plasma':'plasma',
            'human liver cancer plasma':'plasma',
            'human liver cirrhosis plasma':'plasma',
            'pbmc':'PBMC',
            'tumor':'tumor tissue',
            'normal adjacent tumor tissue':'normal adjacent tumor tissue',
            'synthetic srna equimolar pool (synthrna476v1)':np.nan,
            'platelet-poor blood plasma':'plasma',
            'serum':'serum',
            'blood plasma':'plasma'
        }
    }
    # For the sun dataset, both the "biomaterial" and the "source_name" columns contain similar information about the biomaterial. We simply keep the main column and discard the "source_name" column.
    cols = list(mapping.keys())
    main_col = cols[0]
    for col in cols:
        df[col] = df[col].str.lower().replace(mapping[col])
    df[main_col] = df[main_col].fillna(df[cols[1]]).fillna(df[cols[2]])
    df = df.drop(cols[1:], axis=1)
    print(f"df['{main_col}'].unique():\n", df[main_col].unique())

    # Drop empty columns
    df = df.dropna(how='all', axis=1)

    return df
