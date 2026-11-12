import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import PyComplexHeatmap as pch
from sklearn.preprocessing import MinMaxScaler
import os
import matplotlib.cm as cm # For default colormaps
import matplotlib.colors as mcolors # For more detailed color mapping
import numpy as np # Import numpy for checking all zeros

def create_heatmap(# ~ BASIC INPUT PARAMETERS ~ #
                   df: pd.DataFrame,
                   value_columns: str or list,
                   metadata_columns: str or list,
                   output: str = "output/",
                   fig_width: float = 10,
                   rescale_values: bool = False,
                   min_relative_abundance: float = 0.0,
                   min_sample_percentage: float = 0.0,
                   row_cluster: bool = True,
                   remove_all_zero: bool = False,
                   metadata_colors_mapping: dict = None,
                   # ~ HEATMAP LAYOUT ~ #
                   row_height: float = 6,
                   show_legend: bool = True,
                   fig_height_inc: bool = None,
                   fig_width_inc: bool = None,
                   heatmap_label_fontsize: float = 10,
                   heatmap_value_name: str = "Value",
                   # ~ ANNOTATION LAYOUT  ~ #
                   annotation_fontsize:float = 10,
                   # ~ LEGEND ~ #
                   legend_vpad: float = -30,
                   legend_hpad: float = 0,
                   legend_gap: float = 3,
                   legend_border: bool = False,
                   legend_order: list = None,
                   legend_rename: dict = None):
    """
    Generates a heatmap from a pandas DataFrame, displaying multiple
    value columns along with multiple metadata columns.

    Args:
        # ~ BASIC INPUT PARAMETERS ~ #
        df (pd.DataFrame): The input pandas DataFrame.
        value_columns (list of str): A list of column names in `df` to use for the heatmap values.
        metadata_columns (list of str): A list of column names in `df` to use as
            metadata. These columns will form part of the x and y axes.
        output (str, optional): The base path for saving the heatmap. Defaults to "output/".
        fig_width (float, optional): figure wifth. Defaults to 10.
        rescale_values (bool, optional): If True, values in `value_columns` will be
            rescaled to a 0-1 range using Min-Max scaling before plotting.
            Defaults to False.
        min_relative_abundance (float, optional): The minimum relative abundance (between 0 and 1)
            a cell type must have in a sample to be considered. Defaults to 0.0.
        min_sample_percentage (float, optional): The minimum percentage of samples (between 0 and 1)
            in which a cell type must meet `min_relative_abundance` to be included. Defaults to 0.0.
        row_cluster (bool, optional): If True, rows will be clustered. Defaults to True.
        remove_all_zero (bool, optional): If True, value columns that contain only zero
            values across all rows will be removed. Defaults to False.
        metadata_colors_mapping (dict, optional): A dictionary where keys are metadata column names
            and values are dictionaries mapping specific metadata values to their colors.
            If a color is not provided for a value, it will be automatically assigned.
            Defaults to None.

        # ~ HEATMAP LAYOUT ~ #
        row_height (float, optional): The height of the rows (the scale can change based on other elements in the plot). Defaults to 6.
        show_legend (bool, optional): Show heatmap legend, Defaults to True.
        fig_height_inc (float, optional) overwrite the height of the heatmap in inches. Defaults to None
        fig_width_inc (float, optional) overwrite the width of the heatmap in inches. Defaults to None
        heatmap_label_fontsize (float, optional) size of the heatmap names. Defaults to 10.
        heatmap_value_name (str, optional) Name of heatmap value. Default to "Value"

        # ~ ANNOTATION LAYOUT  ~ #
        annotation_fontsize (float, optional) size of the annotation names. Defaults to 10.

        # ~ LEGEND ~ #
        legend_vpad (float, optional) the vertical pad of the legend items. Defaults to -30
        legend_hpad (float, optional) the horizontal pad of the legend items. Defaults to 0
        legend_gap (float, optional) gap between legend. Defaults to 3
        legend_border (bool, optional) Add border to legend box. Defaults to False
        legend_order (list, optional) Specify the order of legend. Defaults to None
        legend_rename (dict, optional) Dictionary to rename columns name for legend visualization. Defaults to None

    Returns:
        A single PyComplexHeatmap heatmap. Returns None if the input DataFrame is
            empty or if the specified columns are not found, or if no value columns
            remain after filtering.
    """

    # ~ Basic input validation ~ #
    if not isinstance(df, pd.DataFrame):
        print("Error: df must be a pandas DataFrame.")
        return None
    if df.empty:
        print("Error: DataFrame is empty.")
        return None
    if not isinstance(value_columns, list) or not all(isinstance(col, str) for col in value_columns):
        print("Error: value_columns must be a list of strings.")
        return None
    if not isinstance(metadata_columns, list) or not all(isinstance(col, str) for col in metadata_columns):
        print("Error: metadata_columns must be a list of strings.")
        return None
    if not isinstance(min_relative_abundance, (int, float)) or not (0.0 <= min_relative_abundance <= 1.0):
        print("Error: min_relative_abundance must be a float between 0 and 1.")
        return None
    if not isinstance(min_sample_percentage, (int, float)) or not (0.0 <= min_sample_percentage <= 1.0):
        print("Error: min_sample_percentage must be a float between 0 and 1.")
        return None
    if not isinstance(row_cluster, bool):
        print("Error: row_cluster must be a boolean.")
        return None
    if not isinstance(remove_all_zero, bool):
        print("Error: remove_all_zero must be a boolean.")
        return None
    if metadata_colors_mapping is not None and not isinstance(metadata_colors_mapping, dict):
        print("Error: metadata_colors_mapping must be a dictionary or None.")
        return None

    # ~ Check for the existence of columns in the DataFrame ~ #
    all_columns = value_columns + metadata_columns
    missing_columns_in_df = [col for col in all_columns if col not in df.columns]
    if missing_columns_in_df:
        print(f"Error: The following columns are not found in the DataFrame: {', '.join(missing_columns_in_df)}")
        return None

    # ~ Filter value_columns based on min_relative_abundance and min_sample_percentage ~ #
    filtered_value_columns = []
    num_samples = len(df)
    if num_samples == 0:
        print("Warning: DataFrame is empty, no filtering possible for value columns.")
        return None

    for col in value_columns:
        samples_meeting_abundance = (df[col] >= min_relative_abundance).sum()
        percentage_of_samples = samples_meeting_abundance / num_samples

        if percentage_of_samples >= min_sample_percentage:
            filtered_value_columns.append(col)
    
    if not filtered_value_columns:
        print(f"Warning: No value columns remain after filtering with min_relative_abundance={min_relative_abundance} and min_sample_percentage={min_sample_percentage}. Heatmap cannot be generated.")
        return None
    
    value_columns = filtered_value_columns # Update value_columns after initial filtering

    # ~ Filter out columns with all zero values if remove_all_zero is True ~ #
    if remove_all_zero:
        non_zero_value_columns = []
        removed_zero_columns = []
        for col in value_columns:
            # Check if all values in the column are zero
            # Ensure the column is numeric before checking for all zeros
            if pd.api.types.is_numeric_dtype(df[col]):
                if not (df[col] == 0).all():
                    non_zero_value_columns.append(col)
                else:
                    removed_zero_columns.append(col)
            else:
                # If it's not numeric, we can't reliably check for all zeros,
                # so we keep it. Or, you might choose to raise an error/warning
                # if a non-numeric column is in value_columns.
                non_zero_value_columns.append(col)

        if removed_zero_columns:
            print(f"Info: Removed the following value columns as they contain all zero values: {', '.join(removed_zero_columns)}")
        
        value_columns = non_zero_value_columns # Update value_columns after all-zero filtering

        if not value_columns:
            print("Warning: No value columns remain after removing all-zero columns. Heatmap cannot be generated.")
            return None


    # ~ Handling missing metadata for annotation and assigning colors ~ #
    df_for_metadata_annotation = df[metadata_columns].copy()
    
    colors_for_annotation = {} 

    default_cmaps = [plt.cm.tab10, plt.cm.Set3, plt.cm.Pastel1, plt.cm.Dark2, plt.cm.Accent]
    cmap_idx = 0

    for col_idx, col in enumerate(metadata_columns):
        # Fill actual NaN values with "Missing Information"
        df_for_metadata_annotation[col] = df_for_metadata_annotation[col].fillna("Missing Information")

        # Ensure the entire column is of string type.
        df_for_metadata_annotation[col] = df_for_metadata_annotation[col].astype(str)
        
        unique_values = df_for_metadata_annotation[col].unique()
        
        col_color_dict = {}
        
        # Always assign gray to "Missing Information"
        if "Missing Information" in unique_values:
            col_color_dict["Missing Information"] = "#BEBEBE"

        # Values to assign colors to (excluding "Missing Information")
        values_to_color = sorted([val for val in unique_values if val != "Missing Information"])
        
        # Apply user-defined colors first
        if metadata_colors_mapping and col in metadata_colors_mapping:
            for val in values_to_color:
                if val in metadata_colors_mapping[col]:
                    col_color_dict[val] = metadata_colors_mapping[col][val]
        
        # Assign automatic colors for values not in user mapping
        unassigned_values = [val for val in values_to_color if val not in col_color_dict]
        
        if len(unassigned_values) > 0:
            current_cmap = default_cmaps[cmap_idx % len(default_cmaps)]
            cmap_idx += 1 
            
            # Normalize for colormap if there's more than one unassigned value
            norm = mcolors.Normalize(vmin=0, vmax=len(unassigned_values) - 1)
            for i, val in enumerate(unassigned_values):
                color = current_cmap(norm(i))
                col_color_dict[val] = mcolors.to_hex(color) 
        
        colors_for_annotation[col] = col_color_dict
            
    if not metadata_columns:
        print("Warning: No metadata columns provided or remaining after processing. Heatmap annotation will be skipped.")
        col_ha = None
    else:
        valid_metadata_columns = [col for col in metadata_columns if not df_for_metadata_annotation[col].empty]
        
        if not valid_metadata_columns:
            print("Warning: All metadata columns are empty or became effectively empty after processing. Heatmap annotation will be skipped.")
            col_ha = None
        else:
            final_colors_for_annotation = {
                col: colors_for_annotation.get(col, {}) for col in valid_metadata_columns
            }

            print(final_colors_for_annotation)

            # ~ Create anno_simple for each column ~ #
            anno_simple_args = {}
            for col_name in df_for_metadata_annotation.columns:
                if col_name in final_colors_for_annotation:
                    # Initialize legend_kws
                    legend_kws = dict(
                        frameon=legend_border,
                        labelcolor='black',
                        color_text=False,
                        labelspacing=0.5
                    )
                    
                    # Add legend title if legend_rename is provided
                    if legend_rename is not None and col_name in legend_rename:
                        legend_kws["title"] = legend_rename[col_name]

                    anno_simple_args[col_name] = pch.anno_simple(
                        df_for_metadata_annotation[col_name],
                        colors=final_colors_for_annotation[col_name],
                        legend_kws=legend_kws
                    )

                else:
                    print(f"Warning: No colors defined for column '{col_name}'. Skipping annotation for this column.")

            col_ha = pch.HeatmapAnnotation(**anno_simple_args,
                                           plot=False,
                                           label_kws={'fontsize':annotation_fontsize},
                                           ticklabels_kws={'fontsize':annotation_fontsize},
                                           legend_gap=0.5,
                                           hgap=0,
                                           axis=1,
                                           legend_vpad=0)

    # ~ Apply rescaling if requested ~ #
    if rescale_values:
        df_for_heatmap_data = df[value_columns].copy()
        scaler = MinMaxScaler()
        df_for_heatmap_data[value_columns] = scaler.fit_transform(df_for_heatmap_data[value_columns])
    else:
        df_for_heatmap_data = df[value_columns]

    # ~ Calculate dynamic figure height based on number of rows (value_columns) ~ #
    base_height = row_height  # Minimum height
    height_per_row = 0.1 # Height contribution per row (value column)
    dynamic_height = max(base_height, len(value_columns) * height_per_row)

    plt.figure(figsize=(fig_width, dynamic_height))

    # ~ Initialize ClusterMapPlotter ~ #
    cm = pch.ClusterMapPlotter(data=df_for_heatmap_data.transpose(),
                               top_annotation=col_ha,
                               col_cluster=True,
                               row_cluster=row_cluster,
                               yticklabels=True,
                               col_split_gap=0.5,
                               row_split_gap=0.8,
                               label=heatmap_value_name,
                               row_dendrogram=True,
                               col_dendrogram=True,
                               show_rownames=True,
                               show_colnames=False,
                               tree_kws={'row_cmap': 'Set1'},
                               verbose=1,
                               plot_legend=show_legend,
                               legend_gap=legend_gap,
                               legend_vpad=legend_vpad,
                               legend_hpad=legend_hpad,
                               yticklabels_kws={'labelsize': heatmap_label_fontsize},
                               cmap='viridis',
                               legend_order= metadata_columns+[heatmap_value_name] if legend_order is None else legend_order+[heatmap_value_name]
)
    
    if show_legend:
        for i in range(len(cm.cbars)):
            print(i)
            print(type(cm.cbars[i]))
            if isinstance(cm.cbars[i], matplotlib.colorbar.Colorbar):
                print(type(cm.cbars[i]))
                cm.cbars[i].set_label('Cell type\nfraction',fontsize=5)
                cm.cbars[i].set_ticklabels([0,0.5,1],fontsize=5)
                cm.cbars[i].ax.tick_params(axis='y',  width=0.25)
                cm.cbars[i].outline.set_linewidth(0.25)
            plt.rcParams["legend.title_fontsize"] = 5
            plt.rcParams["legend.fontsize"] = 5               
        plt.rcParams["legend.title_fontsize"] = 5
        plt.rcParams["legend.fontsize"] = 5

    current_fig = plt.gcf()
    current_fig_width, current_figure_height = current_fig.figure.get_size_inches()

    new_width = current_fig_width
    new_height = current_figure_height

    if fig_height_inc:
        new_height = fig_height_inc
    if fig_width_inc:
        new_width = fig_width_inc

    current_fig.set_size_inches(new_width, new_height)
    
    # ~ Save plots ~ #
    # Create the directory if it doesn't exist
    dir_name = os.path.dirname(output)
    if dir_name and not os.path.exists(dir_name):
        os.makedirs(dir_name)
    
    # Save the plot (PDF)
    filename_pdf = output + ".pdf"
    plt.savefig(filename_pdf, format='pdf', bbox_inches='tight', dpi=600)
    print(f"Heatmap saved to {filename_pdf}")

    # Save the plot (PNG)
    filename_png = output + ".png"
    plt.savefig(filename_png, format='png', bbox_inches='tight', dpi=600)
    print(f"Heatmap saved to {filename_png}")

    # Save the plot (SVG)
    filename_svg = output + ".svg"
    plt.rcParams["svg.fonttype"] = "none"
    plt.savefig(filename_svg, format='svg', bbox_inches='tight', dpi=600)
    print(f"Heatmap saved to {filename_svg}")

    plt.show()

    final_data = df_for_metadata_annotation[valid_metadata_columns]
    final_color = final_colors_for_annotation

    return cm, final_data, final_color

