import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

def create_stacked_barplot(
                           # ~ BASIC INPUT PARAMETERS ~ #
                           dataset: pd.DataFrame,
                           meta_column: str,
                           value_columns: list,
                           value_order: str = 'default',
                           meta_order: list = None,
                           output: str = None,
                           fig_width: float = 10,
                           fig_height: float = 7,
                           ax_width: float = None,
                           aspect: float = None,
                           dpi: float = 600,
                           # ~ SCALING OPTIONS ~ #
                           scaling: str = 'none',
                           # ~ NORMALIZATION OPTIONS ~ #
                           normalize_data: bool = False,
                           # ~ EDIT FOCUS VALUES FEATURES ~ #
                           focus_value: list = None,
                           collapse_focus_values_as: str = None,
                           collapsed_color: str = None,
                           color_map: dict = None,
                           # ~ ERROR BARS PARAMETERS ~ #
                           add_error_bars: bool = False,
                           # ~ CONNECTING SHADES PARAMETERS ~ #
                           add_connecting_shades: bool = False,
                           connecting_shades_alpha: float = 0.15,
                           # ~ CATEGORY BORDER PARAMETERS ~ #
                           add_category_border: bool = False,
                           category_border_width: float = 0.5,
                           # ~ BOXES PARAMETERS ~ #
                           boxes_column: str = None,
                           boxes_color_map: dict = None,
                           boxes_y_position: float = 1.05,
                           boxes_height: float = 0.1,
                           boxes_width: float = 1,
                           boxes_borderwidth: float = 0.5,
                           boxes_legend: bool = True,
                           boxes_legend_pos: str = None,
                           boxes_legend_title: str = None, 
                           boxes_legend_fontsize: float = 15,
                           boxes_legend_y_pos: float = 1,
                           boxes_legend_x_pos: float = 0.5,
                           # ~ METADATA GROUPING ~ #
                           group_by_column: str = None,
                           group_spacing: float = 0,
                           group_label_rotation: float = 0,
                           group_label_fontsize: int = 12,
                           group_label_y_offset: float = 0.0,
                           group_bracket_linewidth: float = 1.0,
                           group_bracket_vertical_line_length: float = 0.05,
                           group_bracket_horizontal_line_length: float = 0.5,
                           group_position: str = 'bottom',
                           show_group_label: bool = True,
                           # ~ FONTSIZE AND LAYOUT PARAMETERS ~ #
                           # ~ axis Label ~ #
                           show_xlabel: bool = True, 
                           xlabel: str = None,                                           
                           xlabel_fontsize: float = 12,
                           show_ylabel: bool = True,
                           ylabel: str = None,
                           ylabel_fontsize: float = 12,
                           # ~ axis tick label ~ #
                           x_tick_label_fontsize: float = 10,
                           y_tick_label_fontsize: float = 10,
                           x_tick_label_rotation: float = 90,
                           x_ticks_label_pad: float = 5,
                           # ~ axis tick ~ #
                           hide_bottom_tick: bool = False,
                           hide_left_tick: bool = False,
                           yticks: list = None,
                           # ~ title ~ #
                           show_title: bool = True,
                           title: str = None,
                           title_fontsize: float = 14,
                           title_pad: float = 10,
                           # ~ legend ~ #
                           legend_title: str = None,
                           legend_title_fontsize: float = 12,
                           legend_y_pos: float = 0.5,
                           legend_x_pos: float = 1.5,
                           legend_fontsize: float = 10,
                           legend_position: str = "center left",
                           # ~ bar ~ #
                           bar_width: float = 0.8,
                           # ~ spine ~ #
                           hide_top_spine: bool = False,
                           hide_right_spine: bool = False,
                           hide_bottom_spine: bool = False,
                           hide_left_spine: bool = False,
                           # ~ others ~ #
                           y_upper_pad: float = 0.05):
    """
    Generates a stacked bar plot with specified columns, collapsing, scaling, and optional error bars
    representing the raw standard deviation.

    Args:

        # ~ BASIC INPUT PARAMETERS ~ #

        dataset (pd.DataFrame): The input pandas DataFrame.
        meta_column (str): The column to be used on the X-axis (only one accepted).
        value_columns (list): A list of columns that represent the values to represent on the Y-axis of the stacked bar.
                               If {value_order} is 'default', the list will be used to order the categories
        value_order (str, optional): Determines the order of the values.
                                     'default': {value_column} list is used to fetch the order
                                     'median_descending': Median value of category across all {meta_column} is used to descending order.
                                     'median_ascending': Median value of category across all {meta_column} is used to ascending order.
        meta_order (list, optional): Determines the order of {meta_column} on X-axis. If None, default order provided in {meta_column} is used.
                                        If grouping is used, make sure {meta_order} do not create conflicts with consistent grouping.
        output (str, optional): The path plus filename (without extension) to save the plot.
                                 If provided, the plot will be saved as a PDF and PNG.
        fig_width (float, optional): The width of the figure in inches. Defaults to 10.
        fig_height (float, optional): The height of the figure in inches. Defaults to 7.
        ax_width (float, optional): Define absolute size of axes width. Defaults to None
        aspect (float, optional): If {ax_width} is specified, define the ratio between {ax_width} to axes height. Defaults to None
        dpi (float, optional): Dots per inch. Defaults to 600

        # ~ SCALING OPTIONS ~ #

        scaling (str, optional): Determines how values are aggregated for each bar.
                                 'none': Values are summed for each meta_column group.
                                 'median': The median of each category's values is used for each meta_column group.
                                 'mean': The mean of each category's values is used for each meta_column group.
                                 Defaults to 'none'.

        # ~ NORMALIZATION OPTIONS ~ #

        normalize_data (bool, optional): If True, the stacked bars will represent proportions (summing to 1).
                                         If False, the raw counts will be plotted. Defaults to True.

        # ~ EDIT FOCUS VALUES FEATURES ~ #

        focus_value (list, optional): A subset of value_columns to focus on, gouping the excluded ones in a "other" category.
                                       If {collapse_focus_values_as} is None, columns in focus_value
                                       will be plotted individually.
                                       If {collapse_focus_values_as} is provided, columns in focus_value
                                       will be summed into a single category.
                                       If None or an empty list, all value_columns are plotted individually
                                       without an "others" category.
        collapse_focus_values_as (str, optional): If provided and {focus_value} is not empty,
                                                   all columns listed in {focus_value} will be summed
                                                   and plotted as a single category with this name.
                                                   Cannot be named 'others'.
        collapsed_color (str, optional): The color to use for the collapsed category specified by
                                          {collapse_focus_values_as}. If None, a default color will be used.
        color_map (dict, optional): A dictionary mapping category names to specific colors.
                                            Colors can be matplotlib color names (strings) or RGBA tuples (e.e., (0.1, 0.2, 0.3, 1.0)).
                                            If a category is not in the map, a default color will be used.

        # ~ ERROR BARS PARAMETERS ~ #

        add_error_bars (bool, optional): If True, error bars representing the raw standard deviation
                                          of the counts will be added to the bars. Defaults to False.

        # ~ CONNECTING SHADES PARAMETERS ~ #

        add_connecting_shades (bool, optional): If True, adds semi-transparent shaded regions connecting
                                                 the same categories across different stacked bars. Defaults to False.
        connecting_shades_alpha (float, optional): Transparency level for connecting shades. Defaults to 0.15.

        # ~ CATEGORY BORDER PARAMETERS ~ #

        add_category_border (bool, optional): If True, adds a thin black line around each category in the stacked bars.
                                               Defaults to True.
        category_border_width (float, optional): The width of the black line around categories if {add_category_border} is True.
                                                 Defaults to 0.5.

        # ~ BOXES PARAMETERS ~ #
        
        #TO DO: add more than one line in top_box
        boxes_column (str, optional): Name of the column in the provided dataset to be reppresented as a box above bars.
        boxes_color_map (dict, optional): Dictionary mapping unique values from {boxes_column} to colors.
                                                If None, defaults colors will be used.
        boxes_y_position (float, optional): Value to select the box position in the y axis. defaults to 1.05
                                                    The value is proportional to the Y-axis scale.
                                                    It might be useful to harmonize this value with the {y_upper_pad} to have a better visualization.
        boxes_height (float, optional) : The height of the top boxes. Defaults to 0.1
                                                The value is proportional to the Y-axis scale.
                                                It might be useful to harmonize this value with the {y_upper_pad} to have a better visualization.                                                
        boxes_width (float, optional): The width of the top boxes. Defaults to 1
        boxes_borderwidth (float, optional): linewidth for boxes. Defaults to 0.5
        boxes_legend (bool, optional): Show top_boxes position. Defaults True
        boxes_legend_pos (str, optional): Custom position for boxes legend.
                                            "bottom": will show a one line legend in the bottom part of the figure.
                                            Defaults to None
        boxes_legend_title (str, ptional): Boxes legend Title. Default to {boxes_column} 
        boxes_legend_fontsize (float, optional): font size of boxes legend element text. Defaults to 15
        boxes_legend_y_pos (float, optional): Position of top_box legend on Y-axis
        boxes_legend_x_pos (float, optional): Position of top_box legend on X-axis

        # ~ METADATA GROUPING ~ #

        group_by_column (str, optional): The name of a column in the provided dataset to use for grouping {meta_column} values.
                                          If provided, bars will be grouped based on this column's values,
                                          and group labels with brackets will be added below the x-axis. Defaults to None.
        group_spacing (float, optional): The extra space to add between groups of bars if {group_by_column} is used.
                                          Defaults to 0.
        group_label_rotation (float, optional): Rotation angle for the group labels in degrees. Defaults to 0 (the one below the brackets).
        group_label_fontsize (int, optional): Font size for the group labels. Defaults to 12.
        group_label_y_offset (float, optional): Vertical offset for the group labels and brackets. A negative value
                                                 moves them further down from the x-axis. Defaults to 0.0.
        group_bracket_linewidth (float, optional): The line width for the square brackets when grouping is enabled.
                                                    Defaults to 1.0.
        group_bracket_vertical_line_length (float, optional): Controls the length of the vertical lines of the square brackets.
                                                                This value will be added to or subtracted from the base
                                                                {bracket_y_level} to define the extent of the vertical lines.
                                                                Defaults to 0.05.
        group_bracket_horizontal_line_length (float, optional): Controls the length of the horizontal lines of the square brackets.
                                                                0.5 mean brackets end will overlap with the next bracket start
                                                                Defaults to 0.5.
        group_position (str, optional): Determines the position of group labels and brackets. default to 'bottom'
                                        'bottom': Below the x-axis. x_tick_label_rotation to 90 suggested for better visualization
                                        'middle': between x-ticks and x-label
                                        'top': Above the plot, near the title.
        show_group_label (bool, optional): If True, the group labels will be displayed. Defaults to True.

        # ~ FONTSIZE AND LAYOUT PARAMETERS ~ #
        # ~ axis Label ~ #
        show_xlabel (bool, optional): If True, the x-axis label will be displayed. Defaults to True.
        xlabel (str, optional): Custom label for the x-axis. If provided, overrides default. Defaults to None.
        xlabel_fontsize (float, optional): The font size for the x-axis label. Defaults to 12.
        show_ylabel (bool, optional): If True, the y-axis label will be displayed. Defaults to True.
        ylabel (str, optional): Custom label for the y-axis. If provided, overrides default. Defaults to None.
        ylabel_fontsize (float, optional): The font size for the y-axis label. Defaults to 12.
        # ~ axis tick label ~ #
        x_tick_label_fontsize (float, optional): The font size for the x-axis tick labels. Defaults to 10.
        y_tick_label_fontsize (float, optional): The font size for the y-axis tick labels. Defaults to 10.
        x_tick_label_rotation (float, optional): Rotation angle for the x label in degrees. Defaults to 90.
        x_ticks_label_pad (float, optional): Distance of x-tick labels from the plot. Defaults to 5.
        # ~ axis tick ~ #
        hide_bottom_tick (bool, optional): Hide bottom ticks. Defaults to False
        hide_left_tick (bool, optional): Hide left ticks. Defaults to False
        yticks (list, optional): list of Y-axis ticks. Defaults to None.
        # ~ title ~ #
        show_title (bool, optional): If True, the plot title will be displayed. Defaults to True.
        title (str, optional): Custom title, if provided, overrides default. Defaults to None.
        title_fontsize (float, optional): The font size for the plot title. Defaults to 14.
        title_pad (float, optional): Distance between title and plot. Default to 10
        # ~ legend ~ #
        legend_title (str, optional): Custom name for legend. Default is None
        legend_title_fontsize (float, optional): legend title fontsize. Default to 12.
        legend_y_pos (float, optional): Position of legend in Y-axis. Default 0.5
        legend_x_pos (float, optional): Position of legend in X-axis. Default 1.5
        legend_fontsize (float, optional): legend text fontsize. Default to 10.
        legend_position (str, optional): position of the main legend, based on matplot and custom values.
                                            Matplot values, 'best','upper right','upper left','lower left','lower right','right','center left','center right','lower center','upper center','center'
                                            Custom values, "custom bottom"
                                            Defaults to 'center left'.
        # ~ bar ~ #
        bar_width (float, optional): width of the bars. Default to 0.8 
        # ~ spine ~ #
        hide_top_spine (bool, optional): Hide the top spine. Default to False
        hide_right_spine (bool, optional): Hide the right spine. Default to False
        hide_bottom_spine (bool, optional): Hide the bottom spine. Default to False
        hide_left_spine (bool, optional): Hide the left spine. Default to False
        # ~ others ~ #
        y_upper_pad (float, optional): Size of the space between the end of the bar and the upper plot margin. Default to 0.05 (i.e. 5%)
    """

    # ~ Input Validation ~ #
    if not isinstance(dataset, pd.DataFrame):
        raise TypeError("Input dataset must be a pandas DataFrame.")
    if meta_column not in dataset.columns:
        raise ValueError(f"Column '{meta_column}' not found in the dataset.")
    if not all(col in dataset.columns for col in value_columns):
        raise ValueError("Not all columns in \{value_columns\} found in the dataset.")

    # Ensure focus_value is a list, defaulting to empty if None
    if focus_value is None:
        focus_value = []

    if not all(col in value_columns for col in focus_value):
        raise ValueError("All columns in \{focus_value\} must also be present in \{value_columns\}.")

    # Validation for collapse_focus_values_as
    if collapse_focus_values_as is not None:
        if not isinstance(collapse_focus_values_as, str):
            raise TypeError("\{collapse_focus_values_as\} must be a string if provided.")
        if not focus_value:
            raise ValueError("\{focus_value' cannot be empty if \{collapse_focus_values_as\} is provided, as there would be nothing to collapse.")
        if collapse_focus_values_as == 'others':
            raise ValueError("\{collapse_focus_values_as\} cannot be named 'others' as it conflicts with an internal category name (non-focus categories).")
        # Check for name conflict with existing columns not being melted
        if collapse_focus_values_as in dataset.columns and \
           collapse_focus_values_as not in value_columns and \
           collapse_focus_values_as != meta_column:
            raise ValueError(f"The chosen name '{collapse_focus_values_as}' for the collapsed category already exists as a column in the original dataset and is not part of the 'value_columns' or 'meta_column'. Please choose a different name.")

    # Validate collapsed_color if provided
    if collapsed_color is not None and not isinstance(collapsed_color, str):
        raise TypeError("\{collapsed_color\} must be a string representing a color if provided.")

    # Validate color_map if provided
    if color_map is not None and not isinstance(color_map, dict):
        raise TypeError("\{color_map\} must be a dictionary if provided.")
    if color_map is None: # Ensure it's a dict to avoid errors later
        color_map = {}
    
    # validate ax_width and aspect
    if (ax_width is None and aspect is not None) or (ax_width is not None and aspect is None):
        raise ValueError("Both 'ax_width' and 'aspect' must be specified together")

    # Validate group_by_column if provided
    if group_by_column:
        if not isinstance(group_by_column, str):
            raise TypeError("\{group_by_column\} must be a string if provided.")
        if group_by_column not in dataset.columns:
            raise ValueError(f"Column '{group_by_column}' not found in the dataset.")
        # Ensure that each meta_column value consistently maps to one group_by_column value
        group_consistency_check = dataset[[meta_column, group_by_column]].drop_duplicates()
        if group_consistency_check.duplicated(subset=[meta_column]).any():
            conflicting_meta_values = group_consistency_check[group_consistency_check.duplicated(subset=[meta_column])][meta_column].tolist()
            raise ValueError(f"Each value in '{meta_column}' must correspond to a single value in '{group_by_column}'. "
                             f"Conflicting '{meta_column}' values: {conflicting_meta_values}")
    
    if group_position not in ['bottom', 'middle', 'top']:
        raise ValueError("Invalid value for \{group_position\}. Choose from 'bottom', 'middle', or 'top'.")

    # Validate scaling parameter
    if scaling not in ['none', 'median', 'mean']:
        raise ValueError("Invalid value for \{scaling\}. Choose from 'none', 'median', or 'mean'.")

    # Validate value_order parameter
    if value_order not in ['default', 'median_descending', 'median_ascending']:
        raise ValueError("Invalid value for \{value_order\}. Choose from 'default', 'median_descending', or 'median_ascending'.")

    # Validate meta_order parameter
    if meta_order is not None and not isinstance(meta_order, list):
        raise TypeError("Meta_order must be a list of meta_column values.")
    if meta_order is not None:
        # Check if all meta_column values in the dataset are present in meta_order
        dataset_meta_values = dataset[meta_column].unique()
        if not set(dataset_meta_values).issubset(set(meta_order)):
            missing_values = set(dataset_meta_values) - set(meta_order)
            raise ValueError(f"Not all unique values from '{meta_column}' are present in {{meta_order}}. Missing: {list(missing_values)}")
        # Check for duplicates in meta_order
        if len(meta_order) != len(set(meta_order)):
            raise ValueError("Meta_order list must not contain duplicate values.")
        
    # validate top boxes column parameters
    if boxes_column:
        if not isinstance(boxes_column, str):
            raise TypeError("\{boxes_column\} must be a string") 
        if boxes_column not in dataset.columns:
            raise ValueError(f"Column {{boxes_column}} not found in the dataset")
        # Ensure consistency between metadata and top boxes column values
        top_boxes_consistency_check = dataset[[meta_column,boxes_column]].drop_duplicates()
        if top_boxes_consistency_check.duplicated(subset=[meta_column]).any():
            conflicting_meta_values = top_boxes_consistency_check[top_boxes_consistency_check.duplicated(subset=[meta_column])][meta_column].to_list()
            raise ValueError(f"Each value in {{meta_column}} should correspond to a unique value in {{boxes_column}}, conflicting {{meta_column}} values: {conflicting_meta_values}")
        # Validate {boxes_color_map}
        if boxes_color_map is not None and not isinstance(boxes_color_map, dict):
            raise TypeError("\{boxes_color_map\} must be a dictionary if provided")
        if boxes_color_map is None:
            boxes_color_map = {}

    # ~ Data Preparation, handle focus values scenarios ~ #
    # Create a copy of the relevant columns to avoid modifying the original DataFrame
    df_plot = dataset[[meta_column] + value_columns].copy()
    if group_by_column:
        df_plot = df_plot.merge(dataset[[meta_column, group_by_column]].drop_duplicates(), on=meta_column, how='left')
    if boxes_column:
        df_plot = df_plot.merge(dataset[[meta_column, boxes_column]].drop_duplicates(), on=meta_column, how='left')


    # List for columns that will be plotted on Y-axis
    cols_to_plot_y = []
    # Dictionary to map category names to their assigned colors
    category_colors_internal_map = {}

    # Use default matplotlib color and allow repeating the cicle
    prop_cycle = plt.rcParams['axes.prop_cycle']
    default_colors = prop_cycle.by_key()['color']
    color_idx = 0 # Index to cycle through default colors

    # Scenario 1: Collapse focus_value into a single new category
    if collapse_focus_values_as is not None and focus_value:
        df_plot[collapse_focus_values_as] = df_plot[focus_value].sum(axis=1)
        cols_to_plot_y.append(collapse_focus_values_as)
        # Assign a color for the collapsed focus_value category
        if collapsed_color: # Prioritize specific collapsed_color if provided
            category_colors_internal_map[collapse_focus_values_as] = collapsed_color
        elif collapse_focus_values_as in color_map: # check general color_map
            category_colors_internal_map[collapse_focus_values_as] = color_map[collapse_focus_values_as]
        else: # Otherwise use default cycle
            category_colors_internal_map[collapse_focus_values_as] = default_colors[color_idx % len(default_colors)]
        color_idx += 1

        # Determine if there are other columns to collapse into 'others'
        remaining_value_columns_for_others = [col for col in value_columns if col not in focus_value]
        if remaining_value_columns_for_others:
            df_plot['others'] = df_plot[remaining_value_columns_for_others].sum(axis=1)
            cols_to_plot_y.append('others')
            if 'others' in color_map:
                category_colors_internal_map['others'] = color_map['others']
            else:
                category_colors_internal_map['others'] = 'gray' # gray color for 'others'

    # Scenario 2: Plot focus_value individually and collapse remaining to 'others'
    elif focus_value:
        cols_to_plot_y.extend(focus_value)
        # Assign individual colors for each focus_value column
        for col in focus_value:
            # Use color from color_map if available, otherwise use default cycle
            if col in color_map:
                category_colors_internal_map[col] = color_map[col]
            else:
                category_colors_internal_map[col] = default_colors[color_idx % len(default_colors)]
                color_idx += 1

        # Identify columns to be collapsed into "others"
        others_cols = [col for col in value_columns if col not in focus_value]
        if others_cols:
            df_plot['others'] = df_plot[others_cols].sum(axis=1)
            cols_to_plot_y.append('others')
            if 'others' in color_map:
                category_colors_internal_map['others'] = color_map['others']
            else:
                category_colors_internal_map['others'] = 'gray' # gray color for 'others'

    # Scenario 3: No focus_value specified, plot all value_columns individually
    else:
        cols_to_plot_y = value_columns
        # Assign individual colors for all value_columns
        for col in value_columns:
            if col in color_map:
                category_colors_internal_map[col] = color_map[col]
            else:
                category_colors_internal_map[col] = default_colors[color_idx % len(default_colors)]
                color_idx += 1

    # Select the columns for the final DataFrame used for grouping
    df_plot_final = df_plot[[meta_column] + cols_to_plot_y]
    if group_by_column:
        df_plot_final = df_plot_final.merge(dataset[[meta_column, group_by_column]].drop_duplicates(), on=meta_column, how='left')
    if boxes_column:
        df_plot_final = df_plot_final.merge(dataset[[meta_column, boxes_column]].drop_duplicates(), on=meta_column, how='left')


    # ~ Calculate error bars (if requested) before any scaling/normalization ~ #
    # The standard deviation should always be calculated from the original "counts" per meta_column,
    grouped_std = None
    if add_error_bars:
        # Calculate standard deviation for each component within each meta_column group.
        grouped_std = df_plot_final.groupby(meta_column)[cols_to_plot_y].std()
        grouped_std = grouped_std.fillna(0) # Fill NaN standard deviations with 0

    dataframe_before_scaling = df_plot_final

    # ~ Apply Scaling (Sum, Median, or Mean) ~ #
    if scaling == 'none':
        grouped_df = df_plot_final.groupby(meta_column)[cols_to_plot_y].sum()
    elif scaling == 'median':
        grouped_df = df_plot_final.groupby(meta_column)[cols_to_plot_y].median()
    elif scaling == 'mean':
        grouped_df = df_plot_final.groupby(meta_column)[cols_to_plot_y].mean()

    dataframe_after_scaling = grouped_df

    # ~ Apply Normalization ~ #
    if normalize_data:
        # Calculate row sums for normalization. Handle cases where a row sum might be zero to avoid division by zero.
        row_sums = grouped_df.sum(axis=1)
        row_sums[row_sums == 0] = 1 # Replace 0 sums with 1 to avoid division by zero, these rows will result in 0 proportions.
        grouped_df_scaled = grouped_df.div(row_sums, axis=0).fillna(0) # Scale and fill NaNs with 0

        # If error bars are enabled, normalize them as well.
        if add_error_bars and grouped_std is not None:
            grouped_std_scaled = grouped_std.div(row_sums, axis=0).fillna(0)
    else:
        grouped_df_scaled = grouped_df.copy()
        if add_error_bars and grouped_std is not None:
            grouped_std_scaled = grouped_std.copy()


    # ~ Determine Plotting Order based on value_order parameter ~ #
    sorted_categories = []
    if value_order == 'default':
        # Use the order as provided in cols_to_plot_y, but ensure 'others' is present
        sorted_categories = [col for col in cols_to_plot_y if col != 'others']
        if 'others' in cols_to_plot_y:
            sorted_categories.append('others')
    elif value_order == 'median_descending':
        # Calculate the median contribution for each category across all groups
        median_contributions = grouped_df_scaled.median()
        # Sort categories based on their median contribution in descending order (for biggest at bottom)
        sorted_categories = median_contributions.sort_values(ascending=False).index.tolist()
        # Ensure 'others' is always at the top position
        if 'others' in sorted_categories:
            sorted_categories.remove('others')
            sorted_categories.append('others')
    elif value_order == 'median_ascending':
        # Calculate the median contribution for each category across all groups
        median_contributions = grouped_df_scaled.median()
        # Sort categories based on their median contribution in ascending order (for smallest at bottom )
        sorted_categories = median_contributions.sort_values(ascending=True).index.tolist()
        # Ensure 'others' is always at the top position
        if 'others' in sorted_categories:
            sorted_categories.remove('others')
            sorted_categories.append('others')

    # Reorder the columns of the scaled DataFrame for plotting
    grouped_df_scaled = grouped_df_scaled[sorted_categories]

    # Reorder the std DataFrame as well if it exists
    if add_error_bars and grouped_std is not None:
        grouped_std_scaled = grouped_std_scaled[sorted_categories]

    # Reorder the colors list to match the sorted categories
    final_colors_for_plotting = [category_colors_internal_map[cat] for cat in sorted_categories]


    # ~ Determine Meta-column Order ~ #
    final_meta_order = []
    if meta_order is not None:
        # Filter meta_order to include only values present in grouped_df_scaled.index
        filtered_meta_order = [m for m in meta_order if m in grouped_df_scaled.index]
        grouped_df_scaled = grouped_df_scaled.reindex(filtered_meta_order)
        if add_error_bars and grouped_std is not None:
            grouped_std_scaled = grouped_std_scaled.loc[filtered_meta_order, sorted_categories]
        final_meta_order = filtered_meta_order
    else: # Default order for meta_column if no meta_order
        grouped_df_scaled = grouped_df_scaled
        if add_error_bars and grouped_std is not None:
            grouped_std_scaled = grouped_std_scaled
        final_meta_order = grouped_df_scaled.index.tolist()

    # ~ Plotting ~ #
    if ax_width:
        # Use axes of absolute size in a large aenough figure canvas
        ax_height = ax_width / aspect
        fig_width = 2*ax_width
        fig_height = 2*ax_height
        fig = plt.figure(figsize=(fig_width, fig_height))
        ax = fig.add_axes([0.25, 0.25, ax_width / fig_width, ax_height / fig_height])
    else:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Calculate x-positions for bars and handle grouping
    x_positions = []
    x_labels = []
    current_x = 0
    group_label_data = [] # Stores (group_name, start_x, end_x) for brackets and labels
    group_member_counts = {} # To count members per group
    
    if group_by_column:
        internal_grouping_map = dataset.set_index(meta_column)[group_by_column].to_dict()
        previous_group = None
        group_start_x = 0

        # Check for consistency if group_by_column is used with meta_order
        if meta_order is not None:
            group_sequence = [internal_grouping_map.get(m) for m in final_meta_order]
            # Check if groups are contiguous
            if group_sequence and any(group_sequence[i] != group_sequence[i-1] and group_sequence[i-1] in group_sequence[i+1:]
                                      for i in range(1, len(group_sequence))):
                warnings.warn(f"The provided 'meta_order' results in non-contiguous groups for '{group_by_column}'. "
                              "Grouping brackets might appear broken or inconsistent. "
                              "Consider reordering 'meta_order' to keep groups together.", UserWarning)

        for i, meta_val in enumerate(final_meta_order):
            current_group = internal_grouping_map.get(meta_val)
            group_member_counts[current_group] = group_member_counts.get(current_group, 0) + 1

            if previous_group is not None and current_group != previous_group:
                # Store end_x for the previous group 
                group_label_data[-1] = (group_label_data[-1][0], group_label_data[-1][1], x_positions[-1] + group_bracket_horizontal_line_length)
                # Add spacing between groups
                current_x += group_spacing
                group_start_x = current_x # Update the start for the new group
                # Start new group data
                group_label_data.append((current_group, group_start_x - group_bracket_horizontal_line_length, None)) # Adjust for bar width
            elif not group_label_data: # First item, start the first group
                group_label_data.append((current_group, current_x - group_bracket_horizontal_line_length, None)) # Adjust for bar width
            previous_group = current_group

            x_positions.append(current_x)
            x_labels.append(meta_val)
            current_x += 1 # Standard bar width + default matplotlib spacing

        if group_label_data:
            # Update end_x for the last group 
            group_label_data[-1] = (group_label_data[-1][0], group_label_data[-1][1], x_positions[-1] + group_bracket_horizontal_line_length)

    else: # No grouping
        x_positions = np.arange(len(grouped_df_scaled.index))
        x_labels = grouped_df_scaled.index.tolist()


    # ~ Add connecting shades if requested ~ #
    if add_connecting_shades:
        cumulative_df_scaled = grouped_df_scaled.cumsum(axis=1)

        for category in sorted_categories:
            y_cumulative = cumulative_df_scaled[category].values

            if category == sorted_categories[0]:
                y_previous_shades = np.zeros(len(x_positions))
            else:
                prev_category_index = sorted_categories.index(category) - 1
                prev_category_name = sorted_categories[prev_category_index]
                y_previous_shades = cumulative_df_scaled[prev_category_name].values

            ax.fill_between(x_positions, y_previous_shades, y_cumulative,
                            color=category_colors_internal_map[category], alpha=connecting_shades_alpha,
                            edgecolor=None, linewidth=0)


    # Determine border properties
    border_kwargs = {}
    if add_category_border:
        border_kwargs['edgecolor'] = 'black'
        border_kwargs['linewidth'] = category_border_width

    # Plot the stacked bars manually to control x-positions and legend
    plotted_categories_for_legend = set() # To ensure each category appears only once in legend

    for i, meta_val in enumerate(grouped_df_scaled.index):
        bottom_val = 0
        current_x_pos = x_positions[i]
        for category in sorted_categories:
            height = grouped_df_scaled.loc[meta_val, category]
            yerr_val = grouped_std_scaled.loc[meta_val, category] if add_error_bars and grouped_std_scaled is not None else 0

            # Add label only for the first bar of this category encountered across all bars
            label = category if category not in plotted_categories_for_legend else "_nolegend_"

            ax.bar(current_x_pos, height=height, bottom=bottom_val,
                   color=category_colors_internal_map[category],
                   yerr=yerr_val if add_error_bars else None,
                   capsize=4,
                   width=bar_width,
                   label=label, # Assign label for legend
                   **border_kwargs)

            if category not in plotted_categories_for_legend:
                plotted_categories_for_legend.add(category)

            bottom_val += height

    ax.set_xticks(x_positions)
    ax.set_xlim(x_positions[0] - (0.5 + (0.5 - (bar_width/2))), x_positions[-1] + (0.5 + (0.5 - (bar_width/2))))

    # ~ ticks label adjustment based on rotation ~ #
    ha_for_set_xticklabels = 'center' # Default
    rotation_mode_for_xticklabels = None # Default

    # Automatic alignment based on x_tick_label_rotation
    if x_tick_label_rotation == 0:
        ha_for_set_xticklabels = 'center'
    elif x_tick_label_rotation > 0 and x_tick_label_rotation < 90:
        ha_for_set_xticklabels = 'right'
        rotation_mode_for_xticklabels = 'anchor'
    elif x_tick_label_rotation < 0 and x_tick_label_rotation > -90:
        ha_for_set_xticklabels = 'left'
        rotation_mode_for_xticklabels = 'anchor'
    else: # For 90 degrees or other angles
        ha_for_set_xticklabels = 'right'
        rotation_mode_for_xticklabels = 'anchor'

    # Set the tick labels using the determined string alignment and rotation mode
    ax.set_xticklabels(x_labels, 
                       rotation=x_tick_label_rotation, 
                       ha=ha_for_set_xticklabels,
                       va="top",
                       fontsize=x_tick_label_fontsize,
                       rotation_mode=rotation_mode_for_xticklabels)
    
    # Set robust y pad
    y_transform = ax.get_xaxis_transform()
    # Tranform points to proper axis coordinates
    y_display_at_zero = y_transform.transform((0, 0))[1] 
    y_display_at_padding = y_transform.transform((0, x_ticks_label_pad))[1]
    # Gget the figure's dpi to make the conversion more accurate
    fig_dpi = fig.dpi
    # Convert pixels to points: pixels * (72 / dpi)
    padding_in_pixels = abs(y_display_at_padding - y_display_at_zero)
    xticks_label_pad_calculated = padding_in_pixels * (72 / fig_dpi)
    ax.tick_params(axis='x', pad=xticks_label_pad_calculated)

    # Add group brackets and labels if group_by_column is provided
    if group_by_column:
        if group_position == 'bottom':
            bracket_horizontal_y_level = -0.05 + group_label_y_offset # Base offset + user-defined offset
            label_y_level = -0.15 + group_label_y_offset  # Base offset + user-defined offset
            bracket_line_offset = -1 * group_bracket_vertical_line_length 
        elif group_position == 'middle':
            # In 'middle' position, group_label_y_offset controls the horizontal line's position
            bracket_horizontal_y_level = group_label_y_offset 
            # The label's position is then relative to this horizontal line
            label_y_level = bracket_horizontal_y_level - group_bracket_vertical_line_length - 0.05
            
            bracket_line_offset = group_bracket_vertical_line_length

        elif group_position == 'top':
            # Position above the plot, near the title.
            bracket_horizontal_y_level = 1.05 + group_label_y_offset
            label_y_level = 1.15 + group_label_y_offset
            bracket_line_offset = -1 * group_bracket_vertical_line_length 

        for group_name, start_x_coord, end_x_coord in group_label_data:

            if start_x_coord is None or end_x_coord is None: # Skip incomplete group data
                continue

            # Only draw brackets if the group has more than one member
            if group_member_counts.get(group_name, 0) > 1:
                
                # Effective x-coordinates adjusted for bar width
                effective_start_x = start_x_coord
                effective_end_x = end_x_coord

                # Create bracket
                # Draw left vertical line of bracket
                ax.plot([effective_start_x, effective_start_x], [bracket_horizontal_y_level, bracket_horizontal_y_level + bracket_line_offset],
                     color='black', transform=ax.get_xaxis_transform(), clip_on=False, linewidth=group_bracket_linewidth)
                # Draw horizontal line of bracket
                ax.plot([effective_start_x, effective_end_x], [bracket_horizontal_y_level + bracket_line_offset, bracket_horizontal_y_level + bracket_line_offset],
                        color='black', transform=ax.get_xaxis_transform(), clip_on=False, linewidth=group_bracket_linewidth)
                # Draw right vertical line of bracket
                ax.plot([effective_end_x, effective_end_x], [bracket_horizontal_y_level, bracket_horizontal_y_level + bracket_line_offset],
                        color='black', transform=ax.get_xaxis_transform(), clip_on=False, linewidth=group_bracket_linewidth)

                # Add group label
                if show_group_label:
                    ax.text((effective_start_x + effective_end_x) / 2, label_y_level,
                            group_name, ha='center', va='top',
                            fontsize=group_label_fontsize, rotation=group_label_rotation,
                            transform=ax.get_xaxis_transform(), clip_on=False)
                    
    # ~ Add boxes  ~ #

    if boxes_column:
        # Determine the order for unique_top_box_values
        if boxes_color_map:
            # Use keys from boxes_color_map for order, then add any remaining unique values
            ordered_unique_top_box_values = list(boxes_color_map.keys())
            # Add any unique values from the column that aren't in the color map
            for val in dataset[boxes_column].unique():
                if val not in ordered_unique_top_box_values:
                    ordered_unique_top_box_values.append(val)
        else:
            # If no color map is provided, use the natural order of unique values from the dataset
            ordered_unique_top_box_values = dataset[boxes_column].unique().tolist()


        # Populate missing colors in boxes_color_map
        current_top_box_color_idx = 0
        final_boxes_color_map = {}
        for val in ordered_unique_top_box_values:
            if val in boxes_color_map:
                final_boxes_color_map[val] = boxes_color_map[val]
            else:
                final_boxes_color_map[val] = default_colors[current_top_box_color_idx % len(default_colors)]
                current_top_box_color_idx += 1
            
        # Create a mapping from meta_column values
        meta_to_top_box_value = dataset.set_index(meta_column)[boxes_column].to_dict()

        for i, meta_val in enumerate(final_meta_order):
            box_value = meta_to_top_box_value.get(meta_val)
            if box_value is not None:
                box_color = final_boxes_color_map.get(box_value, 'grey')
                # Calculate box position relative to the bar X-postion
                box_x_left = x_positions[i] - (boxes_width / 2.0)

                # Add the rectangle patch
                rect = mpatches.Rectangle((box_x_left, boxes_y_position),
                                          boxes_width,
                                          boxes_height,
                                          facecolor=box_color,
                                          edgecolor='black',
                                          linewidth=boxes_borderwidth,
                                          transform=ax.get_xaxis_transform(),
                                          clip_on=False)
                
                ax.add_patch(rect)

        if boxes_legend:
            
            if boxes_legend_title is not None:
                legend_title_str = boxes_legend_title
            else:
                legend_title_str = boxes_column

            # Add a separate legend for the top boxes
            top_box_legend_handles = []
            for val in ordered_unique_top_box_values: # Iterate through ordered values for sorted legend
                color = final_boxes_color_map.get(val, 'grey')
                top_box_legend_handles.append(mpatches.Patch(color=color, label=str(val)))

            # Create the second legend (for top boxes)
            if boxes_legend_pos == "bottom":

                # Initialize lists for all legend elements
                all_legend_handles = []
                all_legend_labels = []

                all_legend_handles.append(mlines.Line2D([], [], color='none', marker='None', linestyle='None'))
                replacements = str.maketrans({"_": "\_", " ": "\ "})
                all_legend_labels.append(r"$\mathbf{" + legend_title_str.translate(replacements) + r"}$")

                # Add the colored square patches and their labels for the actual categories
                for val in ordered_unique_top_box_values:
                    color = final_boxes_color_map.get(val, 'grey')
                    all_legend_handles.append(mpatches.Patch(color=color)) # Patch for the square color
                    all_legend_labels.append(str(val)) # Label for the square

                top_box_legend = ax.legend(
                    handles=all_legend_handles,
                    labels=all_legend_labels,
                    bbox_to_anchor=(boxes_legend_x_pos, boxes_legend_y_pos),
                    loc='lower center', # Position at the bottom center
                    ncol=len(all_legend_labels), # Ensure all elements are on a single row
                    fontsize=boxes_legend_fontsize,
                    title_fontsize=0,
                    frameon=False, # No frame around the legend
                    handlelength=0.7, # Default handle length for patches
                    handletextpad=0.5, # Space between handle (square) and its text
                    columnspacing=0.5 # Compact spacing between legend columns
                )
            else:
            
                top_box_legend = ax.legend(handles=top_box_legend_handles,
                                        title=legend_title_str,
                                        bbox_to_anchor=(boxes_legend_x_pos, boxes_legend_y_pos),
                                        loc="center left",
                                        fontsize=boxes_legend_fontsize,
                                        title_fontsize=boxes_legend_fontsize)
        
            # Manually add the first legend back to the figure
            ax.add_artist(top_box_legend)


    # Set y-axis limits to ensure white space at the top
    max_y_value = grouped_df_scaled.sum(axis=1).max()
    if add_error_bars and grouped_std_scaled is not None:
        # If error bars are present, include their upper bound in the max_y_value calculation
        max_y_value += grouped_std_scaled.max().max() # Add the largest standard deviation
    
    # Add a pad for Y-axis limit
    y_upper_limit = max_y_value * (1+y_upper_pad)
    ax.set_ylim(0, y_upper_limit)

    # Basic label and font-size setting
    title_suffix = ""
    if scaling == 'none':
        title_suffix += "Counts"
    elif scaling == 'median':
        title_suffix += "Medians"
    elif scaling == 'mean':
        title_suffix += "Means"

    if normalize_data:
        title_suffix = "Proportions"

    if show_title:
        # Set titles and labels
        #tranform points to proper axis coordinates
        y_display_at_zero = y_transform.transform((0, 0))[1] 
        y_display_at_padding = y_transform.transform((0, title_pad))[1]
        # Get the figure's dpi to make the conversion more accurate
        fig_dpi = fig.dpi
        # Convert pixels to points: pixels * (72 / dpi)
        padding_in_pixels = abs(y_display_at_padding - y_display_at_zero)
        title_label_pad_calculated = padding_in_pixels * (72 / fig_dpi)
        ax.tick_params(axis='x', pad=xticks_label_pad_calculated)

        if title is not None:
            plt.title(title, fontsize=title_fontsize,pad=title_label_pad_calculated)
        else:
            plt.title(f'Stacked Bar Plot of {title_suffix} by {meta_column}', fontsize=title_fontsize, pad=title_label_pad_calculated)

    ax.set_axisbelow(True)
    plt.grid(axis='y', linestyle='-', alpha=0.7)

    #Hide top spine
    if hide_top_spine:
        ax.spines['top'].set_visible(False)

    #Hide right spine
    if hide_right_spine:
        ax.spines['right'].set_visible(False)

    #Hide bottom spine
    if hide_bottom_spine:
        ax.spines['bottom'].set_visible(False)
        ax.axhline(y=0.0, color='grey', linestyle='-', alpha=0.7, zorder=1.5)

    #Hide left spine
    if hide_left_spine:
        ax.spines['left'].set_visible(False)

    #Hide bottom ticks
    if hide_bottom_tick:
        plt.tick_params(bottom=False)

    #Hide left ticks
    if hide_left_tick:
        plt.tick_params(left=False)

    # Set x-axis label
    if show_xlabel:
        if xlabel is not None:
            plt.xlabel(xlabel, fontsize=xlabel_fontsize)
        else:
            plt.xlabel(meta_column, fontsize=xlabel_fontsize) # Use xlabel_fontsize

    # Set y-axis label
    if show_ylabel:
        if ylabel is not None:
            plt.ylabel(ylabel, fontsize=ylabel_fontsize)
        else:
            ylabel_text = "Proportion" if normalize_data else ( "Count" if scaling == 'none' else ( "Median Value" if scaling == 'median' else "Mean Value" ))
            plt.ylabel(ylabel_text, fontsize=ylabel_fontsize)

    plt.yticks(fontsize=y_tick_label_fontsize)
    plt.xticks(fontsize=x_tick_label_fontsize)

    # Create the first legend (for stacked bars)
    main_legend_title = legend_title if legend_title is not None else 'Category'

    if legend_position != "custom bottom":
        main_legend = ax.legend(title=main_legend_title, bbox_to_anchor=(legend_x_pos, legend_y_pos), loc="center left", 
                                fontsize=legend_fontsize, title_fontsize=legend_title_fontsize)
    else:
        ordered_unique_main_legend_values = list(category_colors_internal_map.keys())
        all_legend_handles = []
        all_legend_labels = []

        all_legend_handles.append(mlines.Line2D([], [], color='none', marker='None', linestyle='None'))
        replacements = str.maketrans({"_": "\_", " ": "\ "})
        all_legend_labels.append(r"$\mathbf{" + main_legend_title.translate(replacements) + r"}$")

        # Add the colored square patches and their labels for the actual categories
        for val in ordered_unique_main_legend_values:
            color = category_colors_internal_map.get(val, 'grey')
            all_legend_handles.append(mpatches.Patch(color=color)) # Patch for the square color
            all_legend_labels.append(str(val)) # Label for the square

            main_legend = ax.legend(
                handles=all_legend_handles,
                labels=all_legend_labels,
                bbox_to_anchor=(0.5, legend_y_pos),
                loc='lower center', # Position at the bottom center
                ncol=len(all_legend_labels), # Ensure all elements are on a single row
                fontsize=legend_fontsize,
                title_fontsize=0,
                frameon=False, # No frame around the legend
                handlelength=0.7, # Default handle length for patches
                handletextpad=0.5, # Space between handle (square) and its text
                columnspacing=0.5 # Compact spacing between legend columns
                )
            
    ax.add_artist(main_legend)
    ax.get_legend().remove()
    ax.grid(axis='x', visible=False)
    if yticks:
        ax.set_yticks(yticks)

    # ~ Saving Plot (if output path is provided) ~ #
    if output:
        dir_name = os.path.dirname(output)
        if dir_name and not os.path.exists(dir_name):
            os.makedirs(dir_name)

    # select legend to plot
    if boxes_legend_pos == "bottom":
        artist_elements = [main_legend,top_box_legend]
    else:
        artist_elements = [main_legend]

    filename_pdf = output + ".pdf"
    plt.savefig(filename_pdf, format='pdf', dpi=dpi, bbox_inches='tight', bbox_extra_artists=artist_elements)
    print(f"Box plot saved to {filename_pdf}")

    filename_png = output + ".png"
    plt.savefig(filename_png, format='png', dpi=dpi, bbox_inches='tight',bbox_extra_artists=artist_elements)
    print(f"Box plot saved to {filename_png}")

    filename_svg = output + ".svg"
    plt.rcParams["svg.fonttype"] = "none"
    plt.savefig(filename_svg, format='svg', dpi=dpi, bbox_inches='tight',bbox_extra_artists=artist_elements)
    print(f"Box plot saved to {filename_svg}")

    return dataframe_before_scaling, dataframe_after_scaling, ax, artist_elements