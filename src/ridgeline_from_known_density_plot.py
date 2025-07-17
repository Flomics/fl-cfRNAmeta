import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
import seaborn as sns

def ridgeline_from_known_density_plot(# ~ Basic input ~ #
                                  dataset: pd.DataFrame,
                                  x_column: str,
                                  density_column: str,
                                  category_column: str,
                                  category_order: list = None,
                                  overlap: float = 0.5,
                                  fig_width: float = 10,
                                  fig_height: float = 7,
                                  output: str = None,
                                  # ~ Layout ~ #
                                  show_xlabel: bool = True,
                                  xlabel: str = None,
                                  xlabel_va: float = 0.5,
                                  xlabel_fontsize: float = 12,
                                  show_ylabel: bool = True,
                                  ylabel: str = None,
                                  ylabel_fontsize: float = 12,
                                  show_title: bool = True,
                                  title: str = None,
                                  title_fontsize: float = 14,
                                  xticks_fontsize: float = 10,
                                  show_yticks: bool = True,
                                  yticks_fontsize: float = 10,
                                  y_ticks_frequency: float = 1,
                                  fill_alpha: float = 0.7,
                                  line_color: str = 'black',
                                  fill_color: str = 'skyblue',
                                  line_width: float = 0.5,
                                  # ~ Normalization ~ #
                                  normalization: str = None,
                                  # ~ Individual sub-plots ~ #
                                  show_individual_yaxis: bool = False,
                                  show_individual_hspace: float = 0.4,
                                  consistent_y_scale: bool = True): # New parameter
    """
    Generate a ridgeline plot from pre-calculated density values.

    Args:

        # ~ Basic input ~ #
        dataset (pd.DataFrame): The input pandas DataFrame containing pre-calculated densities values
        x_column (str): Column containing the coordinates for the density curves.
        density_column (str): Column containing the density values for the curves.
        category_column (str): Column containing categories to group the curves.
        category_order (list, optional): Order in which categories should be plotted from bottom to top.
        overlap (float, optional): The degree of overlap between the lines (from 0 to 1). Defaults to 0.5.
        fig_width (float, optional): The width of the figure in inches. Defaults to 10.
        fig_height (float, optional): The height of the figure in inches. Defaults to 7.
        output (str, optional): The full path filename (without extension) to save the plot.

        # ~ Layout ~ #
        show_xlabel (bool, optional): Display the X-axis label. Defaults to True.
        xlabel (str, optional): Custom label for the X-axis. Defaults to None.
        xlabel_va (float, optional): X-axis label vertical adjustment. Default to 0.5.
        xlabel_fontsize (float, optional): Font size for the X-axis label. Defaults to 12.
        show_ylabel (bool, optional): Display the Y-axis label. Defaults to True
        ylabel (str, optional): Custom label for the Y-axis. Defaults to None.
        ylabel_fontsize (float, optional): Font size for the Y-axis label. Defaults to 12.
        show_title (bool, optional): Display plot title. Defaults to True.
        title (str, optional): Custom title. Defaults to None.
        title_fontsize: float (optional): Font size for title. Defaults to 14.
        xticks_fontsize (float, optional): Font size for the X-axis tick labels. Defaults to 10.
        show_yticks (bool, optional): Show Y-axis ticks and label. Default to
        yticks_fontsize (float, optional): Font size for the Y-axis tick labels. Defaults to 10.
        y_ticks_frequency (float, optional): Frequency of Y-axis ticks. Default to 1
        fill_alpha (float, optional): Transparency of the filled area under the curves. Defaults to 0.7.
        line_color (str, optional): Color of the curve lines. Defaults to 'black'.
        fill_color (str or dict, optional): Color of the filled area under the curves.
                                            Can be a single color string for all categories,
                                            or a dictionary mapping category names
                                            Defaults to 'skyblue'.
        line_width (float, otiopnal): Size of the line in the curve. Default to 0.5.

        # ~ Normalization ~ #

        normalization (str, optional): normalization method for the curve. Default to None.
                                        If "mean", {category_column} and {x_column} will be grouped to calculate the {density_column} mean
                                        If "median", {category_column} and {x_column} will be grouped to calculate the {density_column} median

        # ~ Individual sub-plots ~ #

        show_individual_yaxis (bool, optional): If True, each category will have its own y-axis.
                                                This overrides 'overlap' to 0.
                                                Defaults to False.
        show_individual_hspace (float, optional): The height spacing between subplots when `show_individual_yaxis` is True. Defaults to 0.4.
        consistent_y_scale (bool, optional): If False, The y scale will be different for each subplot. Default to True
    """

    # Input Validation
    if not isinstance(dataset, pd.DataFrame):
        raise TypeError("Input dataset must be a pandas DataFrame.")
    if x_column not in dataset.columns:
        raise ValueError(f"Column '{x_column}' not found in the dataset.")
    if density_column not in dataset.columns:
        raise ValueError(f"Column '{density_column}' not found in the dataset.")
    if category_column not in dataset.columns:
        raise ValueError(f"Column '{category_column}' not found in the dataset.")
    if not pd.api.types.is_numeric_dtype(dataset[x_column]):
        raise TypeError(f"Column '{x_column}' must be numeric.")
    if not pd.api.types.is_numeric_dtype(dataset[density_column]):
        raise TypeError(f"Column '{density_column}' must be numeric.")
    if not (0 <= overlap <= 1):
        raise ValueError("Overlap must be between 0 and 1.")
    if category_order is not None and not isinstance(category_order, list):
        raise TypeError("category_order must be a list.")
    if category_order is not None:
        if not all(cat in dataset[category_column].unique() for cat in category_order):
            raise ValueError("All categories in category_order must be present in category_column.")
    if not isinstance(fill_color, (str, dict)):
        raise TypeError("fill_color must be a string or a dictionary.")
    if normalization is not None and normalization not in ["mean","median"]:
        raise ValueError("Invalid value for `normalization`. Choose from 'mean' or 'median'.")
    if show_individual_yaxis:
        if overlap != 0:
            warnings.warn("`show_individual_yaxis` is True, so `overlap` is being set to 0 to prevent overlap.", UserWarning)
        overlap = 0

    # ~ Normalization ~ #
    if normalization == "mean":
        dataset = dataset.groupby([category_column, x_column])[density_column].mean().reset_index()
    if normalization == "median":
        dataset = dataset.groupby([category_column, x_column])[density_column].median().reset_index()

    # ~ Get unique categories and define order ~ #
    categories = dataset[category_column].unique()
    if category_order is None:
        categories = sorted(categories)
    else:
        # Filter and order categories based on category_order, ensuring they exist in the data
        categories = [cat for cat in category_order if cat in categories]

    # ~ Assigning colors based on fill_color type ~ #
    category_fill_colors = {}
    if isinstance(fill_color, str):
        # If a single color string, apply to all categories
        for cat in categories:
            category_fill_colors[cat] = fill_color
    else: # It's a dictionary
        prop_cycle = plt.rcParams['axes.prop_cycle']
        default_colors = prop_cycle.by_key()['color']
        color_idx = 0
        for cat in categories:
            if cat in fill_color:
                category_fill_colors[cat] = fill_color[cat]
            else:
                # Use default cycle color if not specified in the dictionary
                category_fill_colors[cat] = default_colors[color_idx % len(default_colors)]
                warnings.warn(f"No fill_color specified for category '{cat}'. Using default color: {category_fill_colors[cat]}", UserWarning)
                color_idx += 1

    if show_individual_yaxis:
        # Create subplots
        fig, axes = plt.subplots(len(categories), 1, figsize=(fig_width, fig_height), sharex=True, gridspec_kw={'hspace': show_individual_hspace})
        if len(categories) == 1: # If only one subplot, axes is not an array
            axes = [axes]
    else:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # ~ Calculate the maximum density across all plots to determine the appropriate vertical shift ~ #
    if not show_individual_yaxis:
        max_overall_density = dataset[density_column].max()
        if pd.isna(max_overall_density) or max_overall_density == 0:
            vertical_shift_amount = 0.5 # A default small shift if no density is found
            warnings.warn("No density found in the density_column. Plots may appear flat or empty.", UserWarning)
        else:
            vertical_shift_amount = max_overall_density * (1 - overlap)
    else:
        if consistent_y_scale:
            vertical_shift_amount = 0 # No vertical shift for individual y-axes
            # Determine max density for consistent Y-axis limit across subplots
            global_max_density = dataset[density_column].max()
            if pd.isna(global_max_density) or global_max_density == 0:
                global_max_density = 1.0
                

    # ~ Plot each density curve ~ #
    for i, category in enumerate(categories):
        subset = dataset[dataset[category_column] == category].sort_values(by=x_column)
        x_kde = subset[x_column].values
        y_kde = subset[density_column].values

        if len(y_kde) == 0:
            warnings.warn(f"Category '{category}' has no valid density data. Skipping plot for this category.", UserWarning)
            continue

        if show_individual_yaxis:
            current_ax = axes[i]
            y_shifted = y_kde # No shift needed, each has its own axis
        else:
            current_ax = ax
            y_shifted = y_kde + (i * vertical_shift_amount)

        # Get the appropriate fill color for the current category
        current_fill_color = category_fill_colors.get(category, 'gray') # Fallback to gray if not mapped

        # Plot the filled area
        if show_individual_yaxis:
            current_ax.fill_between(x_kde, 0, y_shifted, # Start fill from 0 for individual y-axis
                                    color=current_fill_color, alpha=fill_alpha)
            if consistent_y_scale:
                current_ax.yaxis.set_ticks(np.arange(0, global_max_density * 1.1, y_ticks_frequency)) # Set consistent Y-axis limits
            else:
                max_density = dataset[dataset[category_column] == category][density_column].max()
                current_ax.yaxis.set_ticks(np.arange(0, max_density * 1.1, y_ticks_frequency)) # Set variable Y-axis limit
        else:
            current_ax.fill_between(x_kde, i * vertical_shift_amount, y_shifted,
                                    color=current_fill_color, alpha=fill_alpha)


        # Plot the outline
        current_ax.plot(x_kde, y_shifted, color=line_color, linewidth=line_width)

        if show_individual_yaxis:
            current_ax.tick_params(axis='x', labelsize=xticks_fontsize)
            if show_yticks:
                current_ax.tick_params(axis='y', labelsize=yticks_fontsize)
            else:
                current_ax.set_yticks([])
            current_ax.spines['right'].set_visible(False)
            current_ax.spines['top'].set_visible(False)
            current_ax.grid(axis='x', linestyle='--', alpha=0.7)
            if show_ylabel:
                current_ax.set_ylabel(category, fontsize=ylabel_fontsize, rotation=0, ha='right', va='center')
                current_ax.yaxis.set_label_coords(-0.05, 0.5)
            else:
                current_ax.set_ylabel('')

        else:
            # Add category label as a custom y-tick for ridgeline plot
            current_ax.text(-0.02, i * vertical_shift_amount, str(category), va='center', ha='right',
                            fontsize=yticks_fontsize, color='black', transform=current_ax.get_yaxis_transform())


    # ~ Set x-axis label ~ #
    if show_xlabel:
        if show_individual_yaxis:
            # Place x-label on the bottom-most subplot
            bottom_ax = axes[-1]
            bottom_ax.set_xlabel(xlabel if xlabel is not None else x_column, fontsize=xlabel_fontsize, x = xlabel_va, ha='center')
            # Ensure ticks are visible on the bottom-most axis
            bottom_ax.tick_params(axis='x', labelbottom=True)
        else:
            ax.set_xlabel(xlabel if xlabel is not None else x_column, fontsize=xlabel_fontsize, x = xlabel_va, ha='center')
    else:
        if show_individual_yaxis:
            axes[-1].fig.supxlabel('')
        else:
            ax.set_xlabel('')

    # ~ Set y-axis label and ticks (hide default numerical y-axis for ridgeline) ~ #
    if not show_individual_yaxis:
        if show_ylabel:
            ax.set_ylabel(ylabel if ylabel is not None else 'Category', fontsize=ylabel_fontsize)
        else:
            ax.set_ylabel('')

        ax.set_yticks([])
        ax.tick_params(axis='y', length=0)
        ax.tick_params(axis='x', labelsize=xticks_fontsize)
        # ~ Remove top, right, and left plot lines for a cleaner look ~ #
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.grid(axis='x', linestyle='--', alpha=0.7)


    # ~ Set title ~ #
    if show_title:
        if show_individual_yaxis:
            fig.suptitle(title if title is not None else f'Density Plots of {x_column} by {category_column}',
                         fontsize=title_fontsize)
        else:
            ax.set_title(title if title is not None else f'Ridgeline Plot of {x_column} Densities by {category_column}',
                         fontsize=title_fontsize)

    #normal layout for regular plot
    if not show_individual_yaxis:
        plt.tight_layout()
    # Use fig.subplots_adjust for individual curve plots
    #TO DO: make it automatic
    if show_individual_yaxis:
        fig.subplots_adjust(left = 0.3, bottom=0.03, top=1)

    # ~ Saving Plot (if output path is provided) ~
    if output:
        dir_name = os.path.dirname(output)
        if dir_name and not os.path.exists(dir_name):
            os.makedirs(dir_name)

        filename_pdf = output + ".pdf"
        plt.savefig(filename_pdf)
        print(f"Plot saved to {filename_pdf}")

        filename_png = output + ".png"
        plt.savefig(filename_png)
        print(f"Plot saved to {filename_png}")

        filename_svg = output + ".svg"
        plt.rcParams["svg.fonttype"] = "none"
        plt.savefig(filename_svg)
        print(f"Box plot saved to {filename_svg}")

    plt.show()