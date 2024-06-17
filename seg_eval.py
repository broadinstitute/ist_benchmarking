"""
Author: 
Huan Wang (whuan@broadinstitute.org)
STP, Broad Institute of MIT and Harvard

Functions:
    calculate_area_bounds()
    calculate_seg_eval_metrics()
    dapi_mask_annotation()
    calculate_tp_fn_fp_from_gdf()
    dapi_mem_mask_annotation()
    gdf_flip()
    load_cell_boundary_parquet()
    plot_cell_filtration()
    plot_cell_transcripts()
    translate_to_bbox()
"""

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
from shapely.affinity import scale
import numpy as np
import matplotlib.pyplot as plt
import glob


def calculate_area_bounds(percentage, total_size):
    """
    Calculate the bounding box coordinates for a given percentage area of a square plot.

    Parameters:
        percentage (float): Desired percentage of the area (e.g., 25, 50, 75).
        total_size (int): Total size of one side of the square (default 4256).

    Returns:
        dict: A dictionary with 'xmin', 'xmax', 'ymin', 'ymax'.
    """
    # Calculate the side length of the square that corresponds to the desired area percentage
    area_side_length = np.sqrt((percentage / 100) * (total_size ** 2))

    # Calculate the center of the original square
    center = total_size / 2

    # Calculate min and max coordinates based on the area_side_length
    xmin = center - (area_side_length / 2)
    xmax = center + (area_side_length / 2)
    ymin = xmin  # Symmetry in square
    ymax = xmax

    return {
        'xmin': int(xmin),
        'xmax': int(xmax),
        'ymin': int(ymin),
        'ymax': int(ymax)
    }


def calculate_seg_eval_metrics(tp, fn, fp):
    """
    Calculate precision, recall, F1 score, and IoU based on true positives, false negatives, and false positives.

    Parameters:
    tp (int): Number of true positive cells
    fn (int): Number of false negative cells
    fp (int): Number of false positive cells

    Returns:
    dict: A dictionary containing precision, recall, F1 score, and IoU
    """
    # Precision (Positive Predictive Value)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    
    # Recall (Sensitivity, True Positive Rate)
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    
    # F1 Score
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    # Intersection over Union (IoU), also known as the Jaccard Index
    iou = tp / (tp + fp + fn) if (tp + fp + fn) > 0 else 0
    
    return {
        'precision': precision,
        'recall': recall,
        'f1_score': f1_score,
        'iou': iou
    }


def calculate_tp_fn_fp_from_gdf(gdf_mask, gdf_gt):
    """
    Calculate true positives, false negatives, and false positives based on overlap
    of mask polygons with ground truth centroids.

    Parameters:
        gdf_mask (GeoDataFrame): A GeoDataFrame containing polygons of cell segmentations.
        gdf_gt (GeoDataFrame): A GeoDataFrame containing points of ground truth centroids.

    Returns:
        tuple: (true positives, false negatives, false positives)
    """
    # Spatial join to find overlaps (true positives)
    tp_gdf = gpd.sjoin(gdf_gt, gdf_mask, how='left', op='within')

    # True Positives (TP): Ground truth points that overlap with a mask polygon
    tp = tp_gdf['index_right'].notna().sum()

    # False Negatives (FN): Ground truth points that do not overlap with any mask polygon
    fn = tp_gdf['index_right'].isna().sum()

    # False Positives (FP): Mask polygons that do not contain any ground truth point
    # Count mask polygons then subtract count of unique mask polygons that overlap with ground truths
    all_masks = len(gdf_mask)
    matched_masks = tp_gdf['index_right'].dropna().unique().shape[0]
    fp = all_masks - matched_masks

    return tp, fn, fp


def dapi_mask_annotation(img_dapi, gdf_mask, gdf_truth, figwidth, 
                         q0, q1, scale_bar, xy_range=False, area_bounds={},
                         lw=1, markersize=1, dapi_cmap='bone', boundary_color='white', 
                         gt_color='red', save=False, fname=''):
    """
    Visualize DAPI stained image with overlaid segmentation and ground truth annotations.
    
    Parameters:
        img_dapi (array): 2D array of the DAPI channel.
        gdf_mask (GeoDataFrame): Geodataframe containing the segmentation mask.
        gdf_truth (GeoDataFrame): Geodataframe containing the ground truth annotations.
        figwidth (float): Width of the figure.
        q0 (float): Lower quantile for DAPI image normalization.
        q1 (float): Upper quantile for DAPI image normalization.
        scale_bar (matplotlib.artist.Artist): Scale bar to be added to the plot.
        xy_range (bool): Whether to apply custom limits to the x and y axes.
        area_bounds (dict): Dictionary with 'xmin', 'xmax', 'ymin', 'ymax' if xy_range is True.
        lw (float): Line width for plotting the segmentation mask.
        markersize (float): Size of the marker for ground truth points.
        dapi_cmap (str): Colormap for DAPI image.
        boundary_color (str): Color for the boundary of the segmentation mask.
        gt_color (str): Color for the ground truth annotations.
        save (bool): Whether to save the figure to a file.
        fname (str): Filename or path to save the figure if save is True.

    Returns:
        None
    """
    figheight = figwidth * img_dapi.shape[0] / img_dapi.shape[1]
    fig, ax = plt.subplots(1, 1, figsize=(figwidth, figheight))
    ax.imshow(
        img_dapi, vmin=np.quantile(img_dapi, q0), vmax=np.quantile(img_dapi, q1), cmap=dapi_cmap)
    gdf_mask.plot(aspect=1, ax=ax, legend=False,
                  edgecolor=boundary_color,
                  facecolor='none',
                  linewidth=lw)
    gdf_truth.plot(
        aspect=1, ax=ax, legend=False,
        edgecolor=gt_color,
        facecolor=gt_color,
        linewidth=markersize)
    ax.add_artist(scale_bar)
    if xy_range:
        ax.set_xlim(area_bounds['xmin'], area_bounds['xmax'])
        ax.set_ylim(area_bounds['ymin'], area_bounds['ymax'])
    plt.axis('off')
    plt.show()
    if save:
        fig.savefig(fname, format='png', dpi=200)
        fig.savefig(fname.replace('png', 'eps'), format='eps', dpi=200)


def dapi_mem_mask_annotation(img_dapi, img_mem, gdf_mask, gdf_truth, figwidth, 
                             q_dapi, q_mem, scale_bar, xy_range=False, area_bounds={},
                             lw=1, markersize=1, boundary_color='white', 
                             gt_color='red', save=False, fname=''):    
    """
    Visualize a composite image combining DAPI and membrane staining with
    overlaid segmentation and ground truth annotations.
    
    Parameters:
        img_dapi (array): 2D array of the DAPI channel.
        img_mem (array): 2D array of the membrane channel.
        gdf_mask (GeoDataFrame): Geodataframe containing the segmentation mask.
        gdf_truth (GeoDataFrame): Geodataframe containing the ground truth annotations.
        figwidth (float): Width of the figure.
        q_dapi (float): Upper quantile for DAPI image normalization.
        q_mem (float): Upper quantile for membrane image normalization.
        scale_bar (matplotlib.artist.Artist): Scale bar to be added to the plot.
        xy_range (bool): Whether to apply custom limits to the x and y axes.
        area_bounds (dict): Dictionary with 'xmin', 'xmax', 'ymin', 'ymax' if xy_range is True.
        lw (float): Line width for plotting the segmentation mask.
        markersize (float): Size of the marker for ground truth points.
        boundary_color (str): Color for the boundary of the segmentation mask.
        gt_color (str): Color for the ground truth annotations.
        save (bool): Whether to save the figure to a file.
        fname (str): Filename or path to save the figure if save is True.

    Returns:
        None
    """
    figheight = figwidth * img_dapi.shape[0] / img_dapi.shape[1]
    fig, ax = plt.subplots(1, 1, figsize=(figwidth, figheight))

    r = np.zeros_like(img_dapi)
    g = img_mem
    b = img_dapi

    g = g / np.percentile(g, q_mem * 100)
    b = b / np.percentile(b, q_dapi * 100)
    stack = np.dstack((r, g, b))
    
    ax.imshow(stack)
    gdf_mask.plot(aspect=1, ax=ax, legend=False,
                  edgecolor=boundary_color,
                  facecolor='none',
                  linewidth=lw)
    gdf_truth.plot(
        aspect=1, ax=ax, legend=False,
        edgecolor=gt_color,
        facecolor=gt_color,
        linewidth=markersize)
    ax.add_artist(scale_bar)
    if xy_range:
        ax.set_xlim(area_bounds['xmin'], area_bounds['xmax'])
        ax.set_ylim(area_bounds['ymin'], area_bounds['ymax'])
    plt.axis('off')
    plt.show()
    if save:
        fig.savefig(fname, format='png', dpi=200)
        fig.savefig(fname.replace('png', 'eps'), format='eps', dpi=200)


def gdf_flip(gdf, direction='lr'):
    """
    Flips the geometries in a GeoDataFrame horizontally or vertically.
    
    Parameters:
        gdf (GeoDataFrame): The GeoDataFrame to be flipped.
        direction (str): Direction of the flip, 'lr' for left-right, 'ud' for up-down.
        
    Returns:
        GeoDataFrame: A new GeoDataFrame with flipped geometries.
    """
    new_gdf = gdf.copy()
    minx, miny, maxx, maxy = new_gdf.total_bounds
    centroid_x = (maxx + minx) / 2
    centroid_y = (maxy + miny) / 2

    # Determine the scale factors based on the direction of the flip
    if direction == 'lr':
        # Flip left-right
        new_gdf['geometry'] = new_gdf['geometry'].apply(
            lambda geom: scale(geom, xfact=-1, yfact=1, origin=(centroid_x, centroid_y))
        )
    elif direction == 'ud':
        # Flip up-down
        new_gdf['geometry'] = new_gdf['geometry'].apply(
            lambda geom: scale(geom, xfact=1, yfact=-1, origin=(centroid_x, centroid_y))
        )
    else:
        raise ValueError("Invalid direction specified. Use 'lr' for left-right or 'ud' for up-down.")

    return new_gdf


def load_cell_boundary_parquet(sample, data_dir='data'):
    """
    Loads cell boundary data from parquet files for different platforms and merges it with cell level information.

    This function handles the loading and processing of cell boundaries differently based on the 'sample' platform
    (e.g., 'xenium', 'merscope', 'cosmx'). It also merges the resulting geometries with cell level data such as 'core'
    and 'tissue_type' for further analyses.

    Parameters:
        sample (str): The sample identifier which includes the platform and potentially other identifiers that
                      dictate the file paths and processing logic.

    Returns:
        GeoDataFrame: A GeoDataFrame containing merged cell geometry and cell level information.
    """

    # Xenium platform processing
    if 'xenium' in sample:
        df_b = pd.read_parquet(f'{data_dir}/{sample}/cell_boundaries.parquet')

        # Decode object dtype columns to UTF-8 if not from the year 2024
        if '2024' not in sample:
            for col in df_b.select_dtypes(include=['object']).columns:
                df_b[col] = df_b[col].str.decode('utf-8')

        grouped = df_b.groupby('cell_id')
        cell_ids = []
        polygons = []

        # Creating polygons from cell boundary vertices
        for cell_id, group in grouped:
            polygon = Polygon(zip(group['vertex_y'], group['vertex_x']))
            polygons.append(polygon)
            cell_ids.append(cell_id)

        gdf = gpd.GeoDataFrame({'cell_id': cell_ids, 'geometry': polygons}, crs='EPSG:4326')

        # Load cell level data and merge
        df_c = pd.read_parquet(f'{data_dir}/cell_level_csv/{sample}_cell_level.parquet.gzip', engine='pyarrow')[['cell_id', 'core', 'tissue_type']]
        gdf_join = pd.merge(gdf, df_c, on='cell_id', how='inner')

    # Merscope platform processing
    elif 'merscope' in sample:
        gdf = gpd.GeoDataFrame()

        # Process each region's cell boundary data
        for i, folder in enumerate(glob.glob(f'{data_dir}/{sample}/region_*')):
            print(f"region_{i}")
            file = f"{folder}/cell_boundaries.parquet"
            gdf_b = gpd.read_parquet(file)
            gdf_b['EntityID'] = gdf_b['EntityID'].apply(lambda x: f"{x}_region_{i}")
            if gdf_b['EntityID'].duplicated().any():
                print(f"region_{i} has duplicated EntityID")
            gdf = pd.concat([gdf, gdf_b], ignore_index=True)

        # Select mid Z level to standardize data layer
        mid_index = len(sorted(gdf.ZLevel.unique())) // 2
        mid_num = sorted(gdf.ZLevel.unique())[mid_index]
        print(f'Mid Z level: {mid_num}')
        gdf = gdf.loc[gdf['ZLevel'] == mid_num].explode()
        
        # Load cell level data and merge
        df_c = pd.read_parquet(f'{data_dir}/cell_level_csv/{sample}_cell_level.parquet.gzip', engine='pyarrow')[['cell_id', 'core', 'tissue_type']]
        gdf_join = pd.merge(gdf, df_c, on='cell_id', how='inner')

    # Cosmx platform processing
    elif 'cosmx' in sample:
        gdf = gpd.GeoDataFrame()
        fovs = [x for x in range(1, 151)] if 'htma' in sample else [x for x in range(1, 177)]

        # Process each field of view (FOV)
        for fov in fovs:
            try:
                original_fov = fov
                fov = str(fov).zfill(3)
                parquet_file = f'{data_dir}/{sample}/segmentation_by_fov_parquet/FOV_{fov}.parquet.gzip'
                gdf_fov = gpd.read_parquet(parquet_file)
                print(fov, len(gdf_fov))
                gdf_fov['fov'] = fov
                gdf_fov['cell_id'] = gdf_fov['value'].apply(lambda x: f"c_1_{original_fov}_{x}")
                gdf_fov['fov'] = gdf_fov['fov'].astype(int)
                gdf = pd.concat([gdf, gdf_fov], ignore_index=True)
            except:
                continue

        # Load cell level data and merge
        df_c = pd.read_parquet(f'{data_dir}/cell_level_csv/{sample}_cell_level.parquet.gzip', engine='pyarrow')[['cell_id', 'core', 'tissue_type', 'fov']]
        gdf_join = pd.merge(df_c, gdf, on=['cell_id', 'fov'], how='left')

    return gdf_join


def plot_cell_filtration(
        gdf_mask, figwidth, scale_bar,
        xy_range=False, area_bounds={},
        markersize=1, boundary_color='white', save=False, fname=''):
    """
    Plot a GeoDataFrame representing cell filtration data with custom color mapping.
    
    Parameters:
        gdf_mask (GeoDataFrame): GeoDataFrame with a 'Keep' column indicating whether to keep (1) or drop (0) a cell.
        figwidth (int): Width of the figure in inches.
        scale_bar (matplotlib.artist.Artist): Scale bar to be added to the plot.
        xy_range (bool): If True, use custom limits for the plot axes.
        area_bounds (dict): Dictionary with keys 'xmin', 'xmax', 'ymin', 'ymax' for setting plot limits.
        markersize (int): Size of the markers in the plot.
        save (bool): If True, save the plot to a file.
        fname (str): Filename or path to save the plot if `save` is True.
    """
    # Define custom colors for the 'Keep' values
    cmap = ListedColormap(['#D55E00', '#009E73'])
    norm = BoundaryNorm([0, 0.5, 1], cmap.N)  # Define the boundaries for "Drop" (0) and "Keep" (1)

    fig, ax = plt.subplots(figsize=(figwidth, figwidth))
    fig.patch.set_facecolor('black')

    # Plot the GeoDataFrame with the custom colormap
    scatter = gdf_mask.plot(column='Keep', aspect=1, markersize=markersize, ax=ax,
                            legend=False, edgecolor=boundary_color, cmap=cmap, norm=norm,
                            linewidth=1.5)

    # Optionally add a scale bar
    ax.add_artist(scale_bar)

    # Set custom axis limits if specified
    if xy_range:
        ax.set_xlim(area_bounds['xmin'], area_bounds['xmax'])
        ax.set_ylim(area_bounds['ymin'], area_bounds['ymax'])

    # Hide axis details
    plt.axis('off')
    plt.show()

    # Save the figure if required
    if save:
        fig.savefig(fname, format='png', dpi=200)
        fig.savefig(fname.replace('png', 'eps'), format='eps', dpi=200)


def plot_cell_transcripts(
        gdf_t, gene, figwidth, scale_bar,
        xy_range=False, area_bounds={},
        markersize=1, only_keep=False, save=False, fname=''):
    """
    Plot transcripts for a specified gene within a geographic data frame.
    
    Parameters:
        gdf_t (GeoDataFrame): GeoDataFrame containing transcript information.
        gene (str): Specific gene to highlight in the plot.
        figwidth (float): Width of the figure in inches.
        scale_bar (ScaleBar): Scale bar to be added to the plot, indicating the scale.
        xy_range (bool): If True, use custom limits for the plot axes.
        area_bounds (dict): Dictionary specifying the 'xmin', 'xmax', 'ymin', 'ymax' for setting plot limits.
        markersize (int): Size of the markers in the plot.
        only_keep (bool): If True, only plot data for transcripts marked as 'Keep'.
        save (bool): If True, save the plot to a file.
        fname (str): Filename or path to save the plot if `save` is True.
    
    Outputs:
        A plot is displayed and optionally saved to a file.
    """
    # Copy the data to avoid changing the original DataFrame
    gdf_t_plot = gdf_t.copy()

    # Filter for only 'Kept' transcripts if only_keep is True
    if only_keep:
        gdf_t_plot = gdf_t.loc[gdf_t['Keep'] == 1]

    # Filter transcripts for the specific gene
    gdf_t_plot_gene = gdf_t_plot.loc[gdf_t_plot['gene'] == gene]

    # Set up the plot
    fig, ax = plt.subplots(figsize=(figwidth, figwidth))
    fig.patch.set_facecolor('black')

    # Plot all transcripts
    gdf_t_plot.plot(aspect=1, markersize=markersize * 0.1, ax=ax, legend=False, color='limegreen')

    # Add scale bar to the plot
    ax.add_artist(scale_bar)

    # Plot specific gene transcripts
    gdf_t_plot_gene.plot(aspect=1, markersize=markersize, ax=ax, legend=True, color='blue')

    # Set custom axis limits if specified
    if xy_range:
        ax.set_xlim(area_bounds['xmin'], area_bounds['xmax'])
        ax.set_ylim(area_bounds['ymin'], area_bounds['ymax'])

    plt.axis('off')
    plt.show()

    # Save the figure if required
    if save:
        fig.savefig(fname, format='png', dpi=200)
        fig.savefig(fname.replace('png', 'eps'), format='eps', dpi=200)


def translate_to_bbox(original_gdf, bound_box_values_um, des_img_shape):
    # Calculate the bounds of the current GeoDataFrame
    minx = bound_box_values_um[0]
    miny = bound_box_values_um[1]
    maxx = bound_box_values_um[2]
    maxy= bound_box_values_um[3]


    # Calculate translation factors to move minx, miny to 0, 0
    trans_x = -minx
    trans_y = -miny

    # Translate the geometries
    translated_gdf = original_gdf.copy()
    translated_gdf['geometry'] = translated_gdf['geometry'].translate(trans_x, trans_y)

    # Calculate the current max extents after translation
    new_maxx, new_maxy = maxx + trans_x, maxy + trans_y

    # Calculate scaling factors
    scale_x = des_img_shape[1] / new_maxx
    scale_y = des_img_shape[0] / new_maxy

    # Apply scaling
    translated_gdf['geometry'] = translated_gdf['geometry'].scale(xfact=scale_x, yfact=scale_y, origin=(0, 0))

    return translated_gdf
