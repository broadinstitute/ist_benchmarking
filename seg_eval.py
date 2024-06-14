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
    translate_to_bbox()
"""

import geopandas as gpd
import numpy as np
from shapely.affinity import scale
import numpy as np
import matplotlib.pyplot as plt


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
