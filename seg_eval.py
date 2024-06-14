"""
Author: 
Huan Wang (whuan@broadinstitute.org)
OPP, Broad Institute of MIT and Harvard

Functions:
    calculate_seg_eval_metrics()
    calculate_tp_fn_fp_from_gdf()
"""

import geopandas as gpd

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
