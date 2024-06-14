
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from scipy.ndimage import binary_dilation
from skimage import measure



def xy_df_to_tif(df, shape, tif_filename):
    
    # Create an empty array with the shape of the output image
    image_array = np.zeros(shape, dtype=np.uint8)

    # Update the array based on x, y coordinates from the CSV, assuming columns are named 'x' and 'y'
    for index, row in df.iterrows():
        x, y = int(row['x']), int(row['y'])
        if 0 <= x < shape[1] and 0 <= y < shape[0]:  # Check bounds
            image_array[y, x] = 1  # Set the value to 1

    # Create a new raster data source
    transform = from_origin(0, 0, 1, 1)  # Top left corner coordinates and cell size
    with rasterio.open(
        tif_filename, 'w',
        driver='GTiff',
        height=shape[0],
        width=shape[1],
        count=1,
        dtype=image_array.dtype,
        crs='+proj=latlong',
        transform=transform,
    ) as dst:
        dst.write(image_array, 1)

    return image_array


def thicker_dot(array, dilation_size=2):

    structuring_element = np.ones((dilation_size, dilation_size), dtype=int)

    # Perform binary dilation
    dilated_array = binary_dilation(array, structure=structuring_element).astype(int)
    return dilated_array


def calculate_tp_fn_fp_from_img(img_mask, img_gt):
    """
    Calculate the number of true positive (TP), false negative (FN), and false positive (FP) cells 
    in segmentation masks.

    Parameters:
    img_mask (numpy.ndarray): A 2D array where each value represents a segmented cell, with 0 meaning no cell.
    img_gt (numpy.ndarray): Ground truth, A 2D array of the same shape as img_mask, 
            where 0 means no cell and 1 means the presence of a cell.

    Returns:
    tuple: A tuple containing the number of true positive cells (int), false negative cells (int), 
           and false positive cells (int).
    
    True Positive (TP): Cells in img_mask that overlap with cells in img_gt.
    False Negative (FN): Cells present in img_gt but missing in img_mask.
    False Positive (FP): Cells present in img_mask but not in img_gt.
    """
    # Set of unique labels in img_mask excluding the background (0)
    mask_labels = set(np.unique(img_mask)) - {0}
    # Set of unique labels in img_mask where img_gt equals 1, excluding the background (0)
    gt_cells = set(np.unique(img_mask[img_gt == 1])) - {0}
    
    # True Positive Cells
    tp_cells = mask_labels & gt_cells
    tp = len(tp_cells)
    
    # False Negative Cells
    fn_cells = gt_cells - mask_labels
    fn = len(fn_cells)
    
    # False Positive Cells
    fp_cells = mask_labels - gt_cells
    fp = len(fp_cells)
    
    return tp, fn, fp


def sc_extractor_skimage(img_mask, img_src):
    """
    Extracts the mean intensity for each cell from a multi-dimensional image (2D or 3D) based on a segmentation mask
    using skimage.measure to efficiently calculate properties.

    Parameters:
        img_mask (ndarray): An array where 0 represents background and other integer values indicate cell IDs.
        img_src (ndarray): A multi-dimensional image array where each band/channel is a different imaging channel if 3D.

    Returns:
        DataFrame: A pandas DataFrame with columns for cell_id and mean intensity for each channel.
    """
    # Check if img_src is 2D and convert to 3D if needed
    if img_src.ndim == 2:
        img_src = np.expand_dims(img_src, axis=-1)  # Convert to a single-channel 3D array

    # Create an empty list to hold data for each channel
    data = []
    for channel in range(img_src.shape[2]):
        # Calculate properties using regionprops_table
        props = measure.regionprops_table(img_mask, intensity_image=img_src[:, :, channel],
                                          properties=['label', 'mean_intensity'])

        # Convert to DataFrame
        df_channel = pd.DataFrame(props)
        df_channel.rename(columns={'mean_intensity': f'channel_{channel + 1}'}, inplace=True)
        
        # Merge or concatenate channel data
        if channel == 0:
            df = df_channel
        else:
            df = df.merge(df_channel, on='label', how='outer')

    # Rename columns appropriately
    df.rename(columns={'label': 'cell_id'}, inplace=True)

    return df