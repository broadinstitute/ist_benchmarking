"""
Author: 
Huan Wang (whuan@broadinstitute.org)
OPP, Broad Institute of MIT and Harvard

Functions:
    calculate_qc_metric()
    convert_to_px()
    correct_tissue_names()
    create_circular_mask()
    data_loader()
    df_2_gdf()
    get_cell_by_gene()
    get_df_color()
    get_gdf_core()
    get_gdf_core_from_polygon()
    get_gene_type()
    get_processed()
    imshow()
    invert_y()
    load_gis_csv()
    load_shape_file()
    name_parser()
    save_html()
    seg_overlay()
    sim_overlay()
    standardize_data()
    subset_on_condition()
    transcript_loader()
    vecterize()


"""
import os
import sys
import pandas as pd
import scanpy as sc
import anndata as ad
import squidpy as sq
import subprocess
from collections import Counter
import json
import glob
from geopandas import GeoDataFrame
from shapely.geometry import Point,Polygon,box
import geopandas as gpd
from shapely.affinity import scale
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.ops import unary_union
from matplotlib import pyplot as plt
import numpy as np
from skimage import io
from skimage.transform import rescale



LEANER_COL_DICT = {'merscope':['transcript_id','cell_id','gene','global_x','global_y','fov'],
                   'xenium':['transcript_id','cell_id','feature_name','x_location','y_location']}


def calculate_qc_metric(adata, percent_top=(50, 100, 200, 250), min_counts=10, min_cells=5):
    """calculate qc metrics using squidpy
    """
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],percent_top=percent_top, inplace=True)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_counts)

    return adata


def convert_to_px(conversion, mat):
    """This function transforms coordinates from um to px given a conversion matrix 
    
    Parameters
    -----------
    conversion: conversion matrix from vizgen the file is called micron_to_mosaic_pixel_transform
    mat: pandas df of x and y coordinates
    Returns
    -----------
    xcoord_transformed: X coordinate with new origin and scaled
    ycoord_transformed: Y coordinate with new origin and scaled
    """
        
    mat = np.array(mat)
    new_col = np.ones((mat.shape[0], 1))
    all_data = np.concatenate((mat, new_col), 1)
    
    # affine transform
    transformed_mat = (conversion @ np.array(all_data.T)).T
    
    return [transformed_mat[[0]], transformed_mat[[1]]]


def correct_platform_panel(df, col_name):
    correct_platform_panel = {'cosmx_multitissue':'CosMx,1k',
                              'merscope_breast':'MERSCOPE,breast',
                              'merscope_lung':'MERSCOPE,lung',
                              'xenium_breast':'Xenium,breast',
                              'xenium_lung':'Xenium,lung',
                              'xenium_panhuman':'Xenium,multi-tissue'}
    df[col_name] = df[col_name].str.lower()
    df[col_name] = df[col_name].replace(correct_platform_panel)
    return df


def correct_tissue_names(sample, df):

    htma_correct_names = {'Bladder':'BlC',
                        'BrC':'BrC', 
                        'CRC':'CRC',
                        'HNSCC':'HNSCC',
                        'Lymph node':'Lymph node',
                        'LN':'Lymph node',
                        'Marker':'Marker',
                        'MARKER':'Marker',
                        'Melanoma':'Mel',
                        'NSCLC':'NSCLC',
                        'OV':'OvC',
                        'Tonsil':'Tonsil'}


    normal_correct_names = {'Bladder':'Bladder', 
                            'Breast':'Breast',
                            'Colon':'Colon',
                            'Heart':'Heart',
                            'Kidney':'Kidney',
                            'Liver':'Liver',
                            'Lung':'Lung',
                            'Lymph node':'Lymph node',
                            'LN':'Lymph node',
                            'Marker':'Marker',
                            'MARKER':'Marker',                            
                            'Ovarian':'Ovary',
                            'Pancreas':'Pancreas',
                            'Prostate':'Prostate',
                            'Renal':'Renal pelvis',
                            'Skin':'Skin',
                            'Spleen':'Spleen',
                            'Thyroid':'Thyroid',
                            'Tonsil':'Tonsil'}
    
    tumor2_correct_names = {
        'Bladder':'BlC',
        'Colon CA':'CRC',
        'DCIS grade 2':'Breast non-invasive DCIS_2',
        'DCIS grade 3':'Breast non-invasive DCIS_3',
        'Invasive breast':'Breast invasive',
        "Hodgkin's Lymphoma":'Lymphoma Hodgkin',
        "Non-Hodgkin's Lymphoma":'Lymphoma non-Hodgkin',    
        'LN B-cell lymphoma':'Lymphoma LN B cell',            
        'Kidney':'Kidney cancer',
        'Liposarcoma':'Liposarcoma',
        'Liver':'Liver cancer',
        'Marker Liver':'Marker normal liver',
        'Melanoma':'Mel',
        'Ovarian':'OvC',
        'Pancreas':'Pancreas cancer',
        'Prostate':'Prostate cancer',
        'Renal':'Renal cancer',
        'Marker Spleen':'Marker normal spleen',
        'Squamous cell carcinoma':'SCC',
        'Testes':'Testes cancer',
        'Thyroid':'Thyroid cancer',
        }
    
    if 'htma' in sample:
        df['tissue_type'] = df['tissue_type'].replace(htma_correct_names)
    elif 'normal' in sample:
        df['tissue_type'] = df['tissue_type'].replace(normal_correct_names)
    elif 'tumor2' in sample:
        df['tissue_type'] = df['tissue_type'].replace(tumor2_correct_names)

    return df


def create_circular_mask(h, w, center=None, radius=None):
    """This function will return a circular mask of a given radius in a h/w defined image array"""    

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask   


def data_loader(SAMPLE, modality, data_type='adata'):

    """Data ingestion for spatial transcriptomic data
    """

    if modality=='xenium':
        
        # Load h5 data
        adata = sc.read_10x_h5(filename=f'{SAMPLE}/cell_feature_matrix.h5')

        # Load cell.parquet
        input_file = f'{SAMPLE}/cells.parquet'
        df = pd.read_parquet(input_file, engine='pyarrow')

        if '2024' not in SAMPLE:
            for col in df.select_dtypes(include=['object']).columns:
                df[col] = df[col].str.decode('utf-8')

        # Rotate 90 and flip up-down
        df['x'] = df['y_centroid']
        df['y'] = df['x_centroid']

        # Populate adata.obs
        df.set_index(adata.obs_names, inplace=True)

    elif modality=='merscope':

        df_cell_by_gene = pd.DataFrame()
        df_cell_metadata = pd.DataFrame()

        for fd in sorted(glob.glob(f'{SAMPLE}/region_*/')):

            print (fd)
            df_cell_by_gene_fd = pd.read_csv(f'{fd}/cell_by_gene.csv')
            df_cell_metadata_fd = pd.read_csv(f'{fd}/cell_metadata.csv')

            # Append region str to avoid duplicated cell id between regions
            df_cell_by_gene_fd['cell'] = df_cell_by_gene_fd['cell'].apply(lambda x: f"{x}_{fd.split('/')[-2]}")
            df_cell_metadata_fd['EntityID'] = df_cell_metadata_fd['EntityID'].apply(lambda x: f"{x}_{fd.split('/')[-2]}")

            df_cell_by_gene = pd.concat([df_cell_by_gene, df_cell_by_gene_fd], ignore_index=True)
            df_cell_metadata = pd.concat([df_cell_metadata, df_cell_metadata_fd], ignore_index=True)

        df_cell_by_gene = df_cell_by_gene.rename(columns={'cell':'cell_id'})
        df_cell_metadata = df_cell_metadata.rename(columns={'EntityID':'cell_id'})

        df_cell_by_gene.to_csv(f'{SAMPLE}/cell_by_gene.csv', index=False)
        df_cell_metadata.to_csv(f'{SAMPLE}/cell_metadata.csv', index=False)

        # Load h5 data
        adata = sc.read_csv(filename=f'{SAMPLE}/cell_by_gene.csv')

        # Load cell.csv
        df = pd.read_csv(f'{SAMPLE}/cell_metadata.csv')
        df = pd.merge(df, df_cell_by_gene, on='cell_id', how='inner')
        df = df.rename(columns={'center_x':'x_centroid',
                                'center_y':'y_centroid',
                                'transcript_count':'transcript_counts'})

        df['x'] = df['x_centroid']
        df['y'] = df['y_centroid']

        # Populate adata.obs
        df['cell_id'] = df['cell_id'].astype(str)
        df = df.set_index('cell_id', drop=False)
        df = df.rename_axis(None)

    elif modality=='cosmx':

        # Load raw csv (previously saved from tile db data)
        df =  pd.read_csv(f'{SAMPLE}/raw_data.csv')

        df['x_centroid'] = df['y']
        df['y_centroid'] = df['x']

        df['x'] = df['x_centroid']
        df['y'] = df['y_centroid']
        df.index=df['cell_id']
   
        df = df.rename(columns={'nCount_RNA':'transcript_counts'})
        adata = ad.AnnData(obs = df)
         
    else:
        print (f'This modality is not supporetd yet: {modality}')

    if data_type == 'adata':
        # Add apatial component
        adata.obs = df.copy()
        adata.obs['cell_id'] = adata.obs.index
        adata.obsm["spatial"] = adata.obs[["x", "y"]].copy().to_numpy()
        return adata
    else:
        return df



def df_2_gdf(df, x_col, y_col, crs="EPSG:4326", drop_xy=False):
    """convert DataFrame to GeoDataFrame
    """
    geometry = [Point(xy) for xy in zip(df[x_col], df[y_col])]
    if drop_xy:
        df = df.drop([x_col, y_col], axis=1)
    gdf = GeoDataFrame(df, crs=crs, geometry=geometry)

    return gdf


def get_cell_by_gene(df_t, sample):

    # Pivot the data to get gene counts
    if 'cosmx' in sample:
        df_pivot = df_t.pivot_table(index='cell_id', columns='gene', values='fov', aggfunc='count', fill_value=0)
    else:
        df_pivot = df_t.pivot_table(index='cell_id', columns='gene', values='transcript_id', aggfunc='count', fill_value=0)

    # Get unique cell_id with their corresponding core and tissue_type
    cell_info = df_t[['cell_id', 'core', 'tissue_type']].drop_duplicates()

    # Merge the cell_info with the pivoted data
    result = cell_info.merge(df_pivot, on='cell_id').set_index('cell_id')

    # Keep cell_id as one column
    result['cell_id'] = result.index

    # Remove UNASSIGNED cells
    result = result.loc[result['cell_id'] != 'UNASSIGNED']

    # Sort the DataFrame by cell_id
    sorted_result = result.sort_index()

    return sorted_result


def get_df_color(df, cmap='coolwarm'):
    """
    Apply a color gradient style to a DataFrame.

    This function takes a DataFrame and applies a color gradient
    based on the values in numeric columns. The default color map used is 'coolwarm', but
    this can be customized with any colormap available in matplotlib.

    Parameters:
    df (pandas.DataFrame): The DataFrame to which the style will be applied.
    cmap (str, optional): The colormap name to use for the gradient. Defaults to 'viridis'.

    Returns:
    pandas.io.formats.style.Styler: A Styler object which holds the styled DataFrame.

  
    """
    return df.style.background_gradient(cmap=cmap)




def get_gdf_core(csv_points, csv_sample_info, scaling_factor, radius_um, points_src, drop_xy=False):

    """get core information from point csv and convert to GeoDataFrame
    """
    df_points = load_gis_csv(csv_points, points_src)
    df_sample = pd.read_csv(csv_sample_info)
    df_join = pd.merge(df_points, df_sample, left_on='id', right_on='core', how='inner')
    df_join['x'] = df_join['x'] * scaling_factor
    df_join['y'] = df_join['y'] * scaling_factor
    df_join['core'] = df_join['core'].apply(lambda x: str(int(x)))
    gdf = df_2_gdf(df_join, 'x', 'y', drop_xy=drop_xy)
    gdf['geometry'] = gdf['geometry'].buffer(radius_um)
    if '2024_cosmx' in csv_points:
        gdf['geometry'] = gdf['geometry'].envelope

    return gdf


def get_gdf_core_from_polygon(shapefile, csv_sample_info, scaling_factor):

    """get core information from poluygon shapefile and convert to GeoDataFrame
    """
    gdf_polygon = gpd.read_file(shapefile)
    df_sample = pd.read_csv(csv_sample_info)
    gdf_join = pd.merge(gdf_polygon, df_sample, left_on='id', right_on='core', how='left')

    # Scale the geometry
    gdf_join['geometry'] = gdf_join['geometry'].apply(
        lambda geom: scale(geom, 
                        xfact=scaling_factor,
                        yfact=scaling_factor,
                        zfact=scaling_factor,
                        origin=(0,0))
        )
    gdf_join['core'] = gdf_join['core'].apply(lambda x: str(int(x)))

    return gdf_join


def get_gene_type(g):
    g = str(g)
    if 'BLANK' in g: # Xenium
        return 'blank'
    elif 'Blank' in g:  # MERSCOPE
        return 'blank'
    elif 'NegControlCodeword' in g: # Xenium
        return 'neg_control_codeword'
    elif 'NegControlProbe' in g: # Xenium
        return 'neg_control_probe'
    elif 'Negative' in g: # CosMX
        return 'neg_control_probe'
    elif 'SystemControl' in g: # CosMX
        return 'sys_control'
    else:
        return 'gene'
    

def get_merscope_region_meta(sample, region):
    folder = f'data/{sample}/region_{region}'
    with open(f'{folder}/manifest.json', 'r') as json_file:
        # Use json.load() to parse the JSON data
        data = json.load(json_file)
    return data['microns_per_pixel'], data['bbox_microns']


def get_processed(sample, data_type, fast=False):
    if data_type == 'cell_by_gene':
        if fast:
            df = pd.read_parquet(f'data/cell_by_gene_csv/{sample}_cell_by_gene.parquet.gzip', engine='pyarrow')
        else:
            df = gpd.read_parquet(f'data/cell_by_gene_csv/{sample}_cell_by_gene.parquet.gzip')
    elif data_type == 'cell_level':
        if fast:
            df = pd.read_parquet(f'data/cell_level_csv/{sample}_cell_level.parquet.gzip', engine='pyarrow')
        else:
            df = gpd.read_parquet(f'data/cell_level_csv/{sample}_cell_level.parquet.gzip')
    elif data_type == 'transcript_level':
        if fast:
            df = pd.read_parquet(f'data/transcript_level_csv/{sample}_transcript_level.parquet.gzip', engine='pyarrow')
        else:
            df = gpd.read_parquet(f'data/transcript_level_csv/{sample}_transcript_level.parquet.gzip')
    elif data_type == 'gene_level':
        df = pd.read_csv(f'data/gene_level_csv/gene_level_csv_{sample}.csv', engine='pyarrow')
    else:
        print(f'{data_type} data or for {sample} is not available')

    return df
        


def get_qced_cell_id(qc_count_threshold, qc_unique_gene_threshold):

    d_perc = {}
    df_m = pd.read_parquet(f'data/single_cell_metrics/all_cell_level.parquet.gzip', engine='pyarrow')
    df_mat = df_m.copy()
    print (f'\nbefore QC: {len(df_mat)}')
    d_before = {**Counter(df_m['sample'])}
    df_mat = df_mat.loc[df_mat['transcript_counts'] > qc_count_threshold]
    print (f'after QC using transcript count per cell: {len(df_mat)}')
    df_mat = df_mat.loc[df_mat['unique_genes'] > qc_unique_gene_threshold]
    print (f'after QC using unique genes per cell: {len(df_mat)}')
    d_after = {**Counter(df_mat['sample'])}
    print (f'Good quality cells: {round(len(df_mat) * 100/len(df_m),1)}')
    for i in d_after.keys():
        d_perc[i] = round(d_after[i] * 100 /d_before[i], 2)
    return df_mat.cell_id.to_list(), d_perc


def imshow(src, resize=True, vlim=True, q0=0.01, q1=0.99, figwidth=10, tt='', cmap='viridis', save =False, fname=''):
    """In-house imshow function to show numpy array as image

    Parameters
    -----------
    src: image data read as numpy array or image file name in string format
    resize: bool, downsample when True
    vlim: bool, set up vis range when True
    figwidth: float
    cmap: color map, default viridis.

    """
    if isinstance(src, str):
        img = io.imread(src)
        print(src)
    else:
        img = src
    print("shape:", img.shape)
    figheight = figwidth * img.shape[0] / img.shape[1]
    if resize:
        img = rescale(img, 0.1, anti_aliasing=False)
        print("resized shape:", img.shape)
    g = plt.figure(figsize=(figwidth, figheight))
    if vlim:
        plt.imshow(
            img, vmin=np.quantile(img, q0), vmax=np.quantile(img, q1), cmap=cmap
        )
    else:
        plt.imshow(img, cmap=cmap)
    plt.title(tt)
    plt.axis('off')
    plt.show()
    if save:
        g.savefig(fname, format='png', dpi=200)



def invert_y(geometry):
    """Function to modify y-values to -y for a Polygon or Point"""
    if geometry.geom_type == 'Polygon':
        exterior = geometry.exterior.coords
        new_exterior = [(x, -y) for x, y in exterior]
        interiors = geometry.interiors
        new_interiors = []
        for interior in interiors:
            new_interior = [(x, -y) for x, y in interior.coords]
            new_interiors.append(new_interior)
        return Polygon(new_exterior, new_interiors)
    elif geometry.geom_type == 'Point':
        x, y = geometry.coords[0]
        return Point(x, -y)
    else:
        # Handle MultiPolygons or other geometries if needed
        return geometry



def load_gis_csv(csv_file, points_src='image'):    
    """This function convert GIS csv file to df with x and y cols.
    """
    
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.lower()
    if points_src=='image':
        df['y'] = df['y'] * -1
    elif points_src=='points':
        df['y'] = df['y']
    df = df.sort_values(by=['id'])
    df = df.reset_index(drop=True)
    df = df[['id','x','y']]
    return df


def load_shape_file(shapefile):
    gdf = gpd.read_file(shapefile)
    gdf['geometry'] = gdf['geometry'].apply(invert_y)
    return gdf

def name_parser(sample):
    platform = sample.split('_')[-3]
    panel = sample.split('_')[-2]
    tma = sample.split('_')[-1]
    return platform, panel, tma


def save_html(ipynb, suffix=''):
    """Save python notebook as html file

    Parameters
    -----------
    ipynb: notebook file name
    suffix: string

    """

    # save html
    cmd = 'jupyter nbconvert --to html ' + ipynb.split('/')[-1]
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    p.communicate()
    # rename
    html = f"{ipynb.split('.ipynb')[0]}.html"
    new_html = f"{html.split('.html')[0]}{suffix}.html"
    if not os.path.exists(new_html):
        os.rename(html, new_html)


def seg_overlay(img_dapi,img_boundary,q0=0, q1=1, figwidth=10,cmap='gist_gray'):
    rgb1 = img_dapi.copy()
    rgb1[img_boundary > 0] = np.max(img_dapi)
    plt.axis('off')
    imshow(rgb1, resize=False, q0=q0, q1=q1, figwidth=figwidth,cmap=cmap)


def show_n(
    n,
    imgs,
    figsize,
    titles=['','',''],
    vlim=True,
    q_low=0.05,
    q_up=0.95,
    cmap='bone',
):
    """Shown two image arrays side by side

    Parameters
    -----------
    n: number of arrays
    imgs: list of image data read as numpy array
    figsize: tuple, such as (10,10)
    titles: list of title of image array
    vlim: bool, True to apply auto vlim
    q_low: quantile to set lower vlim
    q_up: quantile to  set upper vlim
    cmap: color map, default viridis.

    """
    vmins = [0] * n
    vmaxs = [0] * n
    
    fig, ax = plt.subplots(1, n, figsize=figsize)
    for i in range(n):
        
        if vlim:
            vmins[i] = np.quantile(imgs[i], q_low)
            vmaxs[i] = np.quantile(imgs[i], q_up)
            ax[i].imshow(imgs[i], vmin=vmins[i], vmax=vmaxs[i], cmap=cmap)
        else:
            ax[i].imshow(imgs[i], cmap=cmap)
        ax[i].set_title(titles[i])
    plt.show()  


def sim_overlay(img_dapi, vector_mask, figwidth, q0, q1, cmap, scale_bar, xy_range=False, x_range=[], y_range=[],lw=1, save=False, fname=''):    

    figheight = figwidth * img_dapi.shape[0] / img_dapi.shape[1]

    fig, ax = plt.subplots(1, 1, figsize=(figwidth, figheight))

    ax.imshow(
        img_dapi, vmin=np.quantile(img_dapi, q0), vmax=np.quantile(img_dapi, q1), cmap=cmap)

    vector_mask.plot(aspect=1, ax=ax, legend=False,
               edgecolor='white',
               facecolor='none',
               linewidth=lw)

    ax.add_artist(scale_bar)
    if xy_range:
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)
    plt.axis('off')
    plt.show()
    
    if save:
        fig.savefig(fname, format='png', dpi=200)
        fig.savefig(fname.replace('png','eps'), format='eps', dpi=200)


def standardize_data(adata, cluster_key, string1, string2, leiden_resolution=1.0, dpi=200, marker_size=2):
    """standard preprocessing using squidpy
    """

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_resolution)

    sc.pl.umap(
        adata,
        color=[
            "leiden",
        ],
        save=string1,
    )
    plt.show()

    sq.pl.spatial_scatter(
        adata,
        color=cluster_key,
        shape=None,
        size=marker_size,
        dpi=dpi,
        save=string2,
    )
    plt.show()


def subset_on_condition(adata, gdf_core, condition_col, condition_value_list, extra_variable=[],positive=True, xy_flip=False):
    """subset adata by spatial join, requires GeoDataFrame on the right side.
    """
    gdf_adata = adata.obs
    if xy_flip:
        gdf_adata = df_2_gdf(gdf_adata, 'y_centroid', 'x_centroid')
    else:
        gdf_adata = df_2_gdf(gdf_adata, 'x_centroid', 'y_centroid')

    if positive:
        gdf_sub = gdf_core.loc[gdf_core[condition_col].isin(condition_value_list)]
        print (f"condition: {condition_col}, value: {condition_value_list}")
    else:
        gdf_sub = gdf_core.loc[~gdf_core[condition_col].isin(condition_value_list)]
        print (f"condition: {condition_col}, value: not {condition_value_list}")

    adata_sub = gdf_sub[
        ['core','tissue_type','geometry'] + extra_variable
        ].sjoin(gdf_adata, how='right', op='intersects'
                ).drop(columns='index_left').dropna(subset=['core'])

    print (f"original: {len(gdf_adata)}, subset: {len(adata_sub)}, pct: {round(100 * len(adata_sub)/len(gdf_adata), 2)}%")

    return adata_sub, gdf_sub


def transcript_loader(SAMPLE):

    df_t= pd.DataFrame()
    sample = SAMPLE.split('/')[-1]
    modality = sample.split('_')[-3]
    

    if modality=='xenium':
        input_file = f'{SAMPLE}/transcripts.parquet'
        df_t = pd.read_parquet(f'{SAMPLE}/transcripts.parquet')
        # df_t = df_t[LEANER_COL_DICT[modality]]
        df_t = df_t.rename(columns={'feature_name':'gene',
                                    'x_location':'global_x',
                                    'y_location':'global_y'})
        
        for col in df_t.select_dtypes(include=['object']).columns:
            df_t[col] = df_t[col].str.decode('utf-8')

    elif modality=='merscope':

        df_t = pd.DataFrame()
 
        for fd in sorted(glob.glob(f'{SAMPLE}/region_*/')):

            print (fd)
            df_t_fd = pd.read_csv(f'{fd}/detected_transcripts.csv',
                                  engine="pyarrow")
            
            # Append region str to avoid duplicated cell id between regions
            df_t_fd['cell_id'] = df_t_fd['cell_id'].apply(lambda x: f"{x}_{fd.split('/')[-2]}")
            
            df_t = pd.concat([df_t, df_t_fd], ignore_index=True)

        # df_t = df_t[LEANER_COL_DICT[modality]]

    elif modality=='cosmx':

        df_t = pd.read_csv(f'{SAMPLE}/detected_transcripts.csv',
                                engine="pyarrow")
        df_t = df_t.rename(columns={'target':'gene'})

        # df_t = df_t[LEANER_COL_DICT[modality]]
    df_t['modality'] = modality

    print (len(df_t))

    return df_t


def vectorize(arr):
    # Mask out zeros to focus only on non-zero values for processing
    mask = arr != 0
    results = ({'properties': {'value': v}, 'geometry': s}
               for i, (s, v) in enumerate(shapes(arr, mask=mask, connectivity=8)))

    # Convert shapes and values to GeoDataFrame
    geoms = []
    values = []
    for result in results:
        geom = shape(result['geometry'])
        val = result['properties']['value']
        geoms.append(geom)
        values.append(val)

    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame({'value': values, 'geometry': geoms}, crs="EPSG:4326")
    gdf['value'] = gdf['value'].astype(int)

    return gdf
