"""

Author: 
Huan Wang (whuan@broadinstitute.org)
OPP, Broad Institute of MIT and Harvard

Functions:

Distionarys:

"""

from matplotlib.colors import ListedColormap

CB_COLORS_8 = {'grey':'#C0C0C0',
             'blue_light':'#56B4E9',
             'blue_dark':'#0072B2',
             'green':'#009E73',
             'yellow':'#F0E442',
             'orange':'#E69F00',
             'pink':'#CC79A7',
             'red':'#D55E00',
}

cb_8 = ListedColormap(CB_COLORS_8.values())

CB_COLORS_5 = {'grey':'#C0C0C0',
             'blue_dark':'#0072B2',
             'green':'#009E73',
             'yellow':'#F0E442',
             'red':'#D55E00',
}

cb_5 = ListedColormap(CB_COLORS_5.values())

CORRECT_PLATFORM_PANEL = {'cosmx_multitissue':'CosMx,1k',
                            'merscope_breast':'MERSCOPE,breast',
                            'merscope_lung':'MERSCOPE,lung',
                            'xenium_breast':'Xenium,breast',
                            'xenium_lung':'Xenium,lung',
                            'xenium_panhuman':'Xenium,multi-tissue'}

sample_color = {'Xenium,breast':'#0072B2',
              'Xenium,multi-tissue':'#56B4E9',
              'Xenium,lung':'#C0C0C0',
              'MERSCOPE,breast':'#D55E00',
              'MERSCOPE,lung':'#F0E442',
              'CosMx,1k':'#009E73',
              'xenium':'#0072B2',
              'merscope':'#D55E00',
              'cosmx':'#009E73',}

seg_cmap = ListedColormap(['black'] + [(0, 0, i/255) for i in range(256)])

CMAPS = ["CMRmap", "CMRmap_r", "gist_stern_r", "gray_r", "turbo", "viridis", "RdBu", "coolwarm", "PiYG"]


# Unit: um
CORE_RADIUS_DICT = {'xenium_breast_htma':350,
                    'xenium_breast_normal':360,
                    'xenium_panhuman_normal':360,
                    'xenium_panhuman_htma':350,
                    'xenium_lung_htma':350,
                    'xenium_lung_normal':360,
                    'merscope_breast_htma':350,
                    'merscope_breast_normal':360,
                    'merscope_breast_htma_round1':350,
                    'merscope_breast_normal_round1':360,
                    'merscope_lung_normal':360,
                    'merscope_lung_htma':350,
                    'cosmx_multitissue_htma':0.4,
                    'cosmx_multitissue_normal':0.45,
                    '2024_xenium_breast_htma':350,
                    '2024_xenium_breast_tumor2':560,
                     '2024_merscope_breast_htma':350,
                    '2024_merscope_breast_tumor2':560,
                    '2024_cosmx_multitissue_htma':0.35,
                    '2024_cosmx_multitissue_tumor2':0.5,}


modality_palette = dict(zip(['XENIUM','MERSCOPE','COSMX'],
                            [CB_COLORS_8['blue_dark'], CB_COLORS_8['red'], CB_COLORS_8['grey']
                             ]))

modality_panel_palette = dict(zip(['COSMX_Multitissue', 
                                   'MERSCOPE_Breast',
                                   'MERSCOPE_Lung',
                                   'XENIUM_Breast', 
                                   'XENIUM_Lung',
                                   'XENIUM_Panhuman'],
                                   [CB_COLORS_8['grey'],
                                    CB_COLORS_8['orange'],
                                    CB_COLORS_8['red'],
                                    CB_COLORS_8['green'],
                                    CB_COLORS_8['blue_light'],
                                    CB_COLORS_8['blue_dark'],
                             ]))

POINTS_SRC_DICT = {'xenium_breast_htma':'image',
                    'xenium_breast_normal':'image',
                    'xenium_panhuman_normal':'image',
                    'xenium_panhuman_htma':'image',
                    'xenium_lung_htma':'image',
                    'xenium_lung_normal':'image',
                    'merscope_breast_htma':'points',
                    'merscope_breast_normal':'points',
                    'merscope_lung_normal':'points',
                    'merscope_lung_htma':'points',
                    'merscope_breast_htma_round1':'points',
                    'merscope_breast_normal_round1':'points',
                    'cosmx_multitissue_htma':'points',
                    'cosmx_multitissue_normal':'points',
                    '2024_xenium_breast_htma':'image',
                    '2024_xenium_breast_tumor2':'image',
                     '2024_merscope_breast_htma':'points',
                    '2024_merscope_breast_tumor2':'points',
                    '2024_cosmx_multitissue_htma':'points',
                    '2024_cosmx_multitissue_tumor2':'points',}

SAMPLES = ['xenium_breast_htma',
        'xenium_breast_normal',
        'xenium_panhuman_htma',
        'xenium_panhuman_normal',
        'xenium_lung_htma',
        'xenium_lung_normal',
        'merscope_breast_htma',
        'merscope_breast_normal',
        'merscope_lung_htma',
        'merscope_lung_normal',
        'cosmx_multitissue_htma',
        'cosmx_multitissue_normal',
        '2024_xenium_breast_htma',
        '2024_xenium_breast_tumor2',
        '2024_merscope_breast_htma',
        '2024_merscope_breast_tumor2',
        '2024_cosmx_multitissue_htma',
        '2024_cosmx_multitissue_tumor2',
        ]

# Unit: um depretiated
SCALING_FACTOR_DICT = {'xenium_breast_htma':0.2125,
                    'xenium_breast_normal':0.2125,
                    'xenium_panhuman_normal':0.2125,
                    'xenium_panhuman_htma':0.2125,
                    'xenium_lung_htma':0.2125,
                    'xenium_lung_normal':0.2125,
                    'merscope_breast_htma':1,
                    'merscope_breast_normal':1,
                     'merscope_breast_htma_round1':1,
                    'merscope_breast_normal_round1':1,
                    'merscope_lung_normal':1,
                    'merscope_lung_htma':1,
                    'cosmx_multitissue_htma':0.12,
                    'cosmx_multitissue_normal':0.12,
                    '2024_xenium_breast_htma':0.2125,
                    '2024_xenium_breast_tumor2':0.2125,
                    '2024_merscope_breast_htma':1,
                    '2024_merscope_breast_tumor2':1,
                    '2024_cosmx_multitissue_htma':1,
                    '2024_cosmx_multitissue_tumor2':1,
                    }


# Unit: um
PIXEL_TO_UM = {'xenium':0.2125,
               'merscope':0.108,
               'cosmx':0.12}

# Bool
XY_FLIP_DICT = {'xenium_breast_htma':True,
                    'xenium_breast_normal':True,
                    'xenium_panhuman_normal':True,
                    'xenium_panhuman_htma':True,
                    'xenium_lung_htma':True,
                    'xenium_lung_normal':True,
                    'merscope_breast_htma':False,
                    'merscope_breast_normal':False,
                    'merscope_breast_htma_round1':False,
                    'merscope_breast_normal_round1':False,
                    'merscope_lung_normal':False,
                    'merscope_lung_htma':False,
                    'cosmx_multitissue_htma':True,
                    'cosmx_multitissue_normal':True,
                    '2024_xenium_breast_htma':True,
                    '2024_xenium_breast_tumor2':True,
                    '2024_merscope_breast_htma':True,
                    '2024_merscope_breast_tumor2':True,
                    '2024_cosmx_multitissue_htma':True,
                    '2024_cosmx_multitissue_tumor2':True,
                    }

# Unique genes
UNIQUE_GENES_DICT = {'xenium_breast_htma':{'gene': 280, 'blank': 200, 'neg_control_probe': 20, 'neg_control_codeword': 41},
                    'xenium_breast_normal':{'gene': 280, 'blank': 200, 'neg_control_probe': 20, 'neg_control_codeword': 41},
                    'xenium_panhuman_normal':{'gene': 377, 'blank': 103, 'neg_control_probe': 20, 'neg_control_codeword': 41},
                    'xenium_panhuman_htma':{'gene': 377, 'blank': 103, 'neg_control_probe': 20, 'neg_control_codeword': 41},
                    'xenium_lung_htma':{'gene': 289, 'blank': 191, 'neg_control_probe': 20, 'neg_control_codeword': 41}, # TO Update
                    'xenium_lung_normal':{'gene': 289, 'blank': 191, 'neg_control_probe': 20, 'neg_control_codeword': 41}, # TO Update
                    'merscope_breast_htma':{'gene': 255, 'blank': 30},
                    'merscope_breast_normal':{'gene': 255, 'blank': 30},
                    'merscope_lung_normal':{'gene': 220, 'blank': 65},
                    'merscope_lung_htma':{'gene': 220, 'blank': 65},
                    'cosmx_multitissue_htma':{'gene': 1000, 'sys_control': 197, 'neg_control_probe': 10},
                    'cosmx_multitissue_normal':{'gene': 1000, 'sys_control': 197, 'neg_control_probe': 10},
                     '2024_xenium_breast_htma':{'gene': 280, 'blank': 200, 'neg_control_probe': 20, 'neg_control_codeword': 41},
                    '2024_xenium_breast_tumor2':{'gene': 280, 'blank': 200, 'neg_control_probe': 20, 'neg_control_codeword': 41},
                    '2024_merscope_breast_htma':{'gene': 255, 'blank': 30},
                    '2024_merscope_breast_tumor2':{'gene': 255, 'blank': 30},
                    '2024_cosmx_multitissue_htma':{'gene': 1000, 'sys_control': 197, 'neg_control_probe': 10},
                    '2024_cosmx_multitissue_tumor2':{'gene': 1000, 'sys_control': 197, 'neg_control_probe': 10},
                    }



xenium_breast_htma_matching_cores = [
    3, 5, 7, 8, 9, 11, 12, 15, 16, 17, 18, 26, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
    42, 43, 44, 45, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57, 58, 59, 62, 65, 66, 67, 70, 73, 74, 
    75, 78, 79, 82, 84, 88, 89, 93, 94, 95, 96, 98, 99, 103, 104, 105, 106, 109, 113, 114, 115, 
    116, 117, 118, 119, 124, 125, 126, 127, 128, 129, 132, 133, 134, 135, 136, 138, 141, 142, 143, 
    144, 145, 149, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 
    168, 169, 170
    ]


merscope_breast_htma_matching_cores = [
    43, 45, 46, 53, 54, 57, 65, 80, 81, 82, 88, 96, 101, 105, 106, 109, 
    113, 115, 116, 117, 125, 135, 136, 144, 145, 148, 156
    ]

cosmx_multitissue_htma_matching_cores = [
    11, 12, 15, 16, 17, 18, 19, 20, 21, 28, 30, 31, 32, 33, 34,
    36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51,
    52, 53, 54, 55, 56, 58, 59, 60, 61, 62, 65, 66, 67, 68, 69,
    70, 71, 73, 74, 77, 78, 79, 80, 81, 82, 83, 84, 88, 89, 90,
    91, 93, 94, 95, 96, 98, 99, 100, 101, 102, 103, 104, 105, 106, 109,
    110, 111, 112, 113, 115, 116, 117, 118, 119, 120, 121, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 141, 143, 144, 145, 146,
    147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161,
    163, 164, 165, 166, 167, 168, 169, 170, 171, 172
    ]

matching_cores = [
    31, 32, 33, 34, 37, 38, 39, 42, 43, 44, 47, 48, 49, 50, 51,
    52, 53, 54, 58, 60, 65, 66, 67, 70, 71, 72, 73, 74, 75, 76,
    78, 79, 82, 88, 89, 90, 91, 92, 94, 95, 96, 99, 105, 106, 107,
    108, 109, 112, 113, 114, 115, 116, 125, 126, 127, 128, 133,
    134, 135, 136, 143, 144, 145, 146, 147, 148, 155, 156, 157,
    205, 206, 210, 211, 215, 216, 217, 218, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235,
    236, 237, 238, 239, 240, 241, 242, 243, 244, 245,
    301, 302, 303, 304, 305, 307, 308, 309, 310, 311, 312, 313, 314,
    316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 329,
    330, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 345,
    347, 348
    ]

matching_cores_2024 = [ 
    41,  42,  43,  44,  45,  47,  48,  52,  53,  54,  56,  57,  58,
    59,  62,  64,  65,  66,  67,  68,  69,  73,  74,  75,  77,  78,
    79,  81,  82,  83,  88,  96,  97,  98,  99, 100, 101, 103, 104,
    105, 106, 107, 108, 109, 113, 114, 115, 116, 117, 118, 120, 125,
    126, 127, 128, 129, 132, 133, 134, 135, 136, 137, 138, 141, 143,
    144, 145, 146, 147, 148, 149, 152, 153, 154, 155, 156, 157, 158,
    159, 160, 163, 164, 165, 168, 169, 171,
    301, 302, 303, 304, 305, 307, 308, 309, 310, 311, 312, 313, 314,
    316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 329,
    330, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 345,
    347, 348
    ]
