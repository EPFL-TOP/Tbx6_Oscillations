# -*- coding: utf-8 -*-
"""
Created on Mon May  9 17:51:40 2022

@author: Olivier
"""

#%% Importing packages

import pandas as pd
import os
import glob
import matplotlib.pylab as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns
import va_utils
import warnings
from cycler import cycler
from scipy.stats import ttest_ind

warnings.filterwarnings("ignore")


#%% RC params 

def set_rc_params(style = 'light'):
    if style == 'light':
        background_color='1'
        grid_color = '0.85'
        font_color = '0.01'
        line_color='0.01'
        list_cycler = ['sienna',  'royalblue', 'goldenrod', 'deepskyblue', 'lime', 'hotpink', 'bisque','limegreen', ]

        
    else:
        background_color='0.15'
        grid_color = '0.35'
        font_color = 'aqua'
        line_color='aqua'
        list_cycler = ['lime', 'royalblue', 'limegreen', 'bisque', 'goldenrod', 'deepskyblue', 'hotpink', 'lime', 'teal']


          
    mpl.rcParams['lines.linewidth']=1
    mpl.rc('xtick', labelsize=5) 
    mpl.rc('ytick', labelsize=5) 
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.size'] = '7'


    mpl.rcParams['axes.facecolor']=background_color
    mpl.rcParams['figure.facecolor']=background_color
    mpl.rcParams['savefig.facecolor']=background_color
    # mpl.rcParams['figure.edgecolor']=background_color
    
    mpl.rcParams['text.color'] = font_color
    mpl.rcParams['axes.labelcolor']=font_color
    
    mpl.rcParams['axes.edgecolor']=line_color
    mpl.rcParams['xtick.color']=line_color
    mpl.rcParams['ytick.color']=line_color
    mpl.rcParams['grid.color']=line_color
    
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['grid.color'] = grid_color
    mpl.rcParams['grid.linewidth'] = '0.5'
    mpl.rcParams['grid.alpha'] = '0.25'
    
    mpl.rcParams['axes.prop_cycle'] = cycler('color', list_cycler)
    
# Seaborn style
background_color='1'
grid_color = '0.85'
font_color = '0.01'
line_color='0.01'

seaborn_style =  {
 'axes.facecolor': background_color,
 'axes.edgecolor': font_color,
 'axes.grid': False,
 'axes.axisbelow': 'line',
 'axes.labelcolor': font_color,
 'figure.facecolor': background_color,
 'grid.color': '#b0b0b0',
 'grid.linestyle': '-',
 'text.color': font_color,
 'xtick.color': font_color,
 'ytick.color': font_color,
 'xtick.direction': 'out',
 'ytick.direction': 'out',
 'patch.edgecolor': font_color,
 'patch.force_edgecolor': False,
 'image.cmap': 'viridis',
 'font.family': ['sans-serif'],
 'font.sans-serif': ['DejaVu Sans',
  'Bitstream Vera Sans',
  'Computer Modern Sans Serif',
  'Lucida Grande',
  'Verdana',
  'Geneva',
  'Lucid',
  'Arial',
  'Helvetica',
  'Avant Garde',
  'sans-serif'],
 'xtick.bottom': True,
 'xtick.top': False,
 'ytick.left': True,
 'ytick.right': False,
 'axes.spines.left': True,
 'axes.spines.bottom': True,
 'axes.spines.right': True,
 'axes.spines.top': True}
    
# color_2 = '#5C2483'
# color_1 = '#C8D300'

# color_2 = '#5C2483'
# color_1 = 'sienna'

color_1 = 'royalblue'
color_2 = 'darkorange'

color_map = 'hsv'
color_map2 = 'hot'
alpha_value = 0.75
my_color_map =['sienna',  'royalblue', 'goldenrod', 'deepskyblue', 'lime', 'hotpink', 'bisque','limegreen', ]
two_color_map = [color_1, color_2]
color_mesp_single_cells = 'grey'


linewidth_mean = 0.5
linewidth_individual = 0.1
alpha_mean = 1
alpha_individual = 0.25
alpha_std = 0.25

color_looping_mean = 'goldenrod'
color_heidi_mean = 'hotpink'
color_hulk_mean =  'goldenrod'
colors_looping = ['goldenrod',
                  'lime']
colors_heidi = ['magenta',
                  'hotpink']

color_Tbx6 = 'grey'
color_ripply1 = 'magenta'

set_rc_params(style = 'light')


     

#%% Functions


def normalize_df(df_input):
    
    df_output = pd.DataFrame({'temp': df_input.index})
    for i in df_input.columns:
        max_value = df_input[i].max(0)
        min_value = df_input[i].min(0)
        normalized_column = pd.Series(data=(df_input[i]-min_value)/(max_value-min_value), name=i + '_norm')
        df_output = df_output.join(normalized_column)
        df_output.reset_index(inplace=True, drop=True)
    df_output = df_output.drop(['temp', 'Distance_(pixels)_norm'], axis=1)
    return df_output

#%% Defining datasets

n_channels = 3
dt = 2
condition_names = [ '0.0 Posterior', '1.0 Anterior', 
                   '1.0 Posterior', '2.0 Anterior',
                   '2.0 Posterior', '3.0 Anterior', 
                   '4.0 Posterior', '5.0 Anterior', 
                   '6.0 Posterior', '7.0 Anterior', ]
condition_names_label = ['Anterior S1', 'S1 - Anterior',
                         'S1 - Posterior',  'S2 - Anterior',
                         'S2 - Posterior',  'S3 - Anterior',
                         'S4 - Posterior',  'S5 - Anterior',
                         'S6 - Posterior',  'S7 - Anterior',]

condition = 'Somite + AP'


path_figures = '../figures/'
if not os.path.exists(path_figures):
    os.makedirs(path_figures)
    

#%% Defining datasets Looping

path = '../MastodonTracks/Her1-YFP'
n_channels = 3
dt = 2
condition = 'Somite + AP'


list_cells_looping_heidi_20210901, df_looping_heidi_20210901 = va_utils.csv_to_clean_df_with_tag(path + '/Her1-YFP_20210901-edges.csv',  path + '/Her1-YFP_20210901-vertices.csv', 90, 207, n_channels, ['AP', 'Somite'], defined = False)
list_cells_looping_heidi_20210901 = va_utils.combine_categories(list_cells_looping_heidi_20210901[:], 'Somite ', 'AP ', condition)


list_cells_looping_heidi_20220727, df_looping_heidi_20220727 = va_utils.csv_to_clean_df_with_tag(path + '/Her1-YFP_20220727-edges.csv',  path + '/Her1-YFP_20220727-vertices.csv', 80, 180, n_channels, ['AP', 'Somite'], defined = False)
list_cells_looping_heidi_20220727 = va_utils.combine_categories(list_cells_looping_heidi_20220727[:], 'Somite ', 'AP ', condition)


list_cells_looping_heidi_20220805, df_looping_heidi_20220805 = va_utils.csv_to_clean_df_with_tag(path + '/Her1-YFP_20220805-edges.csv', path + '/Her1-YFP_20220805-vertices.csv', 100, 200, n_channels, ['AP', 'Somite'], defined = False)
list_cells_looping_heidi_20220805 = va_utils.combine_categories(list_cells_looping_heidi_20220805[:], 'Somite ', 'AP ', condition)

list_datasets_as_list_looping = [
    list_cells_looping_heidi_20210901,
    list_cells_looping_heidi_20220805,
    list_cells_looping_heidi_20220727

    ]

list_name_datasets_looping = [
    'looping_heidi_20210901',
    'looping_heidi_20220805',
    'looping_heidi_20220727',
    ]

list_datasets_looping_as_df = []
for ind, datasets in enumerate(list_datasets_as_list_looping):
    list_datasets_looping_as_df.append(pd.concat(list_datasets_as_list_looping[ind]))




#%% Defining datasets Hulk


path = '../MastodonTracks/Tbx6-mNG'
n_channels = 3
dt = 2
condition = 'Somite + AP'

list_cells_hulk_heidi_20211026, df_hulk_heidi_20211026 = va_utils.csv_to_clean_df_with_tag(path + '/Tbx6-mNG_20211026-edges.csv', path + '/Tbx6-mNG_20211026-vertices.csv', 80, 175, n_channels, ['AP', 'Somite'], defined = False)
list_cells_hulk_heidi_20211026 = va_utils.combine_categories(list_cells_hulk_heidi_20211026[:], 'Somite ', 'AP ', condition)

list_cells_hulk_heidi_20220315, df_hulk_heidi_20220315 = va_utils.csv_to_clean_df_with_tag(path + '/Tbx6-mNG_20220315-edges.csv', path + '/Tbx6-mNG_20220315-vertices.csv', 85, 145, n_channels, ['AP', 'Somite'], defined = False)
list_cells_hulk_heidi_20220315 = va_utils.combine_categories(list_cells_hulk_heidi_20220315[:], 'Somite ', 'AP ', condition)
 

list_cells_hulk_heidi_20220621, df_hulk_heidi_20220621 = va_utils.csv_to_clean_df_with_tag(path + '/Tbx6-mNG_20220621-edges.csv', path + '/Tbx6-mNG_20220621-vertices.csv', 85, 170, n_channels, ['AP', 'Somite'], defined = False)
list_cells_hulk_heidi_20220621 = va_utils.combine_categories(list_cells_hulk_heidi_20220621[:], 'Somite ', 'AP ', condition)


list_datasets_as_list_hulk = [
    list_cells_hulk_heidi_20211026,
    list_cells_hulk_heidi_20220315,
    list_cells_hulk_heidi_20220621,
    ]

list_name_datasets_hulk = [
    'hulk-heidi 20211026',
    'hulk_heidi_20220315',
    'hulk_heidi_20220621',
    ]


list_datasets_hulk_as_df = []
for ind, datasets in enumerate(list_datasets_as_list_hulk):
    list_datasets_hulk_as_df.append(pd.concat(list_datasets_as_list_hulk[ind]))
    
    
    

list_datasets_as_list = [
    list_cells_hulk_heidi_20211026,
    list_cells_hulk_heidi_20220315,
    list_cells_hulk_heidi_20220621,
    list_cells_looping_heidi_20210901,
    list_cells_looping_heidi_20220805,
    list_cells_looping_heidi_20220727
    ]


list_name_datasets = [
    'hulk-heidi 20211026',
    'hulk_heidi_20220315',
    'hulk_heidi_20220621',
    'looping_heidi_20210901',
    'looping_heidi_20220805',
    'looping_heidi_20220727',
    ]


list_datasets_as_df = []
for ind, datasets in enumerate(list_datasets_as_list):
    list_datasets_as_df.append(pd.concat(list_datasets_as_list[ind]))
    list_datasets_as_df[ind] =  list_datasets_as_df[ind].reset_index()
    


#%% Defining datasets Hulk midSom


path = '../MastodonTracks/Tbx6-mNG_midSom'
n_channels = 2
dt = 2
condition = 'Somite + AP'

list_cells_hulk_midSom_20230119_pos1, df_hulk_midSom_20230119_pos1 = va_utils.csv_to_clean_df_with_tag(path + '/Tbx6-mNG_midSom_20230119_pos1-edges.csv', path + '/Tbx6-mNG_midSom_20230119_pos1-vertices.csv', 20, 120, n_channels, ['AP', 'Somite'], defined = False)
list_cells_hulk_midSom_20230119_pos1 = va_utils.combine_categories(list_cells_hulk_midSom_20230119_pos1[:], 'Somite ', 'AP ', condition)

list_cells_hulk_midSom_20230119_pos2, df_hulk_midSom_20230119_pos2 = va_utils.csv_to_clean_df_with_tag(path + '/Tbx6-mNG_midSom_20230119_pos2-edges.csv', path + '/Tbx6-mNG_midSom_20230119_pos2-vertices.csv', 20, 120, n_channels, ['AP', 'Somite'], defined = False)
list_cells_hulk_midSom_20230119_pos2 = va_utils.combine_categories(list_cells_hulk_midSom_20230119_pos2[:], 'Somite ', 'AP ', condition)



list_datasets_as_list_hulk_midSom = [
    list_cells_hulk_midSom_20230119_pos1,
    list_cells_hulk_midSom_20230119_pos2,
    ]

list_name_datasets_hulk_midSom = [
    'hulk_midSom_20230119_pos1',
    'hulk_midSom_20230119_pos2',

    ]


list_datasets_hulk_midSom_as_df = []
for ind, datasets in enumerate(list_datasets_as_list_hulk_midSom):
    list_datasets_hulk_midSom_as_df.append(pd.concat(list_datasets_as_list_hulk_midSom[ind]))
    list_datasets_hulk_midSom_as_df[ind] = list_datasets_hulk_midSom_as_df[ind].reset_index()



#%% Defining datasets Hulk her1;her7 mutants


path = '../MastodonTracks/Tbx6-mNG_her1_her7_mutants'
n_channels = 2
dt = 2
condition = 'Somite + AP'

list_hulk_gullum_20220902_pos1, df_hulk_gullum_20220902_pos1 = va_utils.csv_to_clean_df_with_tag(path + '/Tbx6-mNG_her1_her7_mutants_20220902_pos1-edges.csv', path +  '/Tbx6-mNG_her1_her7_mutants_20220902_pos1-vertices.csv', 45, 120, n_channels, ['AP', 'Somite'],defined = False)
list_hulk_gullum_20220902_pos1 = va_utils.combine_categories(list_hulk_gullum_20220902_pos1[:], 'Somite ', 'AP ', condition)

list_hulk_gullum_20220902_pos2, df_hulk_gullum_20220902_pos2 = va_utils.csv_to_clean_df_with_tag(path + '/Tbx6-mNG_her1_her7_mutants_20220902_pos2-edges.csv',path +   '/Tbx6-mNG_her1_her7_mutants_20220902_pos2-vertices.csv', 60, 135, n_channels, ['AP', 'Somite'],defined = False)
list_hulk_gullum_20220902_pos2 = va_utils.combine_categories(list_hulk_gullum_20220902_pos2[:], 'Somite ', 'AP ', condition)


list_datasets_hulk_gullum_as_list = [
    list_hulk_gullum_20220902_pos1,
    list_hulk_gullum_20220902_pos2,


    ]

list_name_datasets_hulk_gullum = [
    'hulk_gullum_20220902_pos1',
    'hulk_gullum_20220902_pos2',

    ]

list_datasets_hulk_gullum_as_df = []
for ind, datasets in enumerate(list_datasets_hulk_gullum_as_list):
    list_datasets_hulk_gullum_as_df.append(pd.concat(list_datasets_hulk_gullum_as_list[ind]))
    list_datasets_hulk_gullum_as_df[ind] = list_datasets_hulk_gullum_as_df[ind].reset_index()



    


#%% Aligning signals

# 20211026: epiboly completed at TP 88
# 20220315: epiboly completed at TP 76
# 20220621: epiboly completed at TP 102

# 20210901: epiboly completed at TP 147
# 20220805: epiboly completed at TP 132
# 20220727: epiboly completed at TP 131


# Shift between 20210901 and 20220805: 15

# Shift between 20220621 and 20210901: 41

# Signals aligned with respect to the latest one (20210901).




shift = [59 , 71, 45, 0, 15, 16]


for i, df in enumerate(list_datasets_as_df):

    if 'Spot frame adjusted' in list_datasets_as_df[i].columns: list_datasets_as_df[i].drop(columns = ['Spot frame adjusted'], inplace = True)
    frame_adjusted = pd.DataFrame({'Spot frame adjusted': list_datasets_as_df[i]['Spot frame'] + shift[i]})
    list_datasets_as_df[i].reset_index(drop = True)
    list_datasets_as_df[i] = df.join(frame_adjusted)


    
    

     
#%% Pooled Anterior vs posterior Her1-YFP - 1st to 3rd boundary 1 plot

condition_loop = [
    [ '0.0 Posterior', '1.0 Anterior'],
    [ '1.0 Posterior', '2.0 Anterior'],
    [ '2.0 Posterior', '3.0 Anterior']
    ]


fig, ax = plt.subplots(3,1, figsize=(2, 3),
                         sharex = True, sharey = False)
axs = ax.ravel()
df = pd.concat(list_datasets_as_df[3:])
for ind_boundary, boundary in enumerate(condition_loop):
    
    condition_names = boundary
    
    for ind2, cond in enumerate(condition_names):
        
        count = 0
        for ind_df, df_looping in enumerate(list_datasets_as_df[3:]):
            for j in df_looping.loc[df_looping[condition] == cond]['Cell Idx'].astype(int).drop_duplicates():
                axs[ind_boundary].plot(df_looping.loc[df_looping['Cell Idx'] == j]['Spot frame adjusted']*dt, 
                        df_looping.loc[df_looping['Cell Idx'] == j]['Spot center intensity'],
                        alpha =alpha_individual,
                        color = two_color_map[ind2],
                        linewidth = linewidth_individual)
                count = count +1
        print(cond + ' :' )
        print(count) 
           
    for ind2, cond in enumerate(condition_names):
  
        axs[ind_boundary].fill_between(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                               df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                               df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                               - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                           df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values - df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity'].values, 
                           df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values + df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity'].values,
                           alpha = alpha_std,
                           linewidth = 0, 
                           color =  two_color_map[ind2])
        

            
    for ind2, cond in enumerate(condition_names):
        axs[ind_boundary].plot(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                            df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                            df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                            - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values,
                color =  two_color_map[ind2], 
                linewidth = linewidth_mean,
                alpha = alpha_mean
                )
      
    axs[ind_boundary].set_xlim(left = 0, right =400)
    axs[ind_boundary].spines['bottom'].set_position(('outward', 5))
    axs[ind_boundary].spines['left'].set_position(('outward', 5))
  
    axs[ind_boundary].set_ylabel('Intensity', fontsize = 7)
    
axs[-1].set_xlabel('Time [min]', fontsize = 7)
plt.tight_layout()
fig.savefig(path_figures+ 'Fig1d.pdf', dpi = 450)   



    
    


#%% Pooled Anterior vs posterior Tbx6-mNG - 1st to 3rd boundary 1 plot

condition_loop = [
    [ '0.0 Posterior', '1.0 Anterior'],
    [ '1.0 Posterior', '2.0 Anterior'],
    [ '2.0 Posterior', '3.0 Anterior']
    ]

    
fig, ax = plt.subplots(3,1, figsize=(2,3),
                    sharex = True, sharey = False)
axs = ax.ravel()
df = pd.concat(list_datasets_as_df[:3])
for ind_boundary, boundary in enumerate(condition_loop):

    condition_names = boundary
    
        
    for ind2, cond in enumerate(condition_names):
        
        count = 0
        for ind_df, df_hulk in enumerate(list_datasets_as_df[:3]):
            for j in df_hulk.loc[df_hulk[condition] == cond]['Cell Idx'].astype(int).drop_duplicates():
                axs[ind_boundary].plot(df_hulk.loc[df_hulk['Cell Idx'] == j]['Spot frame adjusted']*dt, 
                        df_hulk.loc[df_hulk['Cell Idx'] == j]['Spot center intensity'],
                        alpha =alpha_individual,
                        color = two_color_map[ind2],
                        linewidth = linewidth_individual)
                count = count + 1
        print(cond + ' :' )
        print(count)

    for ind2, cond in enumerate(condition_names):
  
        axs[ind_boundary].fill_between(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                               df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                               df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                               - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                           df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values - df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity'].values, 
                           df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values + df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity'].values,
                           alpha = alpha_std,
                           linewidth = 0, 
                           color =  two_color_map[ind2])

          
    for ind2, cond in enumerate(condition_names):
            
        axs[ind_boundary].plot(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                            df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                            df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                            - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values,
                color =  two_color_map[ind2], 
                linewidth = linewidth_mean,
                alpha = alpha_mean
                )
      
    axs[ind_boundary].set_xlim(left = 100, right =375)
    axs[ind_boundary].spines['bottom'].set_position(('outward', 5))
    axs[ind_boundary].spines['left'].set_position(('outward', 5))
  
   
    axs[ind_boundary].set_ylabel('Intensity', fontsize = 7)
axs[-1].set_xlabel('Time [min]', fontsize = 7)
plt.tight_layout()
fig.savefig(path_figures+ 'Fig2b.pdf', dpi = 450)   





#%% Pooled Anterior vs posterior hulk - B5 

condition_loop = [
    [ '4.0 Posterior', '5.0 Anterior'],
    ]

name_loop = ['Ant_vs_post_hulk_boundary5_pooled.pdf',
             ]

    
fig, ax = plt.subplots(1,1, figsize=(1.75,1.25),
                    sharex = True, sharey = False)
df = pd.concat(list_datasets_as_df[:3])


for ind_boundary, boundary in enumerate(condition_loop):

    condition_names = boundary
        
    for ind2, cond in enumerate(condition_names):
        
       
        count = 0
        for ind_df, df_hulk in enumerate(list_datasets_as_df[:3]):
        
            for j in df_hulk.loc[df_hulk[condition] == cond]['Cell Idx'].astype(int).drop_duplicates():
                ax.plot(df_hulk.loc[df_hulk['Cell Idx'] == j]['Spot frame adjusted']*dt, 
                        df_hulk.loc[df_hulk['Cell Idx'] == j]['Spot center intensity'],
                        alpha =alpha_individual,
                        color = two_color_map[ind2],
                   linewidth = linewidth_individual)
                count = count +1
        print(cond)
        print(count)
        
    for ind2, cond in enumerate(condition_names):
  
        ax.fill_between(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                                df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                                df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                                - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                            df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values - df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity'].values, 
                            df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values + df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity'].values,
                            alpha = alpha_std,
                            linewidth = 0, 
                            color =  two_color_map[ind2])

          
    for ind2, cond in enumerate(condition_names):
            
       ax.plot(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                            df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                            df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                            - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values,
                color =  two_color_map[ind2], 
                linewidth = linewidth_mean,
                alpha = alpha_mean
                )
      
  
ax.set_xlim(left = 215, right =400)
ax.spines['bottom'].set_position(('outward', 5))
ax.spines['left'].set_position(('outward', 5))
  
ax.set_ylabel('Intensity', fontsize = 7)
ax.set_xlabel('Time [min]', fontsize = 7)
plt.tight_layout()
fig.savefig(path_figures+ 'Fig2b3prime.pdf', dpi = 450)   




  
#%% 20220621 Anterior vs posterior hulk - B7

condition_loop = [
    [ '6.0 Posterior', '7.0 Anterior'],
    ]

fig, ax = plt.subplots(1,1, figsize=(1.75,1.25),
                    sharex = True, sharey = False)
df = list_datasets_as_df[1]


for ind_boundary, boundary in enumerate(condition_loop):

    condition_names = boundary
        
    for ind2, cond in enumerate(condition_names):
        
        count = 0
        df_hulk = df
        for j in df_hulk.loc[df_hulk[condition] == cond]['Cell Idx'].astype(int).drop_duplicates():
            ax.plot(df_hulk.loc[df_hulk['Cell Idx'] == j]['Spot frame adjusted']*dt, 
                    df_hulk.loc[df_hulk['Cell Idx'] == j]['Spot center intensity'],
                    alpha =alpha_individual,
                    color = two_color_map[ind2],
                    linewidth = linewidth_individual)
            count = count +1
        print(cond + ' :' )
        print(count)
        
    for ind2, cond in enumerate(condition_names):
  
        ax.fill_between(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                                df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                                df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                                - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                            df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values - df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity'].values, 
                            df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values + df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity'].values,
                            alpha = alpha_std,
                            linewidth = 0, 
                            color =  two_color_map[ind2])

          
    for ind2, cond in enumerate(condition_names):
            
       ax.plot(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                            df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                            df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                            - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity'].values,
                color =  two_color_map[ind2], 
                linewidth = linewidth_mean,
                alpha = alpha_mean
                )
          
ax.set_xlim(left = 220, right =450)
ax.spines['bottom'].set_position(('outward', 5))
ax.spines['left'].set_position(('outward', 5))
  
ax.set_ylabel('Intensity', fontsize = 7)
ax.set_xlabel('Time [min]', fontsize = 7)
plt.tight_layout()
fig.savefig(path_figures+ 'Fig2c.pdf', dpi = 450)   

  
#%% 20220621 Anterior vs posterior hulk - B7 (Extended Data)

condition_loop = [
    [ '6.0 Posterior', '7.0 Anterior'],
    ]

name_loop = ['Ant_vs_post_hulk_boundary7_pooled.pdf',
             ]

fig, ax = plt.subplots(1,1, figsize=(1.75,1.25),
                    sharex = True, sharey = False)
df = list_datasets_as_df[2]


for ind_boundary, boundary in enumerate(condition_loop):

    condition_names = boundary
        
    for ind2, cond in enumerate(condition_names):
        
        count = 0
        df_hulk = df
        for j in df_hulk.loc[df_hulk[condition] == cond]['Cell Idx'].astype(int).drop_duplicates():
            ax.plot(df_hulk.loc[df_hulk['Cell Idx'] == j]['Spot frame']*dt, 
                    df_hulk.loc[df_hulk['Cell Idx'] == j]['Spot center intensity'],
                    alpha =alpha_individual,
                    color = two_color_map[ind2],
                    linewidth = linewidth_individual)
            count = count +1
        print(cond + ' :' )
        print(count)
        
    for ind2, cond in enumerate(condition_names):
  
        ax.fill_between(np.linspace(df.loc[df[condition]==cond]['Spot frame'].min(axis=0),
                                df.loc[df[condition]==cond]['Spot frame'].max(axis=0), 
                                df.loc[df[condition]==cond]['Spot frame'].max(axis=0)
                                - df.loc[df[condition]==cond]['Spot frame'].min(axis=0) +1)*dt,
                            df.loc[df[condition] == cond].groupby('Spot frame').mean()['Spot center intensity'].values - df.loc[df[condition] == cond].groupby('Spot frame').std()['Spot center intensity'].values, 
                            df.loc[df[condition] == cond].groupby('Spot frame').mean()['Spot center intensity'].values + df.loc[df[condition] == cond].groupby('Spot frame').std()['Spot center intensity'].values,
                            alpha = alpha_std,
                            linewidth = 0, 
                            color =  two_color_map[ind2])

          
    for ind2, cond in enumerate(condition_names):
            
        ax.plot(np.linspace(df.loc[df[condition]==cond]['Spot frame'].min(axis=0),
                            df.loc[df[condition]==cond]['Spot frame'].max(axis=0), 
                            df.loc[df[condition]==cond]['Spot frame'].max(axis=0)
                            - df.loc[df[condition]==cond]['Spot frame'].min(axis=0) +1)*dt,
                df.loc[df[condition] == cond].groupby('Spot frame').mean()['Spot center intensity'].values,
                color =  two_color_map[ind2], 
                linewidth = linewidth_mean,
                alpha = alpha_mean
                )
          
# ax.set_xlim(left = 220, right =450)
ax.spines['bottom'].set_position(('outward', 5))
ax.spines['left'].set_position(('outward', 5))
  
ax.set_ylabel('Intensity', fontsize = 7)
ax.set_xlabel('Time [min]', fontsize = 7)
plt.tight_layout()
fig.savefig(path_figures+ 'ExtendedFig3a.pdf', dpi = 450)   


  
  
    
#%% Boundary 16


condition_loop = [
    [ '16.0 Posterior', '17.0 Anterior'],
    ]


fig, ax = plt.subplots(1,1, figsize=(1.75,1.25),
                    sharex = True, sharey = False)

index_dataset = 1
df = list_datasets_hulk_midSom_as_df[index_dataset]


for ind_boundary, boundary in enumerate(condition_loop):

    condition_names = boundary
        
    for ind2, cond in enumerate(condition_names):
        
        count = 0
        for j in df.loc[df[condition] == cond]['Cell Idx'].astype(int).drop_duplicates():
            ax.plot(df.loc[df['Cell Idx'] == j]['Spot frame']*dt, 
                    df.loc[df['Cell Idx'] == j]['Spot center intensity'],
                    alpha =alpha_individual,
                    color = two_color_map[ind2],
                    linewidth = linewidth_individual)
            count = count +1
        print(cond + ' :' )
        print(count)
    
            
    for ind2, cond in enumerate(condition_names):
  
        ax.fill_between(np.linspace(df.loc[df[condition]==cond]['Spot frame'].min(axis=0),
                                df.loc[df[condition]==cond]['Spot frame'].max(axis=0), 
                                df.loc[df[condition]==cond]['Spot frame'].max(axis=0)
                                - df.loc[df[condition]==cond]['Spot frame'].min(axis=0) +1)*dt,
                            df.loc[df[condition] == cond].groupby('Spot frame').mean()['Spot center intensity'].values - df.loc[df[condition] == cond].groupby('Spot frame').std()['Spot center intensity'].values, 
                            df.loc[df[condition] == cond].groupby('Spot frame').mean()['Spot center intensity'].values + df.loc[df[condition] == cond].groupby('Spot frame').std()['Spot center intensity'].values,
                            alpha = alpha_std,
                            linewidth = 0, 
                            color =  two_color_map[ind2])

          
    for ind2, cond in enumerate(condition_names):
            
        ax.plot(np.linspace(df.loc[df[condition]==cond]['Spot frame'].min(axis=0),
                            df.loc[df[condition]==cond]['Spot frame'].max(axis=0), 
                            df.loc[df[condition]==cond]['Spot frame'].max(axis=0)
                            - df.loc[df[condition]==cond]['Spot frame'].min(axis=0) +1)*dt,
                df.loc[df[condition] == cond].groupby('Spot frame').mean()['Spot center intensity'].values,
                color =  two_color_map[ind2], 
                linewidth = linewidth_mean,
                alpha = alpha_mean
                )
       
# ax.set_xlim(left = 0, right =120)
ax.spines['bottom'].set_position(('outward', 5))
ax.spines['left'].set_position(('outward', 5))
  
ax.set_ylabel('Intensity', fontsize = 7)
ax.set_xlabel('Time [min]', fontsize = 7)
plt.tight_layout()
fig.savefig(path_figures+ 'Fig2d.pdf', dpi = 450)   
  

    
#%% Boundary 16



condition_loop = [
    [ '16.0 Posterior', '17.0 Anterior'],
    ]


fig, ax = plt.subplots(1,1, figsize=(1.75,1.25),
                    sharex = True, sharey = False)

index_dataset = 0
df = list_datasets_hulk_midSom_as_df[index_dataset]


for ind_boundary, boundary in enumerate(condition_loop):

    condition_names = boundary
        
    for ind2, cond in enumerate(condition_names):
        
        count = 0
        for j in df.loc[df[condition] == cond]['Cell Idx'].astype(int).drop_duplicates():
            ax.plot(df.loc[df['Cell Idx'] == j]['Spot frame']*dt, 
                    df.loc[df['Cell Idx'] == j]['Spot center intensity'],
                    alpha =alpha_individual,
                    color = two_color_map[ind2],
                    linewidth = linewidth_individual)
            count = count +1
        print(cond + ' :' )
        print(count)
    
            
    for ind2, cond in enumerate(condition_names):
  
        ax.fill_between(np.linspace(df.loc[df[condition]==cond]['Spot frame'].min(axis=0),
                                df.loc[df[condition]==cond]['Spot frame'].max(axis=0), 
                                df.loc[df[condition]==cond]['Spot frame'].max(axis=0)
                                - df.loc[df[condition]==cond]['Spot frame'].min(axis=0) +1)*dt,
                            df.loc[df[condition] == cond].groupby('Spot frame').mean()['Spot center intensity'].values - df.loc[df[condition] == cond].groupby('Spot frame').std()['Spot center intensity'].values, 
                            df.loc[df[condition] == cond].groupby('Spot frame').mean()['Spot center intensity'].values + df.loc[df[condition] == cond].groupby('Spot frame').std()['Spot center intensity'].values,
                            alpha = alpha_std,
                            linewidth = 0, 
                            color =  two_color_map[ind2])

          
    for ind2, cond in enumerate(condition_names):
            
        ax.plot(np.linspace(df.loc[df[condition]==cond]['Spot frame'].min(axis=0),
                            df.loc[df[condition]==cond]['Spot frame'].max(axis=0), 
                            df.loc[df[condition]==cond]['Spot frame'].max(axis=0)
                            - df.loc[df[condition]==cond]['Spot frame'].min(axis=0) +1)*dt,
                df.loc[df[condition] == cond].groupby('Spot frame').mean()['Spot center intensity'].values,
                color =  two_color_map[ind2], 
                linewidth = linewidth_mean,
                alpha = alpha_mean
                )
       
# ax.set_xlim(left = 0, right =120)
ax.spines['bottom'].set_position(('outward', 5))
ax.spines['left'].set_position(('outward', 5))
  
ax.set_ylabel('Intensity', fontsize = 7)
ax.set_xlabel('Time [min]', fontsize = 7)
plt.tight_layout()
fig.savefig(path_figures+ 'ExtendedFig3b.pdf', dpi = 450)   




#%% Hulk traces in her1;her7 mutants

condition_loop = [
    '1.0 Anterior',
    '2.0 Anterior',
    '3.0 Anterior'
    ]

index_dataset = 0
df = list_datasets_hulk_gullum_as_df[index_dataset]


fig, ax = plt.subplots(3,1, figsize=(2,3),
                    sharex = True, sharey = False)
axs = ax.ravel()

for ind_boundary, boundary in enumerate(condition_loop):

   
    
    count = 0
   
    for j in df.loc[df[condition] == boundary]['Cell Idx'].astype(int).drop_duplicates():
        axs[ind_boundary].plot(df.loc[df['Cell Idx'] == j]['Spot frame']*dt, 
                df.loc[df['Cell Idx'] == j]['Spot center intensity'],
                alpha =alpha_individual*4,
            
                linewidth = linewidth_individual*4)
        count = count + 1
    print(boundary + ' :' )
    print(count)
    
    axs[ind_boundary].spines['bottom'].set_position(('outward', 5))
    axs[ind_boundary].spines['left'].set_position(('outward', 5))
      
   
    axs[ind_boundary].set_ylabel('Intensity', fontsize = 7)
axs[-1].set_xlabel('Time [min]', fontsize = 7)
plt.tight_layout()
fig.savefig(path_figures+ 'Fig2f.pdf', dpi = 450)   



#%% Hulk traces in her1;her7 mutants

     
index_dataset = 1
df = list_datasets_hulk_gullum_as_df[index_dataset]


fig, ax = plt.subplots(3,1, figsize=(2,3),
                    sharex = True, sharey = False)
axs = ax.ravel()

for ind_boundary, boundary in enumerate(condition_loop):

   
    
    count = 0
   
    for j in df.loc[df[condition] == boundary]['Cell Idx'].astype(int).drop_duplicates():
        axs[ind_boundary].plot(df.loc[df['Cell Idx'] == j]['Spot frame']*dt, 
                df.loc[df['Cell Idx'] == j]['Spot center intensity'],
                alpha =alpha_individual*4,
            
                linewidth = linewidth_individual*4)
        count = count + 1
    print(boundary + ' :' )
    print(count)
    
    axs[ind_boundary].spines['bottom'].set_position(('outward', 5))
    axs[ind_boundary].spines['left'].set_position(('outward', 5))
      
   
    axs[ind_boundary].set_ylabel('Intensity', fontsize = 7)
axs[-1].set_xlabel('Time [min]', fontsize = 7)
plt.tight_layout()
fig.savefig(path_figures+ 'ExtendedFig3c.pdf', dpi = 450)   

        

#%% Defining datasets snapshots B3

    
path = '../ROI_profiles/Snapshots_timelapse20220621_B3/'
list_df = []
list_df_norm = []
list_filenames = []
norm_factor = 2


for filename in glob.glob(os.path.join(path, '*.csv')):

    with open(os.path.join(os.getcwd(), filename), 'r') as f: 

        df = pd.read_csv(filename)
        df = df.rename(columns={"Gray_Value" : "Tbx6"})
        list_df.append(df)
        df_norm = normalize_df(df)
        list_df_norm.append(df_norm)
        list_filenames.append(filename[len(path):-4])
        

for ind, df in enumerate(list_df):
        
    fig, ax = plt.subplots(1,1, figsize=(1.5,1))
    ax.plot(df['Distance_(pixels)']*norm_factor, df['Tbx6'], color = color_Tbx6)
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['left'].set_position(('outward', 5))
    ax.set_ylabel('Intensity', fontsize = 7)
    ax.set_xlabel('Posterior to anterior', fontsize = 7)
    ax.set_xlim(left = 0, right =250)
    plt.tight_layout()
    fig.savefig(path_figures + list_filenames[ind] + '_ExtendedFig4.pdf', dpi = 450)   

#%% Defining datasets snapshots B5


path = '../ROI_profiles/Snapshots_timelapse20220621_B5/'
list_df = []
list_df_norm = []
list_filenames = []
norm_factor = 2


for filename in glob.glob(os.path.join(path, '*.csv')):

    with open(os.path.join(os.getcwd(), filename), 'r') as f:

        df = pd.read_csv(filename)
        df = df.rename(columns={"Value" : "Tbx6"})
        list_df.append(df)
        df_norm = normalize_df(df)
        list_df_norm.append(df_norm)
        list_filenames.append(filename[len(path):-4])
        

for ind, df in enumerate(list_df):
        
    fig, ax = plt.subplots(1,1, figsize=(1.5,1))
    ax.plot(df['Distance_(pixels)']*norm_factor, df['Tbx6'], color = color_Tbx6)
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['left'].set_position(('outward', 5))
    ax.set_ylabel('Intensity', fontsize = 7)
    ax.set_xlabel('Posterior to anterior', fontsize = 7)
    ax.set_xlim(left = 0, right =250)
    plt.tight_layout()
    fig.savefig(path_figures + list_filenames[ind] + '_Fig3ace.pdf', dpi = 450)   
    


#%% Defining datasets Tbx6 + ripply1


path = '../ROI_profiles/Tbx6_ripply1/'
list_df = []
list_df_norm = []
list_filenames = []



for filename in glob.glob(os.path.join(path, '*.csv')):

    with open(os.path.join(os.getcwd(), filename), 'r') as f:

        df = pd.read_csv(filename)
        # df = df.rename(columns={"Value" : "Tbx6"})
        list_df.append(df)
        df_norm = normalize_df(df)
        list_df_norm.append(df_norm)
        list_filenames.append(filename[len(path):-4])
        

for ind, df in enumerate(list_df_norm):
        
    fig, ax = plt.subplots(1,1, figsize=(1.5,1))
    ax.plot(df['Tbx6_norm'], color = color_Tbx6)
    ax.plot(df['ripply1_norm'], color = color_ripply1)
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['left'].set_position(('outward', 5))
    ax.set_ylabel('Intensity', fontsize = 7)
    ax.set_xlabel('Posterior to anterior', fontsize = 7)
    ax.set_xlim(left = 0, right =250)
    plt.tight_layout()
    fig.savefig(path_figures + list_filenames[ind] + '_Fig3bdf.pdf', dpi = 450)   
    

  
#%% Pooled Anterior vs posterior Tbx6-mNG (mesp) - 1st to 3rd boundary 1 plot (Extended Data)

condition_loop = [
    [ '0.0 Posterior', '1.0 Anterior'],
    [ '1.0 Posterior', '2.0 Anterior'],
    [ '2.0 Posterior', '3.0 Anterior']
    ]

last_tp_df = [468, 432, 436]
list_df_data_temp = []
    
for ind_df, df in enumerate(list_datasets_as_df[:3]):
    
    print(list_name_datasets[ind_df])
    fig, ax = plt.subplots(3,1, figsize=(2,3),
                        sharex = True, sharey = False)
    axs = ax.ravel()
    for ind_boundary, boundary in enumerate(condition_loop):
    
        condition_names = boundary
        
        for ind2, cond in enumerate(condition_names):
            
            count = 0
           
            for j in df.loc[df[condition] == cond]['Cell Idx'].astype(int).drop_duplicates():
                axs[ind_boundary].plot(df.loc[df['Cell Idx'] == j]['Spot frame adjusted']*dt, 
                        df.loc[df['Cell Idx'] == j]['Spot center intensity.1'],
                        alpha =alpha_individual,
                        color = two_color_map[ind2],
                        linewidth = linewidth_individual)
                count = count + 1
            print(cond + ' :' )
            print(count)
    
        for ind2, cond in enumerate(condition_names):
      
            axs[ind_boundary].fill_between(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                                   df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                                   df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                                   - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                               df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity.1'].values - df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity.1'].values, 
                               df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity.1'].values + df.loc[df[condition] == cond].groupby('Spot frame adjusted').std()['Spot center intensity.1'].values,
                               alpha = alpha_std,
                               linewidth = 0, 
                               color =  two_color_map[ind2])
    
              
        for ind2, cond in enumerate(condition_names):
                
            axs[ind_boundary].plot(np.linspace(df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0),
                                df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0), 
                                df.loc[df[condition]==cond]['Spot frame adjusted'].max(axis=0)
                                - df.loc[df[condition]==cond]['Spot frame adjusted'].min(axis=0) +1)*dt,
                    df.loc[df[condition] == cond].groupby('Spot frame adjusted').mean()['Spot center intensity.1'].values,
                    color =  two_color_map[ind2], 
                    linewidth = linewidth_mean,
                    alpha = alpha_mean
                    )
          
        # axs[ind_boundary].set_xlim(left = 100, right =375)
        axs[ind_boundary].spines['bottom'].set_position(('outward', 5))
        axs[ind_boundary].spines['left'].set_position(('outward', 5))
      
       
        axs[ind_boundary].set_ylabel('Intensity', fontsize = 7)
    axs[-1].set_xlabel('Time [min]', fontsize = 7)
    plt.tight_layout()
    fig.savefig(path_figures+ 'ExtendedFig2_'+ list_name_datasets[ind_df] + '.pdf', dpi = 450)   
    
    
    
    df_last_tp = df.loc[df['Spot frame adjusted'] == last_tp_df[ind_df]/dt]
    
    
    print(ttest_ind(df_last_tp.loc[(df_last_tp[condition]== '0.0 Posterior')]['Spot center intensity.1'], 
                df_last_tp.loc[(df_last_tp[condition]== '1.0 Anterior')]['Spot center intensity.1']))
    
    
    print(ttest_ind(df_last_tp.loc[(df_last_tp[condition]== '1.0 Posterior')]['Spot center intensity.1'], 
                df_last_tp.loc[(df_last_tp[condition]== '2.0 Anterior')]['Spot center intensity.1']))
    
     
    
    print(ttest_ind(df_last_tp.loc[(df_last_tp[condition]== '2.0 Posterior')]['Spot center intensity.1'], 
                df_last_tp.loc[(df_last_tp[condition]== '3.0 Anterior')]['Spot center intensity.1']))


    data_temp = {'Spot center intensity.1': df_last_tp['Spot center intensity.1'],
                  condition:  df_last_tp[condition],
                'dataset': list_name_datasets[ind_df]
                              }

    df_data_temp = pd.DataFrame(data = data_temp)
    list_df_data_temp.append(df_data_temp)
                                                                                                                                                                                                                                                                   

data = pd.concat(list_df_data_temp)
sns.set_style("white", rc=seaborn_style)

p = sns.catplot(data= data.loc[(data[condition]== '0.0 Posterior') | (data[condition]=='1.0 Anterior')], x="dataset", y="Spot center intensity.1", hue=condition, palette=two_color_map, kind="box",linewidth =0.25,)
p.despine(offset = 5, trim=True)
p.set_axis_labels("Datasets", "Intensity", labelpad=10, fontsize=7)
p.legend.set_title("Cell population")
p.figure.set_size_inches(2, 1.5)
p.set_xticklabels(['1','2','3'])
sns.move_legend(p, "upper left", bbox_to_anchor=(1, 1))
p.tight_layout()
p.savefig(path_figures+ 'ExtendedFig2c.pdf', dpi=450)

p = sns.catplot(data= data.loc[(data[condition]== '1.0 Posterior') | (data[condition]=='2.0 Anterior')], x="dataset", y="Spot center intensity.1", hue=condition, palette=two_color_map, kind="box",linewidth =0.25 )
p.despine(offset = 5, trim=True)
p.set_axis_labels("Datasets", "Intensity", labelpad=10, fontsize=7)
p.legend.set_title("Cell population")
p.figure.set_size_inches(2, 1.5)
p.set_xticklabels(['1','2','3'])
sns.move_legend(p, "upper left", bbox_to_anchor=(1, 1))
p.tight_layout()
p.savefig(path_figures+ 'ExtendedFig2cprime.pdf', dpi=450)


p = sns.catplot(data= data.loc[(data[condition]== '2.0 Posterior') | (data[condition]=='3.0 Anterior')], x="dataset", y="Spot center intensity.1", hue=condition, palette=reversed(two_color_map), kind="box", linewidth =0.25)
p.despine(offset = 5, trim=True)
p.set_axis_labels("Datasets", "Intensity", labelpad=10, fontsize=7)
p.legend.set_title("Cell population")
p.figure.set_size_inches(2, 1.5)
p.set_xticklabels(['1','2','3'])
sns.move_legend(p, "upper left", bbox_to_anchor=(1, 1))
p.tight_layout()
p.savefig(path_figures+ 'ExtendedFig2c2prime.pdf', dpi=450)

  
    

