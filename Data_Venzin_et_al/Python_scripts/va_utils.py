9# -*- coding: utf-8 -*-
"""
Created on Thu May  5 11:32:06 2022

@author: Olivier
"""
"""
Utilities for vortex analysis
"""


#%% Importing packaging 
import pandas as pd
import numpy as np

#%% csv to df


    
    
def csv_to_df(filename):
    """Getting DataFrame out of a csv file"""
    df = pd.read_csv(filename)
    df = pd.DataFrame(df)
    df = df.drop([df.index[0], df.index[1]])
    df = df.drop(['Label'], axis=1)
    df.to_numpy().astype(float)
    return df


def check_df_e_ok(df_e):
    if (df_e['Link target IDs'] == df_e['Link target IDs.1']).any():
        print('Error with CSV file: Link target IDs = Link target IDs.1')
        return False
    else:
        return True


def df_to_mother_cells_tracks(df_e, df_v):
    """Getting mother cells tracks from the Dataframes"""
    list_mother_cells_v = []
    list_mother_cells_e = []
    
    df_e['ID'] = df_e['ID'].to_numpy().astype(int)
    df_v['ID'] = df_v['ID'].to_numpy().astype(int)


    for i in range(df_v['Spot track ID'].to_numpy().astype(int).max() + 1):
        df_temp_v = df_v[df_v['Spot track ID'].to_numpy().astype(int) == i];
        df_temp_v['Spot frame'] = df_temp_v['Spot frame'].apply(lambda x: x if pd.isnull(x) else int(int(x)))
        df_temp_v = df_temp_v.sort_values(by='Spot frame', ignore_index=True);                                # Sorting by Spot frame
        df_temp_e = df_e[df_e['Link target IDs'].isin(df_temp_v['ID'].to_numpy().astype(int))]
        list_mother_cells_v.append(df_temp_v)
        list_mother_cells_e.append(df_temp_e)   
    
    return list_mother_cells_e, list_mother_cells_v


def is_source_id(df_e, df_v):
    
    """Check source is in link IDs or IDs.1"""
    
    min_frame = df_v['Spot frame'].loc[df_v['Spot track ID'].astype(int)==0].min()
    first_id = df_v['ID'].loc[df_v['Spot track ID'].astype(int)==0].loc[df_v['Spot frame'] == min_frame].item()
    if first_id in df_e['Link target IDs'].to_numpy().astype(int) and first_id not in df_e['Link target IDs.1'].to_numpy().astype(int):
        return True
    elif first_id not in df_e['Link target IDs'].to_numpy().astype(int) and first_id in df_e['Link target IDs.1'].to_numpy().astype(int):
        return False
        

def df_to_daughter_cells_tracks_source_id(df_e, df_v):
    """Getting daughter cells tracks from the DataFrames"""
    list_mother_cells_e, list_mother_cells_v = df_to_mother_cells_tracks(df_e, df_v)
    list_daughter_cells = []
    for i in range(int(df_v['Spot track ID'].to_numpy().astype(int).max()) + 1):
        if(np.any(list_mother_cells_v[i]['Spot N links'].to_numpy().astype(int)==3)):    #Check if this cell divides in a track ID  
            n_divisions= len(np.where(list_mother_cells_v[i]['Spot N links'].
                                      to_numpy().astype(int)==3)[0])
            n_cells = n_divisions + 1
            list_daughter_cells_temp = [];
            for j in range(n_divisions+1): #Create N divisions + 1 cells and assign them the first time frame
               
                list_daughter_cells_temp.append(pd.DataFrame({'A' : []}))
                list_daughter_cells_temp[-1] = (list_daughter_cells_temp[-1].
                                                append(list_mother_cells_v[i].iloc[0]))

    
            index_id = list_mother_cells_v[i]['ID'][1]
            index_div = 0
            index_cell = 0
            list_mother_cells_e[i].sort_values(by='Link target IDs', 
                                               ignore_index=True); # sort edges by link target IDs.1. When a cell divides, it has two Link target IDs.1
            switches=[]
            for j in range(n_divisions):
                switches.append(0)
                 
            while index_cell < n_cells:
                list_daughter_cells_temp[index_cell] = (list_daughter_cells_temp[index_cell].
                                                        append(list_mother_cells_v[i].
                                                               loc[list_mother_cells_v[i]['ID'] == index_id])) # Append the next frame based on ID. 
              
                if (list_mother_cells_v[i]['Spot N links'].
                    loc[list_mother_cells_v[i]['ID'] == index_id].
                    to_numpy().astype(int) == 3):
                    # division
                    index_id = (list_mother_cells_e[i]['Link target IDs.1'].
                                loc[list_mother_cells_e[i]['Link target IDs'] == index_id]
                                .iloc[switches[index_div]])
                    index_div = index_div + 1
                 
                    
                elif (list_mother_cells_v[i]['Spot N links'].
                      loc[list_mother_cells_v[i]['ID'] == index_id]
                      .to_numpy().astype(int) == 1):
                    # end of cell track

                    if switches[index_div - 1] == 1:
                        switches[index_div - 2] = 1
                        switches[index_div - 1] = 0
                    else:
                        switches[index_div - 1 ] = 1
                        
                        
                    index_div=0
                    index_id =list_mother_cells_v[i]['ID'][1]
                    
                    list_daughter_cells.append(pd.DataFrame({'A' : []}))
                    list_daughter_cells[-1] = list_daughter_cells[-1].append(list_daughter_cells_temp[index_cell])
                    list_daughter_cells[-1] = list_daughter_cells[-1].drop(['A'], axis=1)
                    index_cell = index_cell+1
                    
                else:
                    index_id = list_mother_cells_e[i]['Link target IDs.1'].loc[list_mother_cells_e[i]['Link target IDs'] == index_id].to_numpy().astype(int)[0]

        else:
            list_daughter_cells.append(list_mother_cells_v[i])

    return list_daughter_cells  


def df_to_daughter_cells_tracks_source_id1(df_e, df_v):
    """Getting daughter cells tracks from the DataFrames"""
    list_mother_cells_e, list_mother_cells_v = df_to_mother_cells_tracks(df_e, df_v)
    list_daughter_cells = []

    for i in range(int(df_v['Spot track ID'].to_numpy().astype(int).max()) + 1):
        if(np.any(list_mother_cells_v[i]['Spot N links'].to_numpy().astype(int)==3)):    #Check if this cell divides in a track ID     
            n_divisions= len(np.where(list_mother_cells_v[i]['Spot N links'].
                                      to_numpy().astype(int)==3)[0])
            n_cells = n_divisions+1
            list_daughter_cells_temp = [];
                
            for j in range(n_divisions+1): #Create N divisions + 1 cells and assign them the first time frame
               
                list_daughter_cells_temp.append(pd.DataFrame({'A' : []}))
                list_daughter_cells_temp[-1] = (list_daughter_cells_temp[-1].
                                                append(list_mother_cells_v[i].iloc[0]))
            
            index_id = list_mother_cells_v[i]['ID'][1]
            index_div = 0
            index_cell = 0
            list_mother_cells_e[i].sort_values(by='Link target IDs.1', 
                                               ignore_index=True)
            switches=[]
            for j in range(n_divisions):
                switches.append(0)
                 
            while index_cell < n_cells:
                list_daughter_cells_temp[index_cell] = (list_daughter_cells_temp[index_cell].
                                                        append(list_mother_cells_v[i].
                                                               loc[list_mother_cells_v[i]['ID'] == index_id])) # Append the next frame based on ID. 
                if (list_mother_cells_v[i]['Spot N links'].
                    loc[list_mother_cells_v[i]['ID'] == index_id].
                    to_numpy().astype(int) == 3):
                    # division
                    index_id = (list_mother_cells_e[i]['Link target IDs'].
                                loc[list_mother_cells_e[i]['Link target IDs.1'] == index_id]
                                .iloc[switches[index_div]])
                    index_div = index_div + 1
                 
                    
                elif (list_mother_cells_v[i]['Spot N links'].
                      loc[list_mother_cells_v[i]['ID'] == index_id]
                      .to_numpy().astype(int) == 1):
                    # end of cell track
                    
                    if switches[index_div - 1] == 1:
                        switches[index_div - 2] = 1
                        switches[index_div - 1] = 0
                    else:
                        switches[index_div - 1 ] = 1
                        
                        
                    index_div=0
                    index_id =list_mother_cells_v[i]['ID'][1]
                    
                    list_daughter_cells.append(pd.DataFrame({'A' : []}))
                    list_daughter_cells[-1] = list_daughter_cells[-1].append(list_daughter_cells_temp[index_cell])
                    list_daughter_cells[-1] = list_daughter_cells[-1].drop(['A'], axis=1)
                    index_cell = index_cell+1
                else:
                    index_id = list_mother_cells_e[i]['Link target IDs'].loc[list_mother_cells_e[i]['Link target IDs.1'] == index_id].to_numpy().astype(int)[0]
            

        else:
            list_daughter_cells.append(list_mother_cells_v[i])

    return list_daughter_cells  


def df_to_daughter_cells_tracks(df_e, df_v):
    
    if is_source_id(df_e, df_v):
        try:
            list_daughter_cells = df_to_daughter_cells_tracks_source_id(df_e, df_v)
        except:
            list_daughter_cells = df_to_daughter_cells_tracks_source_id1(df_e, df_v)
        finally: 
            return list_daughter_cells
    else:
        try:
            list_daughter_cells = df_to_daughter_cells_tracks_source_id1(df_e, df_v)
        except:
            list_daughter_cells = df_to_daughter_cells_tracks_source_id(df_e, df_v)
        finally: 
            return list_daughter_cells

        
def select_cells_timepoints(list_cells_as_df, from_timepoint, to_timepoint):
    
    """ Select cells and trim them so that their timepoints are defined"""
    selected_cells = []
    for i, _ in enumerate(list_cells_as_df):
        cell_temp = list_cells_as_df[i].loc[list_cells_as_df[i]['Spot frame'].astype('int') >= from_timepoint]
        selected_cells.append(cell_temp.loc[cell_temp['Spot frame'].astype('int') <= to_timepoint])
        if(selected_cells[-1].shape[0] != (to_timepoint - from_timepoint + 1)):
            del selected_cells[-1]
    return selected_cells


def select_cells_timepoints_undefined(list_cells_as_df, from_timepoint, to_timepoint):
    
    """ Select cells and trim them so that their timepoints are defined"""
    selected_cells = []
    for i, df_cells in enumerate(list_cells_as_df):
        if df_cells['Spot frame'].astype('int').min() <= from_timepoint and df_cells['Spot frame'].astype('int').max() >= to_timepoint:
            selected_cells.append(df_cells)
    return selected_cells


def select_cells_tag(list_cells_as_df, tag):
    selected_cells = []
    for i, _ in enumerate(list_cells_as_df):
        if(list_cells_as_df[i][tag].tail(1).astype('int').item() == 1):
            selected_cells.append(list_cells_as_df[i])
    return selected_cells


def select_cells_positive_for_one_cat_inst(list_cells_as_df, cat_name, list_inst_tag):
    
    """Select cells that are positive for one instance 
    of the category cat_name """
    selected_cells = []
    df = pd.concat(list_cells_as_df)

    for ind, df in enumerate(list_cells_as_df):
        for inst_tag in list_inst_tag:
            if df[inst_tag].tail(1).astype('int').item() == 1:
                selected_cells.append(list_cells_as_df[ind])
    return selected_cells


def trim_df_v(df_v, n_channels):
    """Remove unnecessary columns from dataframes"""
    list_to_drop = ['Detection quality', 'Spot intensity', 'Spot quick mean', 'Spot quick mean.1', 'Spot radius']
    for i in list_to_drop:
        if i in df_v.columns: 
            df_v = df_v.drop([i], axis=1)
    for i in range(1,6*n_channels):
        df_v = df_v.drop(['Spot intensity.'+str(i)], axis=1)
    return df_v


def add_cell_idx_to_df(list_cells_as_df):
    for i, _ in enumerate(list_cells_as_df):
        df_temp = pd.DataFrame({'Cell Idx': np.ones(list_cells_as_df[i].shape[0])*i})
        list_cells_as_df[i].reset_index(inplace=True, drop=True)
        list_cells_as_df[i] = list_cells_as_df[i].join(df_temp)
    return list_cells_as_df


def making_df_v_unit_coherent(df_v):
    
    """Making the variable types coherent, for 2 channels"""
    
    df_v['Spot center intensity'] = df_v['Spot center intensity'].astype(float)
    df_v['Spot center intensity.1'] = df_v['Spot center intensity.1'].astype(float)
    if 'Spot center intensity.2' in df_v.columns:
        df_v['Spot center intensity.2'] = df_v['Spot center intensity.2'].astype(float)
    df_v['Spot N links'] = df_v['Spot N links'].astype(int)
    df_v['Spot frame'] = df_v['Spot frame'].astype(int)
    df_v['Spot position'] = df_v['Spot position'].astype(float)
    df_v['Spot position.1'] = df_v['Spot position.1'].astype(float)
    df_v['Spot position.2'] = df_v['Spot position.2'].astype(float)
    return df_v


def making_list_cells_df_unit_coherent(list_cells_df):
    
    """Making the variable types coherent, for 2 channels"""
    
    for ind, df_v in enumerate(list_cells_df):
        list_cells_df[ind]['Spot center intensity'] = df_v['Spot center intensity'].astype(float)
        list_cells_df[ind]['Spot center intensity.1'] = df_v['Spot center intensity.1'].astype(float)
        if 'Spot center intensity.2' in df_v.columns:
            list_cells_df[ind]['Spot center intensity.2'] = df_v['Spot center intensity.2'].astype(float)
            list_cells_df[ind].loc[list_cells_df[ind]['Spot center intensity.2'] == 0.0]['Spot center intensity.2']  = 1.0
        list_cells_df[ind]['Spot N links'] = df_v['Spot N links'].astype(int)
        list_cells_df[ind]['Spot frame'] = df_v['Spot frame'].astype(int)
        list_cells_df[ind]['Spot position'] = df_v['Spot position'].astype(float)
        list_cells_df[ind]['Spot position.1'] = df_v['Spot position.1'].astype(float)
        list_cells_df[ind]['Spot position.2'] = df_v['Spot position.2'].astype(float)
        list_cells_df[ind]['Cell Idx'] = df_v['Cell Idx'].astype(int)
    return  list_cells_df


def making_df_e_unit_coherent(df_e):
    
    """Making the variable types coherent"""
    df_e['Link target IDs'] = df_e['Link target IDs'].to_numpy().astype(float).astype(int)
    df_e['Link target IDs.1'] = df_e['Link target IDs.1'].to_numpy().astype(float).astype(int)
    df_e['ID'] = df_e['ID'].to_numpy().astype(int)
    return df_e


def interpolate_missing_data(list_cells_as_df):
    for ind, df_cells in enumerate(list_cells_as_df):
        if df_cells['Spot frame'].max() - df_cells['Spot frame'].min() + 1 > len(df_cells):
            list_missing_timepoint = []
            timepoints = np.arange(df_cells['Spot frame'].min(), df_cells['Spot frame'].max()+1, dtype=int)

            for i in timepoints:

                if i not in df_cells['Spot frame'].to_numpy().astype(int):
                    list_missing_timepoint.append(i)
        
             
            for i in list_missing_timepoint:
                df_temp =  list_cells_as_df[ind]
                series_temp = df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i - 1]
                series_temp['Spot center intensity'] = np.mean([df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i - 1]['Spot center intensity'], 
                                                                df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i + 1]['Spot center intensity']])
                series_temp['Spot center intensity.1'] = np.mean([df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i - 1]['Spot center intensity.1'], 
                                                                  df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i + 1]['Spot center intensity.1']])
                if 'Spot center intensity.2' in df_temp.columns:
                    series_temp['Spot center intensity.2'] = np.mean([df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i - 1]['Spot center intensity.2'], 
                                                                      df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i + 1]['Spot center intensity.2']])
                series_temp['Spot position'] = np.mean([df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i - 1]['Spot position'], 
                                                        df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i + 1]['Spot position']])
                series_temp['Spot position.1'] = np.mean([df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i - 1]['Spot position.1'], 
                                                          df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i + 1]['Spot position.1']])
                series_temp['Spot position.2'] = np.mean([df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i - 1]['Spot position.2'], 
                                                          df_temp.loc[df_temp['Spot frame'].to_numpy().astype(int) == i + 1]['Spot position.2']])
                series_temp['Spot frame'] = i
                list_cells_as_df[ind] = df_temp.append(series_temp)
                list_cells_as_df[ind].reset_index(drop = True)
                list_cells_as_df[ind] = list_cells_as_df[ind].sort_values('Spot frame')
    return list_cells_as_df
            
             
def csv_to_clean_df(filename_e, filename_v, from_timepoint, to_timepoint, n_channels, tag = None, defined = True):
    
    """Function to process csv files from Mastodon into usable Dataframes"""
    df_e = csv_to_df(filename_e)
    df_e = making_df_e_unit_coherent(df_e)
    if check_df_e_ok(df_e):
        
        df_v = csv_to_df(filename_v)
        df_v = trim_df_v(df_v, n_channels)
        df_v = making_df_v_unit_coherent(df_v)
        list_cells_as_df = df_to_daughter_cells_tracks(df_e, df_v)
        df_v = pd.concat(list_cells_as_df)
        if tag:
            list_cells_as_df = select_cells_tag(list_cells_as_df, tag)
        list_cells_as_df = interpolate_missing_data(list_cells_as_df)
        if defined:
            list_cells_as_df = select_cells_timepoints(list_cells_as_df, from_timepoint, to_timepoint)
        else:
            list_cells_as_df = select_cells_timepoints_undefined(list_cells_as_df, from_timepoint, to_timepoint)
        list_cells_as_df = add_cell_idx_to_df(list_cells_as_df)
        list_cells_as_df = making_list_cells_df_unit_coherent(list_cells_as_df)
        return list_cells_as_df, pd.concat(list_cells_as_df)
    else:
        print('Error with CSV file: ' +filename_e)
  
        
def csv_to_clean_df_with_tag(filename_e, filename_v, from_timepoint, to_timepoint, n_channels, list_cat = None, defined = True):
    
    """Function to process csv files from Mastodon into usable Dataframes"""
    df_e = csv_to_df(filename_e)
    df_e = making_df_e_unit_coherent(df_e)
    if check_df_e_ok(df_e):
        df_v = pd.read_csv(filename_v)
        df_v = pd.DataFrame(df_v)
        if list_cat:
            list_inst_name = []
            list_inst_tag = []
            for ind, cat_name in enumerate(list_cat):
                name_cat_instances, tag_cat_instances = get_cat_instances_name_and_tag(df_v, cat_name)
                list_inst_name.append(name_cat_instances)
                list_inst_tag.append(tag_cat_instances)

        df_v = csv_to_df(filename_v)
        df_v = trim_df_v(df_v, n_channels)
        df_v = making_df_v_unit_coherent(df_v)
        df_temp_file_id = pd.DataFrame({'source file': [filename_v] * len(df_v)})
        df_v.reset_index(inplace = True, drop = True)
        df_v = df_v.join(df_temp_file_id)
        list_cells_as_df = df_to_daughter_cells_tracks(df_e, df_v)
        df_v = pd.concat(list_cells_as_df)
        if list_cat:
            for ind, cat_name in enumerate(list_cat):
                list_cells_as_df = select_cells_positive_for_one_cat_inst(list_cells_as_df, cat_name, list_inst_tag[ind])
                
        list_cells_as_df = interpolate_missing_data(list_cells_as_df)
        if defined:
            list_cells_as_df = select_cells_timepoints(list_cells_as_df, from_timepoint, to_timepoint)
        else:
            list_cells_as_df = select_cells_timepoints_undefined(list_cells_as_df, from_timepoint, to_timepoint)
        list_cells_as_df = add_cell_idx_to_df(list_cells_as_df)
        list_cells_as_df = making_list_cells_df_unit_coherent(list_cells_as_df)
        
        for ind, cat_name in enumerate(list_cat):
            list_cells_as_df = tag_to_category_name(list_cells_as_df[:], cat_name, list_inst_name[ind])
            
        return list_cells_as_df, pd.concat(list_cells_as_df)
    else:
        print('Error with CSV file: ' +filename_e)
        
def get_number_cat_instances(df, cat_name):
    
    count = 1
    for i in range(1, df.shape[1]):
        if cat_name+'.'+str(i) in df.columns.values:
            count += 1
    return count
 
def get_cat_instances_name_and_tag(df, cat_name):
    
    list_name_cat_instances = [str(df[cat_name].iloc[0])]
    list_tag_cat_instances = [cat_name]
    n_cond = get_number_cat_instances(df, cat_name)

    for i in range(1, n_cond):
        list_tag_cat_instances.append(cat_name +'.'+str(i))
        list_name_cat_instances.append(str(df[cat_name +'.'+str(i)].iloc[0]))

    return list_name_cat_instances, list_tag_cat_instances


def tag_to_category_name(list_df, cat_name,  list_name_cat_instances):
    
    """ Transform the tag system from Mastodon(X.1 == 1, ...)
    to a system with category names (e.g. Somite 18)"""
    
    df = pd.concat(list_df)
    n_cond = len(list_name_cat_instances)
    list_tag_cat_instances = [cat_name]
    for i in range(1, n_cond):
        list_tag_cat_instances.append(cat_name +'.'+str(i))
    
    for ind, df in enumerate(list_df):
        for inst_tag, inst_name in zip(list_tag_cat_instances, list_name_cat_instances):
            if df[inst_tag].astype(int).tail(1).item() ==1 :
                df_temp = pd.DataFrame({cat_name + ' ': [inst_name] * len(df)})
                list_df[ind].reset_index(inplace=True, drop=True)
                list_df[ind] = list_df[ind].join(df_temp)
    return list_df


def combine_categories(list_df, cat_name_1, cat_name_2, cat_name_combined):
    for ind, df in enumerate(list_df):
        df_temp = pd.DataFrame({cat_name_combined : df[cat_name_1].astype(str) + ' ' +  df[cat_name_2].astype(str)})
        list_df[ind].reset_index(drop = True)
        list_df[ind] = list_df[ind].join(df_temp)
    return list_df




def df_to_list_df(df, tag = 'Cell Idx'):
    
    list_df = []
    for tag_index in df[tag].astype(int).drop_duplicates():
        df_temp = df.loc[df[tag] == tag_index]
        df_temp = df_temp.reset_index()
        list_df.append(df_temp)
        
    return list_df

    