#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import copy



######### FUNCTIONS #########
def fix_layer_top(layer_top):
    
    if isinstance(layer_top, str):
        if len(layer_top) > 4:
            return int(layer_top[0:2])
        else:
            return int(layer_top)
    else:
        return layer_top


def desig_fraction(profile):
    #function to retrieve min/max depth for profile + mineral & organic fractions 

    prf_id = profile.sp_name.values[0]
    des = profile.INFER_DESIG.values
    top = profile.layer_top.values
    bot = profile.layer_bot.values
                
    num_layers = len(top)
        
    mineral_index = [i for i in range(0,num_layers) if des[i] == "mineral"]
    organic_index = [i for i in range(0,num_layers) if des[i] == "organic"]
            
    layers_height = [b-t for b,t in zip(bot,top)]

    #values to return

    profile_top = min(top)
    profile_bot = max(bot)
        
    profile_height = profile_bot - profile_top
    
    mineral_fraction = 0.0 if len(mineral_index)==0 else sum([layers_height[i] for i in mineral_index])/profile_height
    
    organic_fraction = 0.0 if len(organic_index)==0 else sum([layers_height[i] for i in organic_index])/profile_height
        
    return [prf_id, profile_top, profile_bot, mineral_fraction, organic_fraction]
            


def orgc_stock(profile):

    prf_id = profile.sp_name.values[0]
    top = profile.layer_top.values
    bot = profile.layer_bot.values
    des = profile.INFER_DESIG.values
    
    blk_org = profile.bd_organic.values  #[g/cm³]
    blk_min = profile.bd_mineral.values

    orgc = profile.oc_percent.values #[%]
            
    layer_height = [b-t for b,t in zip(bot,top)]
    num_layers = len(top)
    
    orgc_stock_mineral = 0.0
    orgc_stock_organic = 0.0
    
    for l in range(0,num_layers):
        if des[l] == "organic":
            orgc_stock_organic += (100*100*layer_height[l])*(blk_org[l]*orgc[l]/100)/1000
            
        if des[l] == "mineral":
            orgc_stock_mineral += (100*100*layer_height[l])*(blk_min[l]*orgc[l]/100)/1000
    
    return [prf_id, orgc_stock_mineral, orgc_stock_organic]
    
    

def depth_clip(profile, depth_min=0, depth_max=30):
    #takes in a profile with N layers and returns a clipped profile
    #the function assumes layers continuity
        
    top = profile.layer_top.values
    bot = copy.deepcopy(profile.layer_bot.values) #create copy of the original profile bot, to be edited later
    
    num_layers = len(top)
    
    #init array of 'true' -> keep all layers unless proven guilty
    layers_keep = [True for i in range(0,num_layers)]
    
    #if profile doesn't reach threshold, or starts after min, ditch the whole profile
    if (max(bot) < depth_max) or (min(top) > depth_min):    
                
        layers_keep = [False for i in range(0,num_layers)]
        
    #othewise mark layers to drop above threshold, and clip the overlapping one
    else:
        for y in range(0, num_layers):
            #simple case, both top and bot below threshold
            if top[y] >= depth_max:
                layers_keep[y] = False
    
                #layers overlapping threshold
            if (top[y] < depth_max) and (bot[y] > depth_max):
                bot[y] = depth_max
            
                            
        profile.layer_bot = bot
    
    return profile.loc[layers_keep,:]
    


def continuity_check(profile):
    #the function assumes layers are ordered (e.g.  0 -> 20 -> 40)
        
    top = profile.layer_top.values
    bot = profile.layer_bot.values
    
    num_layers = len(top)
    
    #1st layer is always good
    mask = [True]
    
    #from 2nd layer downward, we check for matching depths at bot/top
    mask.extend([top[y] == bot[y-1] for y in range(1,num_layers)])
        
    #if we run into a non-matching layer, we mark as 'false' all following ones
    first_noMatch_id = next((i for (i,x) in enumerate(mask) if x == False), False)
    
    if first_noMatch_id != False: 
        for i in range(first_noMatch_id, num_layers):
            mask[i] = False
    
    return mask    



def infer_designation(old_desig, blk):
    #takes layer and returns its inferred soil designation, either 'organic' or 'mineral'
        
    new_desig = "no_desig" #value to be eventually returned
    threshold = 0.5
    
    control_desig = 'OoHh' #possible chars in existing desig that indicate presence of organic material in the layer
    
    new_desig_organic = "organic"
    new_desig_mineral = "mineral"
        
    #if designation data is missing (set to '' for NA, N/A, and other variants depending on the dataset)
    if old_desig == '':
        #check on bulk density value to decide which designation applies 
        new_desig = new_desig_organic if (blk < threshold) else new_desig_mineral
        
    #designation not missing            
    else    :
        check_desig = any(c in old_desig for c in control_desig) #check if any char in the designation matches control
        check_desig_BD = check_desig and (blk < threshold) #true if any char matches and BD < treshold
        new_desig = new_desig_organic if check_desig_BD else new_desig_mineral
    
        
    return new_desig


######### DATA IMPORT AND GENERAL CLEAN-UP #########

#importing soil horizon data (all layers)
df_layers_orig_c1 = pd.read_csv('ISCN_ALL_DATA_LAYER_C1_1-1.csv', sep=';', decimal=',')
df_layers_orig_c2 = pd.read_csv('ISCN_ALL_DATA_LAYER_C2_1-1.csv', sep=';', decimal=',')
df_layers_orig_c3 = pd.read_csv('ISCN_ALL_DATA_LAYER_C3_1-1.csv', sep=';', decimal=',')
df_layers_orig_c4 = pd.read_csv('ISCN_ALL_DATA_LAYER_C4_1-1.csv', sep=';', decimal=',')

#generate unique df with all imported layers
df_all_orig = df_layers_orig_c1.append(df_layers_orig_c2, ignore_index=True)
df_all_orig = df_all_orig.append(df_layers_orig_c3, ignore_index=True)
df_all_orig = df_all_orig.append(df_layers_orig_c4, ignore_index=True)


#define an unique identifier for each profile
df_all = copy.deepcopy(df_all_orig)
df_all['sp_name'] =  [row['site_name'] + '_' + row['profile_name'] for index, row in df_all.iterrows()]


#pick relevant columns
df_all = df_all.loc[:,['sp_name','layer_name','lat (dec. deg)','long (dec. deg)','observation_date (YYYY-MM-DD)','hzn','hzn_desgn','layer_top (cm)','layer_bot (cm)','oc (percent)']]

#renaming
df_all.rename(columns={'layer_bot (cm)':'layer_bot'}, inplace=True)
df_all.rename(columns={'layer_top (cm)':'layer_top'}, inplace=True)
df_all.rename(columns={'lat (dec. deg)':'lat_dd'}, inplace=True)
df_all.rename(columns={'long (dec. deg)':'long_dd'}, inplace=True)
df_all.rename(columns={'oc (percent)':'oc_percent'}, inplace=True)
df_all.rename(columns={'observation_date (YYYY-MM-DD)':'obs_date_yyyy'}, inplace=True)


#drop layers (rows) for which we can't have missing values
df_all.dropna(subset=['layer_top','layer_bot','oc_percent','lat_dd', 'long_dd'], inplace=True)
df_all.reset_index(drop=True, inplace=True)


#some layer_top values are wrongly represented as "30/01/1900 0:00" instead of just 30
#these values (strings) are longer than 4 chars, keep only first two
df_all['layer_top'] = df_all.apply(lambda x: fix_layer_top(x['layer_top']), axis=1)

#keeping only layers below ground
df_all = df_all.query('layer_bot > layer_top')

#keeping only non-negative ORGC
df_all = df_all.query('oc_percent >= 0')


#dropping layers without date information
df_all['obs_date_yyyy'].replace('', np.nan, inplace=True)
df_all.dropna(subset=['obs_date_yyyy'], inplace=True)

#uniforming date information to year-only
df_all['obs_date_yyyy'] = df_all.apply(lambda x: x['obs_date_yyyy'][-4:], axis=1)

#keep only unique rows
df_all.drop_duplicates(inplace=True)


######### CONTINUITY CHECK #########

#generate continuity mask for the whole dataframe
layers_keep = []
unique_profiles = list(dict.fromkeys(df_all.sp_name.values)) #this method preserves order

for profile in unique_profiles:
    
    layers_keep.extend(continuity_check(df_all.loc[df_all['sp_name'] == profile,:]))
        
df_all_cont = df_all.loc[layers_keep,:]

del profile
del layers_keep
del unique_profiles


######### CLIPPING PROFILES AT DEPTH 30 CM #########

unique_profiles_cont = list(dict.fromkeys(df_all_cont.sp_name.values))
df_all_clip = pd.DataFrame(columns=list(df_all_cont.columns))


for profile in unique_profiles_cont:
    
    df_all_clip = df_all_clip.append(depth_clip(df_all_cont.loc[df_all_cont['sp_name'] == profile,:]), ignore_index=True)
        
    
del profile
del unique_profiles_cont


######### IMPORTING AND CLEANING BULK_DENSITY ESTIMATIONS #########

df_bd_pred_orig = pd.read_csv('ISCNTemplate_NRCS_BD_predictions.csv', sep=';', decimal=',')

#define an unique identifier for each profile
df_bd = copy.deepcopy(df_bd_pred_orig)
df_bd['sp_name'] =  [row['site_name'] + '_' + row['profile_name'] for index, row in df_bd.iterrows()]

#pick relevant columns
df_bd = df_bd.loc[:,['sp_name','layer_name','bd_pred_1','bd_pred_2']]

#drop layers (rows) for which we can't have missing values
df_bd.dropna(subset=['bd_pred_1','bd_pred_2'], inplace=True)
df_bd.reset_index(drop=True, inplace=True)

#renaming
df_bd.rename(columns={'bd_pred_1':'bd_organic'}, inplace=True)
df_bd.rename(columns={'bd_pred_2':'bd_mineral'}, inplace=True)

#keep only unique rows
df_bd.drop_duplicates(inplace=True)



######### MERGING LAYERS' DATA WITH BULK_DENSITY ESTIMATIONS #########

df_merged = df_all_clip.merge(df_bd, on=['sp_name','layer_name'], how='inner')

#fixing null values in layer designation columns
df_merged.loc[df_merged.hzn=='?'] = df_merged.replace(to_replace={ 'hzn' : { '?' : '' }})
df_merged['hzn'].fillna('',inplace = True)

df_merged.loc[df_merged.hzn_desgn=='unknown'] = df_merged.replace(to_replace={ 'hzn_desgn' : { 'unknown' : '' }})
df_merged['hzn_desgn'].fillna('',inplace = True)

#create new column with concatenated designations, but leave row empty if both are empty
df_merged['hzn_merged'] =  [row['hzn'] + '_' + row['hzn_desgn'] for index, row in df_merged.iterrows()]

#using bd_organic as it's more conservative
df_merged['INFER_DESIG'] = df_merged.apply(lambda x: infer_designation(x['hzn_merged'], x['bd_organic']), axis=1)



    
######### COMPUTING ORGC STOCK FOR CLIPPED PROFILES #########

unique_profiles_clip = list(dict.fromkeys(df_merged.sp_name.values))
orgc_stock_list = []


for profile in unique_profiles_clip:
    
    orgc_stock_list.append(orgc_stock(df_merged.loc[df_merged['sp_name'] == profile,:]))
        
    
df_orgc_stock_clip = pd.DataFrame.from_records(orgc_stock_list)
df_orgc_stock_clip.columns = ['profile_id','orgc_stock_mineral','orgc_stock_organic']
    
del profile
del orgc_stock_list 


######### COMPUTING ORG/MIN FRACTIONS #########

fractions_list = []

for profile in unique_profiles_clip:
    
    fractions_list.append(desig_fraction(df_merged.loc[df_merged['sp_name'] == profile,:]))
        
    
df_fractions_clip = pd.DataFrame.from_records(fractions_list)
df_fractions_clip.columns = ['profile_id','top_cm','bot_cm','mineral_fraction','organic_fraction']
    
del profile
del unique_profiles_clip
del fractions_list

######### PROFILE DATA EXTRACTION #########

#pick relevant columns
df_profile_info = df_merged.loc[:,['sp_name','obs_date_yyyy','lat_dd','long_dd']]
        
df_profile_info.rename(columns={'sp_name':'profile_id'}, inplace=True)
        
#keep only unique rows
df_profile_info.drop_duplicates(subset='profile_id', inplace=True)

######### FINAL MERGE AND CSV EXPORT #########

df_export = df_orgc_stock_clip.merge(df_fractions_clip, on='profile_id', how='inner').merge(df_profile_info, on='profile_id', how='left')

df_export.to_csv('orgc_stock_ISCN_clip.csv')