#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import copy


######### DATA IMPORT AND GENERAL CLEAN-UP #########

#importing soil horizon data
df_horizon_orig = pd.read_csv('WISE3_HORIZON.csv', sep=';', decimal=',')

#pick relevant columns
df_hz = df_horizon_orig.loc[:,['WISE3_ID','HONU','DESIG','TOPDEP','BOTDEP','ORGC','BULKDENS']]


#drop layers (rows) for which we can't have missing values
df_hz.dropna(subset=['BULKDENS', 'ORGC','TOPDEP','BOTDEP'], inplace=True)
df_hz.reset_index(drop=True, inplace=True)


#keeping layers below ground
df_hz = df_hz.query('BOTDEP > TOPDEP')


#keeping only non-negative ORGC and BULKDENS
df_hz = df_hz.query('ORGC >= 0')
df_hz = df_hz.query('BULKDENS > 0') 
    

#replace missing DESIG with empty string
df_hz['DESIG'].fillna('', inplace = True) 
#apply infer_desig function
df_hz['INFER_DESIG'] = df_hz.apply(lambda x: infer_designation(x['DESIG'], x['BULKDENS']), axis=1)


#keep only unique rows
df_hz.drop_duplicates()


######### CONTINUITY CHECK #########

#generate continuity mask for the whole dataframe
layers_keep = []
unique_profiles = list(dict.fromkeys(df_hz.WISE3_ID.values)) #this method preserves order

for profile in unique_profiles:
    
    layers_keep.extend(continuity_check(df_hz.loc[df_hz['WISE3_ID'] == profile,:]))
        
df_hz_cont = df_hz.loc[layers_keep,:]

del profile
del layers_keep
del unique_profiles


######### CLIPPING PROFILES AT DEPTH 1M #########

unique_profiles_cont = list(dict.fromkeys(df_hz_cont.WISE3_ID.values))
df_hz_clip = pd.DataFrame(columns=list(df_hz_cont.columns))


for profile in unique_profiles_cont:
    
    df_hz_clip = df_hz_clip.append(depth_clip(df_hz_cont.loc[df_hz_cont['WISE3_ID'] == profile,:]), ignore_index=True)
        
    
del profile
del unique_profiles_cont
    
######### COMPUTING ORGC STOCK FOR CLIPPED PROFILES #########

unique_profiles_clip = list(dict.fromkeys(df_hz_clip.WISE3_ID.values))
orgc_stock_list = []


for profile in unique_profiles_clip:
    
    orgc_stock_list.append(orgc_stock(df_hz_clip.loc[df_hz_clip['WISE3_ID'] == profile,:]))
        
    
df_orgc_stock_clip = pd.DataFrame.from_records(orgc_stock_list)
df_orgc_stock_clip.columns = ['profile_id','orgc_stock_mineral','orgc_stock_organic']
    
del profile
del orgc_stock_list 

######### PROFILE DATA IMPORT AND GENERAL CLEAN-UP #########

#importing data from csv
df_site_orig = pd.read_csv("WISE3_SITE.csv", sep=';', decimal=',')

#pick relevant columns
df_site = df_site_orig.loc[:,['WISE3_id','DATEYR','LATDD','LONDD']]

#drop layers (rows) for which we can't have missing values
df_site.dropna(subset=['LATDD', 'LONDD'], inplace=True)
        
df_site.rename(columns={'WISE3_id':'profile_id'}, inplace=True)
        


######### COMPUTING ORG/MIN FRACTIONS #########

fractions_list = []

for profile in unique_profiles_clip:
    
    fractions_list.append(desig_fraction(df_hz_clip.loc[df_hz_clip['WISE3_ID'] == profile,:]))
        
    
df_fractions_clip = pd.DataFrame.from_records(fractions_list)
df_fractions_clip.columns = ['profile_id','top_cm','bot_cm','mineral_fraction','organic_fraction']
    
del profile
del unique_profiles_clip
del fractions_list

######### FINAL MERGE AND CSV EXPORT #########

df_export = df_orgc_stock_clip.merge(df_fractions_clip, on='profile_id').merge(df_site, on='profile_id')

df_export.to_csv('orgc_stock_ISRIC_WISE_clip.csv')


######### FUNCTIONS #########


def desig_fraction(profile):
#function to retrieve min/max depth for profile + mineral & organic fractions 

    prf_id = profile.WISE3_ID.values[1]
    des = profile.INFER_DESIG.values
    top = profile.TOPDEP.values
    bot = profile.BOTDEP.values
                
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

    prf_id = profile.WISE3_ID.values[1]
    top = profile.TOPDEP.values
    bot = profile.BOTDEP.values
    des = profile.INFER_DESIG.values
    
    blk = profile.BULKDENS.values
    orgc = profile.ORGC.values
            
    layer_height = [b-t for b,t in zip(bot,top)]
    num_layers = len(top)
    
    orgc_stock_mineral = 0.0
    orgc_stock_organic = 0.0
    
    for l in range(0,num_layers):
        if des[l] == "organic":
            orgc_stock_organic += blk[l]*orgc[l]*layer_height[l]/100
            
        if des[l] == "mineral":
            orgc_stock_mineral += blk[l]*orgc[l]*layer_height[l]/100   
    
    return [prf_id, orgc_stock_mineral, orgc_stock_organic]
    
    

def depth_clip(profile, depth_min=0, depth_max=100):
    #takes in a profile with N layers and returns a clipped profile, by default from 0 to 100cm
    #the function assumes layers continuity
        
    top = profile.TOPDEP.values
    bot = copy.deepcopy(profile.BOTDEP.values) #create copy of the original profile bot, to be edited later
    
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
            
                            
        profile.BOTDEP = bot
    
    return profile.loc[layers_keep,:]
    


def continuity_check(profile):
#the function assumes layers are ordered (e.g.  0 -> 20 -> 40)
        
    top = profile.TOPDEP.values
    bot = profile.BOTDEP.values
    
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