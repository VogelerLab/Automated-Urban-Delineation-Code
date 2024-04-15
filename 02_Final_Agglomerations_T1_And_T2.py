# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 12:06:48 2022

@author: orioncr

Description: The purpose of this script is to produce the final agglomerations
for the country of interest for both t1 and t2. The data produced in this script
will be used for calcluations of LCR, PGR, LCRPGR, and other development metrics.
"""
#%% Import Libraries and Initialize QGIS
import os
#import math
#import plotly
import rasterio
#from pathlib import Path
from qgis.core import (
    QgsVectorLayer,
    QgsRasterLayer,
    QgsApplication,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransformContext,
    QgsLayerMetadata,
    QgsMapLayer,
    QgsExpression,
    QgsExpressionContext,
    QgsExpressionContextUtils,
    edit,
    QgsVectorDataProvider, 
    QgsDataProvider,
    QgsField,
    QgsFeatureRequest,
    QgsVectorFileWriter)
#from qgis.PyQt.QtCore import QVariant
#from qgis.utils import iface

#Initialize QGIS
#Supply path to qgis install location
#You can get this using QgsApplication.prefixPath() in QGIS Python console
QgsApplication.setPrefixPath('C:/PROGRA~1/QGIS32~1.4/apps/qgis-ltr', True)
# Create a reference to the QgsApplication.  Setting the
# second argument to False disables the GUI.
qgs = QgsApplication([], False)
# Load providers
qgs.initQgis()

#Add provider to import native algorithms and initialize processing
from qgis.analysis import QgsNativeAlgorithms
import processing
from processing.core.Processing import Processing
Processing.initialize()
#QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())

#%% Set General Variables
#Enter the abbreviation of your country or region of interest here (as a string)
country_abbrv = "XX"

#Enter code for project as a string
#Example projection: Africa Albers Equal Area
projection = 'ESRI:102022'

# %% Directories
#Directory for testing outputs if desired
dir_test = r"" + country_abbrv +"_temp"

dir_00_outputs_t1 = r"Path to 00_outputs_t1 folder"

dir_00_outputs_t2 =  r"Path to 00_outputs_t2 folder"

dir_01_outputs_t1 = r"Path to 01_outputs_t1 folder/"

dir_01_outputs_t2 = r"Path to 01_outputs_t2 folder"

dir_02_outputs = r"Path to new 02_outputs folder" + country_abbrv

if not os.path.exists(dir_02_outputs):
    os.makedirs(dir_02_outputs)
    
# Use this line of code to delete everything in a specific directory.
# You will likely only need to do this for the outputs of the specific script 
# you are working in. Removing files is the best way to 'overwrite'. 
# https://www.techiedelight.com/delete-all-files-directory-python/

#for f in os.listdir(dir_02_outputs):
#    os.remove(os.path.join(dir_02_outputs, f))
#%% Set Characteristic Variables
developed = os.path.join(dir_00_outputs_t1, "Developed_" + country_abbrv + "_t1.tif")

raster = rasterio.open(developed)

cell_size = raster.res[0]

cell_area = cell_size**2    
# %% Assign Filepaths To Variables
unconfigured_aggloms_t1 = os.path.join(dir_01_outputs_t1, "unconfigured_aggloms.shp")

unconfigured_aggloms_t2 = os.path.join(dir_01_outputs_t2, "unconfigured_aggloms.shp")

pg_t1 = os.path.join(dir_01_outputs_t1, "pop_grid_" + country_abbrv + "_t1.tif")

pg_t2 = os.path.join(dir_01_outputs_t2, "pop_grid_" + country_abbrv + "_t2.tif")

developed_pixels_t1 = os.path.join(dir_00_outputs_t1, "Developed_" + country_abbrv + "_t1.tif")

developed_pixels_t2 = os.path.join(dir_00_outputs_t2, "Developed_" + country_abbrv + "_t2.tif") 

clusters_less_threshold_t1 = os.path.join(dir_01_outputs_t1, 'clusters_less_threshold.shp') 

clusters_less_threshold_t2 = os.path.join(dir_01_outputs_t2, 'clusters_less_threshold.shp') 
#%% Processing
# The agglomerations from t1 may not match up with the agglomerations for t2.
# This may be due to cores growing and connecting over time or from clusters that
# grow larger and are closer in distance to other clusters. Since our method uses
# a distance analysis to associate clusters, the associated clusters may change 
# from t1 to t2 thereby changing the entire configuration of the agglomerations between periods
# To address this, we will analyze where overlap occurs between the t1 and t2 agglomerations.
# Single agglomerations from time 2 that overlap with multiple agglomerations from t1 will be
# dissolved. We will then use a spatial join to join the newly dissolved layer for t1 to the
# t2 layer and dissolve.

# Find matching agglomerations for time 1
sj_t1_to_t2_fp = os.path.join(dir_02_outputs, "sj_t1_to_t2.shp")

sj_t1_to_t2 = processing.run("native:joinattributesbylocation", 
               {'INPUT':unconfigured_aggloms_t1,
                'JOIN':unconfigured_aggloms_t2,
                'PREDICATE':[0,1,2,3,4,5,6],
                'JOIN_FIELDS':[],
                'METHOD':2,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'t2_',
                'OUTPUT':sj_t1_to_t2_fp})

dissolve_sj_t1_to_t2_fp = os.path.join(dir_02_outputs, "dsslv_sj_t1_to_t2.shp")

processing.run("native:dissolve",
               {'INPUT':sj_t1_to_t2_fp,
                'FIELD':['t2_Id'],
                'OUTPUT':dissolve_sj_t1_to_t2_fp})

# Find matching agglomerations for time 2
sj_overlap_to_t2_fp = os.path.join(dir_02_outputs, "sj_overlap_to_t2.shp")

sj_t1_to_t2 = processing.run("native:joinattributesbylocation", 
               {'INPUT':unconfigured_aggloms_t2,
                'JOIN':dissolve_sj_t1_to_t2_fp,
                'PREDICATE':[0,1,2,3,4,5,6],
                'JOIN_FIELDS':[],
                'METHOD':2,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'OL_',
                'OUTPUT':sj_overlap_to_t2_fp})

dissolve_sj_overlap_to_t2_fp = os.path.join(dir_02_outputs, "dsslv_sj_overlap_to_t2.shp")

processing.run("native:dissolve",
               {'INPUT':sj_overlap_to_t2_fp,
                'FIELD':['OL_Id'],
                'OUTPUT':dissolve_sj_overlap_to_t2_fp})
#%%Create final agglomerations with full data

# Remove nulls for OL_ID or t2_ID. If these are null, that means that there is 
# not an agglomeration present in one of the years. This is likely due to lack of
# built up data or misclassification of built up data across years.

# time 1 removal of null values
sj_t1_to_t2_vctr = QgsVectorLayer(dissolve_sj_t1_to_t2_fp)

rm_null_t1_aggloms = os.path.join(dir_02_outputs,"t1_rmvd_null_aggloms.shp")

rm_null_t1 = processing.run("qgis:selectbyexpression", 
               {'INPUT':sj_t1_to_t2_vctr,
                'EXPRESSION':' "t2_Id" IS NOT NULL',
                'METHOD':0})

rm_null_t1_output = rm_null_t1['OUTPUT']

processing.run("native:saveselectedfeatures", 
               {'INPUT':rm_null_t1_output,
                'OUTPUT':rm_null_t1_aggloms})

#time 2 removal of null values
sj_ol_to_t2_vctr = QgsVectorLayer(dissolve_sj_overlap_to_t2_fp)

rm_null_t2_aggloms = os.path.join(dir_02_outputs,
                                         "t2_rmvd_null_aggloms.shp")

rm_null_t2 = processing.run("qgis:selectbyexpression", 
               {'INPUT':sj_ol_to_t2_vctr,
                'EXPRESSION':' "OL_Id" IS NOT NULL',
                'METHOD':0})

rm_null_t2_output = rm_null_t2['OUTPUT']

processing.run("native:saveselectedfeatures", 
               {'INPUT':rm_null_t2_output,
                'OUTPUT':rm_null_t2_aggloms})

#%%
# Find non-core clusters from time 1 that overlap with the Final Agglomerations of time 2
# This will ensure that we are not missing clusters that exist but were removed earlier in the analysis
# due to not meeting a threshold
sj_t2_to_t1_ncc_fp = os.path.join(dir_02_outputs, "sj_t2_to_t1_ncc.shp")

sj_t2_to_t1_ncc = processing.run("native:joinattributesbylocation", 
               {'INPUT':clusters_less_threshold_t1,
                'JOIN':rm_null_t2_aggloms,
                'PREDICATE':[0,1,2,3,4,5,6],
                'JOIN_FIELDS':['Id','OL_Id'],
                'METHOD':0,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'NCC_',
                'OUTPUT':sj_t2_to_t1_ncc_fp})

#sj_t2_to_t1_ncc =sj_t2_to_t1_ncc['OUTPUT']

# Extract features that are not NULL
rm_null_fp = os.path.join(dir_02_outputs, "rm_null.shp")

rm_null = processing.run("native:extractbyexpression", 
               {'INPUT':sj_t2_to_t1_ncc_fp,
                'EXPRESSION':'NCC_OL_Id IS NOT NULL',
                'OUTPUT':rm_null_fp})

#rm_null = rm_null['OUTPUT']

# Dissolve by ID field that matches agglomeration and non-core clusters
dissolve_sj_t2_to_t1_ncc = processing.run("native:dissolve",
               {'INPUT':rm_null_fp,
                'FIELD':['NCC_OL_Id'],
                'OUTPUT':'memory:'})

dissolve_sj_t2_to_t1_ncc = dissolve_sj_t2_to_t1_ncc['OUTPUT']

# Merge layers
merge_layers_fp = os.path.join(dir_02_outputs, "merge_layers_t1.shp")

merge_layers = processing.run("native:mergevectorlayers", 
               {'LAYERS':[rm_null_t1_aggloms, dissolve_sj_t2_to_t1_ncc],
                'CRS':, # Fill value
                'OUTPUT':merge_layers_fp})

#merge_layers = merge_layers['OUTPUT']

# Add field to be dissolved
add_field_fp = os.path.join(dir_02_outputs, "add_Field_t1.shp")

add_field = processing.run("native:fieldcalculator", 
               {'INPUT':merge_layers_fp,
                'FIELD_NAME':'New_ID',
                'FIELD_TYPE':1,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'CASE WHEN "t2_Id" IS NULL\r\nTHEN "NCC_Id"\r\nELSE "t2_Id"\r\nEND',
                'OUTPUT':add_field_fp})

add_field = add_field['OUTPUT']

# Dissolve new field
dissolve_field_fp = os.path.join(dir_02_outputs, "dissolve_field_t1.shp")

dissolve_field = processing.run("native:dissolve",
               {'INPUT':add_field_fp,
                'FIELD':['New_ID'],
                'OUTPUT':dissolve_field_fp})

#dissolve_field = dissolve_field['OUTPUT']
# Remove duplicate rows

t1_units_without_calcs = os.path.join(dir_02_outputs, "t1_units_wo_calcs.shp")

rm_rows = rm_null = processing.run("native:extractbyexpression", 
               {'INPUT':dissolve_field_fp,
                'EXPRESSION':' "layer"  != \'output\'',
                'OUTPUT':t1_units_without_calcs})

# os.remove(t1_units_without_calcs)

#%%
# Now we need to do the same for the final time with the time 1 units we just created 
sj_t2_fp = os.path.join(dir_02_outputs, "sj_t2.shp")

sj_t2 = processing.run("native:joinattributesbylocation", 
               {'INPUT':clusters_less_threshold_t2,
                'JOIN':t1_units_without_calcs,
                'PREDICATE':[0,1,2,3,4,5,6],
                'JOIN_FIELDS':['Id','t2_Id'],
                'METHOD':0,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'OUTPUT':sj_t2_fp})

#sj_t2 =sj_t2['OUTPUT']

# Dissolve by ID field that matches agglomeration and non-core clusters
dissolve_sj_t2_fp = os.path.join(dir_02_outputs, "dissolve_sj_t2.shp")

dissolve_sj_t2 = processing.run("native:dissolve",
               {'INPUT':sj_t2_fp,
                'FIELD':['t2_Id'],
                'OUTPUT':dissolve_sj_t2_fp})

#dissolve_sj_t2 = dissolve_sj_t2['OUTPUT']

# Extract features that are not NULL
rm_null_fp = os.path.join(dir_02_outputs, "rm_null_t2.shp")

rm_null = processing.run("native:extractbyexpression", 
               {'INPUT':dissolve_sj_t2_fp,
                'EXPRESSION':'t2_Id IS NOT NULL',
                'OUTPUT':rm_null_fp})

#rm_null = rm_null['OUTPUT']

# Merge layers
merge_layers_fp = os.path.join(dir_02_outputs, "merge_layers_t2.shp") 

merge_layers = processing.run("native:mergevectorlayers", 
               {'LAYERS':[rm_null_t2_aggloms, rm_null_fp],
                'CRS':, # Fill value
                'OUTPUT':merge_layers_fp})

#merge_layers = merge_layers['OUTPUT']

# Add field to be dissolved
add_field = processing.run("native:fieldcalculator", 
               {'INPUT':merge_layers_fp,
                'FIELD_NAME':'New_ID',
                'FIELD_TYPE':1,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'CASE WHEN "OL_t2_Id" IS NULL\r\nTHEN "t2_Id"\r\nELSE "OL_t2_Id"\r\nEND',
                'OUTPUT':'memory:'})

add_field = add_field['OUTPUT']

# Dissolve new field

dissove_field_fp = os.path.join(dir_02_outputs, "dissolve_field_t2.shp")

processing.run("native:dissolve",
               {'INPUT':add_field,
                'FIELD':['New_ID'],
                'OUTPUT':dissove_field_fp})

# Remove duplicate rows if necessary

t2_units_without_calcs = os.path.join(dir_02_outputs, "t2_units_wo_calcs.shp")

rm_rows = processing.run("native:extractbyexpression", 
               {'INPUT':dissove_field_fp,
                'EXPRESSION':' "layer"  != \'output\'',
                'OUTPUT':t2_units_without_calcs})
#%%
# We now have the final agglomerations to which we can calculate population
# and built-up statistics for
# =============================================================================
# 
# uncleaned_final_agglomerations_t1 = os.path.join(dir_02_outputs, country_abbrv + 
#                                              "_t1_Uncleaned_Agglomerations.shp")
# 
# #Use this code if need to remove the file to rerun
# if os.path.isfile(uncleaned_final_agglomerations_t1):
#     os.remove(uncleaned_final_agglomerations_t1)
#    
# =============================================================================

# time 1
# Recalculate population zonal statistics for each unique agglomeration
pop_zonal_stats_t1 = processing.run("native:zonalstatisticsfb", 
               {'INPUT':t1_units_without_calcs,
                'INPUT_RASTER':pg_t1,
                'RASTER_BAND':1,
                'COLUMN_PREFIX':'t1_pop_',
                'STATISTICS':[0,1,2],
                'OUTPUT':'memory:'})

pzs_output_t1 = pop_zonal_stats_t1['OUTPUT']

# Calculate built up area for each unique agglomeration

#Recalculate built-up zonal statistics for each unique agglomeration
bu_zonal_stats_t1 = processing.run("native:zonalstatisticsfb", 
               {'INPUT':pzs_output_t1,
                'INPUT_RASTER':developed_pixels_t1,
                'RASTER_BAND':1,
                'COLUMN_PREFIX':'t1_bu_',
                'STATISTICS':[0,1,2],
                'OUTPUT':'memory:'})

buzs_output_t1 = bu_zonal_stats_t1['OUTPUT']

# Calculate built up area
bu_area_t1 = processing.run("native:fieldcalculator",
               {'INPUT':buzs_output_t1,
                'FIELD_NAME':'t1_bua(m²)',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':' "t1_bu_sum" * {}'.format(cell_area) ,
                'OUTPUT':'memory:'})

bua_output_t1 = bu_area_t1['OUTPUT']

# Zonal stats for t2 in t1 boundaries
bu_zonal_stats_t2_in_t1 = processing.run("native:zonalstatisticsfb", 
               {'INPUT':bua_output_t1,
                'INPUT_RASTER':developed_pixels_t2,
                'RASTER_BAND':1,
                'COLUMN_PREFIX':'t2_in_t1_',
                'STATISTICS':[0,1,2],
                'OUTPUT':'memory:'})

buzs_output_t1 = bu_zonal_stats_t1['OUTPUT']

# Calculate the built up area of t2 in t1 boundaries
bu_area_of_t2_in_t1 = processing.run("native:fieldcalculator",
               {'INPUT':bua_output_t1,
                'FIELD_NAME':'t2_bua_in_t1',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':' "t2_in_t1_sum" * {}'.format(cell_area) ,
                'OUTPUT':'memory:'})

bu_area_of_t2_in_t1_output = bu_area_of_t2_in_t1['OUTPUT']

#Calculate area field for each unique agglomeration

area_final_agglomeration_t1 = processing.run("native:fieldcalculator",
               {'INPUT':bu_area_of_t2_in_t1_output,
                'FIELD_NAME':'area(m)',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'$area',
                'OUTPUT':'memory:'})

afa_output_t1 = area_final_agglomeration_t1['OUTPUT']

ucfa_t1 = processing.run("native:fieldcalculator", 
               {'INPUT':afa_output_t1,
                'FIELD_NAME':'unit_name',
                'FIELD_TYPE':2,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'CASE\r\nWHEN "name_en" IS NULL THEN "name"\r\nELSE "name_en"\r\nEND\r\n\r\n',
                'OUTPUT':'memory:'})

ucfa_t1_output = ucfa_t1['OUTPUT']

#Check field names
field_names_t1 = ucfa_t1_output.fields().names()
print(field_names_t1)

uncleaned_final_agglomerations_t1 = os.path.join(dir_02_outputs,
                                             country_abbrv + "_t1_Uncleaned_Agglomerations.shp")

                                           
processing.run("native:retainfields",
               {'INPUT':ucfa_t1_output,
                'FIELDS':['Id',
                          'New_ID',
                          'full_id',
                          'osm_id',
                          'place',
                          'alt_name',
                          'name_en',
                          'name',
                          't1_pop_count',
                          't1_pop_sum',
                          't1_pop_mean',
                          't1_bu_count',
                          't1_bu_sum',
                          't1_bu_mean',
                          't1_bua(m²)',
                          'area(m)',
                          'unit_name',
                          't2_bua_in_t1'],
                'OUTPUT':uncleaned_final_agglomerations_t1})

clean_final_agglomerations_t1 = os.path.join(dir_02_outputs, country_abbrv + '_t1_Agglomerations.shp')

processing.run("native:refactorfields", 
               {'INPUT':uncleaned_final_agglomerations_t1,
                'FIELDS_MAPPING':[{'expression': '"New_ID"','length': 10,'name': 'ID','precision': 0,'type': 4},
                                  {'expression': '"t1_pop_sum"','length': 23,'name': 't1_Pop','precision': 15,'type': 4},
                                  {'expression': '"t1_bu_sum"','length': 23,'name': 't1_BU_Sum','precision': 15,'type': 4},
                                  {'expression': '"t1_bua(m²)"','length': 23,'name': 't1_BUA(m²)','precision': 15,'type': 4},
                                  {'expression': '"area(m)"','length': 23,'name': 'Area(m²)','precision': 15,'type': 4}],
                'OUTPUT':clean_final_agglomerations_t1})

# Might need this expression for keeping t2_bua_in_t1 but we don't need it for now
# {'expression': '"t2_bua_in_t1"','length': 23,'name': 't2_BUA_t1','precision': 15,'type': 4}

#%%
# We now have the final agglomerations to which we can calculate population
# and built-up statistics for

# time 2
# Recalculate population zonal statistics for each unique agglomeration
pop_zonal_stats_t2 = processing.run("native:zonalstatisticsfb", 
               {'INPUT':t2_units_without_calcs,
                'INPUT_RASTER':pg_t2,
                'RASTER_BAND':1,
                'COLUMN_PREFIX':'t2_pop_',
                'STATISTICS':[0,1,2],
                'OUTPUT':'memory:'})

pzs_output_t2 = pop_zonal_stats_t2['OUTPUT']

# Calculate built up area for each unique agglomeration

#Recalculate built-up zonal statistics for each unique agglomeration
bu_zonal_stats_t2 = processing.run("native:zonalstatisticsfb", 
               {'INPUT':pzs_output_t2,
                'INPUT_RASTER':developed_pixels_t2,
                'RASTER_BAND':1,
                'COLUMN_PREFIX':'t2_bu_',
                'STATISTICS':[0,1,2],
                'OUTPUT':'memory:'})

buzs_output_t2 = bu_zonal_stats_t2['OUTPUT']

# Calculate built up area
bu_area_t2 = processing.run("native:fieldcalculator",
               {'INPUT':buzs_output_t2,
                'FIELD_NAME':'t2_bua(m²)',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':' "t2_bu_sum" * {}'.format(cell_area) ,
                'OUTPUT':'memory:'})

bua_output_t2 = bu_area_t2['OUTPUT']
#Calculate area field for each unique agglomeration
   
area_final_agglomeration_t2 = processing.run("native:fieldcalculator",
               {'INPUT':bua_output_t2,
                'FIELD_NAME':'area(m)',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'$area',
                'OUTPUT':'memory:'})

afa_output_t2 = area_final_agglomeration_t2['OUTPUT']

ucfa_t2 = processing.run("native:fieldcalculator", 	
               {'INPUT':afa_output_t2,	
                'FIELD_NAME':'unit_name',	
                'FIELD_TYPE':2,	
                'FIELD_LENGTH':0,	
                'FIELD_PRECISION':0,	
                'FORMULA':'CASE\r\nWHEN "name_en" IS NULL THEN "name"\r\nELSE "name_en"\r\nEND\r\n\r\n',	
                'OUTPUT':'memory:'})

ucfa_t2_output = ucfa_t2['OUTPUT']

#Check field names	
field_names_t2 = ucfa_t2_output.fields().names()	
print(field_names_t2)	


uncleaned_final_agglomerations_t2 = os.path.join(dir_02_outputs, country_abbrv + 
                                             "_t2_Uncleaned_Agglomerations.shp")

processing.run("native:retainfields",	
               {'INPUT':ucfa_t2_output,	
                'FIELDS':['Id',	
                          'New_ID',	
                          'full_id',	
                          'osm_id',	
                          'place',	
                          'alt_name',	
                          'name_en',	
                          'name',	
                          't2_pop_count',	
                          't2_pop_sum',	
                          't2_pop_mean',	
                          't2_bu_count',	
                          't2_bu_sum',	
                          't2_bu_mean',	
                          't2_bua(m²)',	
                          'area(m)',	
                          'unit_name'],	
                'OUTPUT':uncleaned_final_agglomerations_t2})	

clean_final_agglomerations_t2 = os.path.join(dir_02_outputs, country_abbrv + "_t2_Agglomerations.shp")

processing.run("native:refactorfields", 		
               {'INPUT':uncleaned_final_agglomerations_t2,		
                'FIELDS_MAPPING':[{'expression': '"New_ID"','length': 10,'name': 'ID','precision': 0,'type': 4},		
                                  {'expression': '"t2_pop_sum"','length': 23,'name': 't2_Pop','precision': 15,'type': 4},		
                                  {'expression': '"t2_bu_sum"','length': 23,'name': 't2_BU_Sum','precision': 15,'type': 4},		
                                  {'expression': '"t2_bua(m²)"','length': 23,'name': 't2_BUA(m²)','precision': 15,'type': 4},		
                                  {'expression': '"area(m)"','length': 23,'name': 'Area(m²)','precision': 15,'type': 4}],		
                'OUTPUT':clean_final_agglomerations_t2})
