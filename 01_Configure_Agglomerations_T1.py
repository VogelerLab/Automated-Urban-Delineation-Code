# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 10:57:28 2022

@author: orioncr

Description: This script aims to automate delineation of urban units using a
plethora of open source datasets and software (QGIS, OpenStreetMaps, openrouteservice,etc).
QGIS is the main asset being utilized to conduct geospatial analyses.
The main input for this script is the urban clusters derived from scipt 00_Derive_Urban_Extent.

This python file is for the initial year (t1)
"""
#%% Import Libraries and Initialize QGIS
import os
import json
import csv
import shutil
import glob
import math
import plotly
import itertools
import fiona
import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd
from functools import reduce
from time import sleep
from pathlib import Path
from osgeo import gdal
from osgeo import ogr
from shapely.geometry import shape
from qgis.core import (
    QgsProject,
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
from qgis.PyQt.QtCore import QVariant
from qgis.utils import iface

#Initialize QGIS
#Supply path to qgis install location
#You can get this using QgsApplication.prefixPath() in QGIS Python console
QgsApplication.setPrefixPath('C:/OSGeo4W/apps/qgis-ltr', True)
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

#Print list of algorithms in registry
# =============================================================================
# for alg in QgsApplication.processingRegistry().algorithms():
#     print(alg.id(), "->", alg.displayName())
# =============================================================================

# Set up OpenRouteService for use in Python
import openrouteservice as ors
import requests

#An API key will need to be obtained from openrouteservice to utilize this script
#Insert your API key below
ors_key = ''
client = ors.Client(ors_key)
#%% Set General Variables
#Enter the abbreviation of your country or region of interest here (as a string)
country_abbrv = ""

#Enter the year of interest here (as a string)
#We default to referring to initial year in the study period as t1
year = "t1"

#Enter code for project as a string
#Example projection: Africa Albers Equal Area
projection = 'ESRI:102022'

#This variable determines the number of seconds a function will be suspended for
#This is necessary when using the OpenRouteServiceAPI which has restrictions
#on the number of requests that can be sent to the API per minute. You may have
#to adjust this number to satisfy the restrictions of the API.
sleep_count = 10

#%% Set Threshold Variables
#Minimum population size threshold to classify a core cluster
core_cluster_threshold = 5000

#Minimum population density for a core threshold
core_density_threshold = 1000

#Minimum population size threshold for a non core cluster
noncore_threshold_size = 500

#Minimum population density threshold for a non core cluster
noncore_threshold_density = 500

#Isochrone distance threshold to associate cores that likely interact but are fragmented in meters
polycore_threshold = 5000

#Isochrone distance threshold to associate non core clusters to cores in meters
nc_to_core_threshold = 5000
#%% Directories And Data Acquisition
#Manually obtain land cover rasters, gridded population rasters, 
#country boundary layers, and locality point data.

dir_00_outputs = r"Path to folder" + country_abbrv + "/" +  year

dir_01 = r"Path to folder" + country_abbrv + "/" +  year

if not os.path.exists(dir_01):
    os.makedirs(dir_01)

dir_01_outputs = r"Path to folder" + country_abbrv + "/" +  year + "/outputs"

if not os.path.exists(dir_01_outputs):
    os.makedirs(dir_01_outputs)

dir_popdata = r"Path to folder"

dir_country_boundaries = r"Path to folder"

dir_places = r"Path to folder"

dir_iso_outputs = r"Path to folder" + country_abbrv + "/" +  year + "/iso"

#The outputs are saved to files for quality assessment purposes.
#The following line of code removes output files from previous attempts
if os.path.exists(dir_iso_outputs):
    shutil.rmtree(dir_iso_outputs)

if not os.path.exists(dir_iso_outputs):
    os.makedirs(dir_iso_outputs)

dir_ncc = os.path.join(dir_01, "ncc")

if not os.path.exists(dir_ncc):
    os.mkdir(dir_ncc)

if not os.path.exists(dir_ncc):
    os.makedirs(dir_ncc)

split_vctr_dir = os.path.join(dir_01, "split_vctr")

if os.path.exists(split_vctr_dir):
    shutil.rmtree(split_vctr_dir)

if not os.path.exists(split_vctr_dir):
    os.mkdir(split_vctr_dir)

pop_raster_dir = os.path.join(dir_01, "pop_rasters")

if os.path.exists(pop_raster_dir):
    shutil.rmtree(pop_raster_dir)

if not os.path.exists(pop_raster_dir):
    os.mkdir(pop_raster_dir)

meancoords_dir = os.path.join(dir_01, "meancoords")

if os.path.exists(meancoords_dir):
    shutil.rmtree(meancoords_dir)

if not os.path.exists(meancoords_dir):
    os.mkdir(meancoords_dir)

nc_mc_dir = os.path.join(dir_01, "nc_mc")

if os.path.exists(nc_mc_dir):
    shutil.rmtree(nc_mc_dir)

if not os.path.exists(nc_mc_dir):
    os.mkdir(nc_mc_dir)
    
nc_pr_dir = os.path.join(dir_01, "nc_pr")

if os.path.exists(nc_pr_dir):
    shutil.rmtree(nc_pr_dir)

if not os.path.exists(nc_pr_dir):
    os.mkdir(nc_pr_dir)
    
nc_iso_dir = os.path.join(dir_01, "nc_iso")

if os.path.exists(nc_iso_dir):
    shutil.rmtree(nc_iso_dir)

if not os.path.exists(nc_iso_dir):
    os.mkdir(nc_iso_dir)
    
#Delete all outputs from previous attempts if desired 
#for f in os.listdir(dir_01_outputs):
#    os.remove(os.path.join(dir_01_outputs, f))

#%% Variables Containing File Paths to Data
#File path for raw country boundary
ctry_bndry_fp = os.path.join(dir_country_boundaries, 'Name of Country Boundary File')

#File path for urban clusters derived from Atlas of Urban Expansion 2016 Urban Extent method
raw_urbanextent_fp = os.path.join(dir_00_outputs, 'Urban clusters name' + country_abbrv + "_" +  year + '.shp')

#File path for raw population gridded raster
raw_pop_grid_fp = os.path.join(dir_popdata, "Gridded pop data name")

#File path for all raw locality data from OSM
raw_city_data_fp = os.path.join(dir_places, 'Name of raw city data')

raw_town_data_fp = os.path.join(dir_places, 'Name of raw town data')
#%% Processing 
#Reproject country boundary and locality data

country_boundary = os.path.join(dir_01_outputs, "Country_Boundary_" + country_abbrv + "_" +  year + ".shp")

processing.run("native:reprojectlayer", {'INPUT':ctry_bndry_fp,
                                         'TARGET_CRS':QgsCoordinateReferenceSystem(projection),
                                         'OPERATION':'+proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad +step +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84',
                                         'OUTPUT':country_boundary})

#reproj_city_data = os.path.join(dir_01_outputs, "ET_OSM_cities.shp")
reproj_city_data = processing.run("native:reprojectlayer", {'INPUT':raw_city_data_fp,
                                         'TARGET_CRS':QgsCoordinateReferenceSystem('ESRI:102022'),
                                         'OPERATION':'+proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad +step +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84',
                                         'OUTPUT':'memory:'})

reproj_city_data = reproj_city_data['OUTPUT']

#reproj_town_data = os.path.join(dir_01_outputs, "ET_OSM_towns.shp")
reproj_town_data = processing.run("native:reprojectlayer", {'INPUT':raw_town_data_fp,
                                         'TARGET_CRS':QgsCoordinateReferenceSystem('ESRI:102022'),
                                         'OPERATION':'+proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad +step +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84',
                                         'OUTPUT':'memory:'})

reproj_town_data = reproj_town_data['OUTPUT']



#Clip city and town data to the country's boundary

city_data_full = os.path.join(dir_01_outputs, country_abbrv + "_OSM_cities_full.shp")

processing.run("native:clip", 
               {'INPUT':reproj_city_data,
                'OVERLAY':country_boundary,
                'OUTPUT':city_data_full})

town_data_full = os.path.join(dir_01_outputs, country_abbrv + "_OSM_towns_full.shp")

processing.run("native:clip", 
               {'INPUT':reproj_town_data,
                'OVERLAY':country_boundary,
                'OUTPUT':town_data_full})
#%% Calculate population and built up characteristics for all clusters
#Load urban extent polygon and gridded pop data as QGS layers

raw_urbanextent = QgsVectorLayer(raw_urbanextent_fp, "", "ogr")
raw_pop_grid = QgsRasterLayer(raw_pop_grid_fp, "")

#Fix geometry of polygon layer
#This may be necessary if switching between GIS programs
fg_ue_t1 = processing.run('native:fixgeometries',
               {'INPUT':raw_urbanextent,
                'OUTPUT':'memory:'})

fg_ue_t1 = fg_ue_t1['OUTPUT']

#Reproject population data
pg_t1_fp = os.path.join(dir_01_outputs, "pop_grid_" + country_abbrv + "_" +  year + ".tif")

processing.run("gdal:warpreproject", 
                              {'INPUT':raw_pop_grid,
                               'SOURCE_CRS':, #Fill value
                               'TARGET_CRS':, #Fill value
                               'RESAMPLING':0,
                               'NODATA':None,
                               'TARGET_RESOLUTION':None,
                               'OPTIONS':'',
                               'DATA_TYPE':0,
                               'TARGET_EXTENT':None,
                               'TARGET_EXTENT_CRS':None,
                               'MULTITHREADING':False,
                               'EXTRA':'',
                               'OUTPUT':pg_t1_fp})

#Clip cluster data to country boundaries
#Subsitute output for variable below if wanting to save permanent file
#urban_clusters_t1 = os.path.join(dir_01_outputs, "ue_clusters_" + country_abbrv + "_" +  year + ".shp")

urban_clusters_t1 = processing.run("native:clip", 
               {'INPUT':fg_ue_t1,
                'OVERLAY':country_boundary,
                'OUTPUT':'memory:'})
                
urban_clusters_t1 = urban_clusters_t1['OUTPUT']

#Zonal statistics to determine population within each urban cluster

ue_with_popdata_t1 = processing.run("native:zonalstatisticsfb", 
               {'INPUT':urban_clusters_t1,
                'INPUT_RASTER':pg_t1_fp,
                'RASTER_BAND':1,
                'COLUMN_PREFIX':'_',
                'STATISTICS':[0,1,2],
                'OUTPUT':'memory:'}) 

ue_with_popdata_t1 = ue_with_popdata_t1['OUTPUT']

# Calculate area field for polygons

ue_area_m = processing.run("native:fieldcalculator",
               {'INPUT':ue_with_popdata_t1,
                'FIELD_NAME':'area(m²)',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'$area',
                'OUTPUT':'memory:'})

ue_area_m_output = ue_area_m['OUTPUT']

ue_area_km = processing.run("native:fieldcalculator",
               {'INPUT':ue_area_m_output,
                'FIELD_NAME':'area(km²)',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'"area(m²)" / 1000000',
                'OUTPUT':'memory:'})

ue_area_km_output = ue_area_km['OUTPUT']

ue_density = processing.run("native:fieldcalculator",
               {'INPUT':ue_area_km_output,
                'FIELD_NAME':'density',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'"_sum"/"area(km²)"',
                'OUTPUT':'memory:'})

ue_density_output = ue_density['OUTPUT']

#This shapefile will contain all of the urban extent clusters with calculations
ue_clusters_with_calcs = os.path.join(dir_01_outputs, "ue_clusters_with_calcs.shp")

ue_rm_null = processing.run("native:extractbyexpression", 
                {'INPUT':ue_density_output,
                 'EXPRESSION':' "_sum" >= 1',
                 'OUTPUT':ue_clusters_with_calcs})

#os.remove(ue_clusters_with_calcs)
#%% Separate core clusters from non-core clusters

#Select polygons with a population value greater than the core threshold

ue_clusters_with_calcs_vctr = QgsVectorLayer(ue_clusters_with_calcs,"", "ogr")

processing.run("qgis:selectbyexpression", 
                      {'INPUT':ue_clusters_with_calcs_vctr,
                       'EXPRESSION':' "_sum" >= {}'.format(core_cluster_threshold),
                       'METHOD':0})

#ue_pop_t1_greater_than_threshold = os.path.join(dir_01_outputs, 'UE_Pop_ET_t1_Greater_Than_Threshold.shp') 
clusters_greater_threshold = processing.run("native:saveselectedfeatures", 
               {'INPUT':ue_clusters_with_calcs_vctr,
                'OUTPUT':'memory:'})

clusters_greater_threshold = clusters_greater_threshold['OUTPUT']
#Select polygons with a population value less than the core threshold

processing.run("qgis:selectbyexpression", 
                      {'INPUT':ue_clusters_with_calcs_vctr,
                       'EXPRESSION':' "_sum" <= {}'.format(core_cluster_threshold),
                       'METHOD':0})

clusters_less_threshold = os.path.join(dir_01_outputs, 'clusters_less_threshold.shp')

processing.run("native:saveselectedfeatures", 
               {'INPUT':ue_clusters_with_calcs_vctr,
                'OUTPUT':clusters_less_threshold})

#Remove non-core clusters with a density less than the threshold we've set
clusters_less_threshold_vctr = QgsVectorLayer(clusters_less_threshold,"", "ogr")

processing.run("qgis:selectbyexpression", 
                      {'INPUT':clusters_less_threshold_vctr,
                       'EXPRESSION':' "_sum" >= {} AND "density" >= {}'.format(noncore_threshold_size, noncore_threshold_density),
                       'METHOD':0})

non_core_clusters = os.path.join(dir_01_outputs, 'non_core_clusters.shp') 

processing.run("native:saveselectedfeatures", 
               {'INPUT':clusters_less_threshold_vctr,
                'OUTPUT':non_core_clusters})

#os.remove(non_core_clusters)

#%% Clean place data
city_data = os.path.join(dir_01_outputs, country_abbrv + "_OSM_cities.shp")

processing.run("native:retainfields",
               {'INPUT':city_data_full,
                'FIELDS':['full_id','osm_id','place','name_en','name','alt_name'],
                'OUTPUT':city_data})
 
town_data = os.path.join(dir_01_outputs, country_abbrv + "_OSM_towns.shp")

processing.run("native:retainfields",
               {'INPUT':town_data_full,
                'FIELDS':['full_id','osm_id','place','name_en','name','alt_name'],
                'OUTPUT':town_data})

# =============================================================================
# # Generate a layer of merged locality vector layers if wanted
# all_locality_data = os.path.join(dir_01_outputs, "ET_OSM_all_localities.shp")
# 
# processing.run("native:mergevectorlayers", 
#                {'LAYERS':[city_data,town_data],
#                 'CRS':QgsCoordinateReferenceSystem('ESRI:102022'),
#                 'OUTPUT':all_locality_data})
# =============================================================================

#%% Associate localities with cores
# Cities are joined first

city_within_core = processing.run("native:joinbynearest", 
               {'INPUT':clusters_greater_threshold,
                'INPUT_2':city_data,
                'FIELDS_TO_COPY':[],
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'NEIGHBORS':1,
                'MAX_DISTANCE':5000,
                'OUTPUT':'memory:'})

city_within_core_output = city_within_core['OUTPUT']

#city_cores = os.path.join(dir_01_outputs, 'city_cores.shp')
city_cores = processing.run("native:extractbyexpression", {'INPUT':city_within_core_output,
                                              'EXPRESSION':' "full_id" IS NOT NULL',
                                              'OUTPUT':'memory:'})

city_cores = city_cores['OUTPUT']

ncc = processing.run("native:extractbyexpression", {'INPUT':city_within_core_output,
                                              'EXPRESSION':' "full_id" IS NULL',
                                              'OUTPUT':'memory:'})

ncc_output = ncc['OUTPUT']

cores_noid = processing.run("native:retainfields", 
               {'INPUT':ncc_output,
                'FIELDS':['Id','gridcode','_count','_sum','_mean','area(m²)','area(km²)','density'],
                'OUTPUT':'memory:'})

cores_noid_output = cores_noid['OUTPUT']

# Now we join towns

town_within_core = processing.run("native:joinbynearest", 
               {'INPUT':cores_noid_output,
                'INPUT_2':town_data,
                'FIELDS_TO_COPY':[],
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'NEIGHBORS':1,
                'MAX_DISTANCE':5000,
                'OUTPUT':'memory:'})

town_within_core_output = town_within_core['OUTPUT']

#town_cores = os.path.join(dir_01_outputs, 't1_Town_Cores.shp')

town_cores = processing.run("native:extractbyexpression", {'INPUT':town_within_core_output,
                                              'EXPRESSION':' "full_id" IS NOT NULL',
                                              'OUTPUT':'memory:'})

town_cores = town_cores['OUTPUT']

ntc = processing.run("native:extractbyexpression", {'INPUT':town_within_core_output,
                                              'EXPRESSION':' "full_id" IS NULL',
                                              'OUTPUT':'memory:'})

ntc_output = ntc['OUTPUT']

cores_noid = processing.run("native:retainfields", 
               {'INPUT':ntc_output,
                'FIELDS':['Id','gridcode','_count','_sum','_mean','area(m²)','area(km²)','density'],
                'OUTPUT':'memory:'})

cores_noid_output = cores_noid['OUTPUT']

core_with_locality = os.path.join(dir_01_outputs, "core_with_locality.shp")

processing.run("native:mergevectorlayers", 
               {'LAYERS':[city_cores, town_cores],
                'CRS':QgsCoordinateReferenceSystem('ESRI:102022'),
                'OUTPUT':core_with_locality})

combine_cores = processing.run("native:mergevectorlayers", 
                       {'LAYERS':[core_with_locality, cores_noid_output],
                        'CRS':None,
                        'OUTPUT':'memory:'})
 
combine_cores= combine_cores['OUTPUT']

#Delete duplicates. Some cores have multiple localities located within their polygon.
# Since the distance of the localities from the core are the same, the join by nearest tool
# automatically creates multiple outptuts for these cores. We have to remove these
# in order to contineu. The following is not the best method for getting around duplications
# but it works for now. This is just to associate localities. We can also join 
# all of the localities associated with a single agglomeration at the end of our analysis

rm_duplicates = processing.run("native:removeduplicatesbyattribute", 
               {'INPUT':combine_cores,
                'FIELDS':['Id'],
                'OUTPUT':'memory:'})

rm_duplicates = rm_duplicates ['OUTPUT']

all_cores = os.path.join(dir_01_outputs, "all_cores.shp")

processing.run("native:addautoincrementalfield", {'INPUT':rm_duplicates,
                                                  'FIELD_NAME':'Unique_ID',
                                                  'START':0,
                                                  'MODULUS':0,
                                                  'GROUP_FIELDS':[],                 
                                                  'SORT_EXPRESSION':'',
                                                  'SORT_ASCENDING':True,
                                                  'SORT_NULLS_FIRST':False,
                                                  'OUTPUT':all_cores})


#%%
# Extract by mask and split population raster for each polygon
all_cores_vctr = QgsVectorLayer(all_cores, '', 'ogr')

# Iterate over features in shapefile
# Split core layer by ID attribute
split_vctr = processing.run("native:splitvectorlayer", 
               {'INPUT':all_cores_vctr,
                'FIELD':'Unique_ID',
                'FILE_TYPE':0,
                'OUTPUT':split_vctr_dir})

# Get list of filepaths for each shapefile in directory
for filename in os.listdir(split_vctr_dir):
    # Index filename to pull ID to assign to raster
    # print(filename[7:-5])
    
    f = os.path.join(split_vctr_dir, filename)

    all_core_rasters = os.path.join(pop_raster_dir, str(filename[7:-5]) + "_raster.tif")

    processing.run("gdal:cliprasterbymasklayer",
                       {'INPUT':pg_t1_fp,
                        'MASK':f,
                        'SOURCE_CRS':None,
                        'TARGET_CRS':None,
                        'NODATA':None,
                        'ALPHA_BAND':False,
                        'CROP_TO_CUTLINE':True,
                        'KEEP_RESOLUTION':True,
                        'SET_RESOLUTION':False,
                        'X_RESOLUTION':None,
                        'Y_RESOLUTION':None,
                        'MULTITHREADING':False,
                        'OPTIONS':'',
                        'DATA_TYPE':0,
                        'EXTRA':'',
                        'OUTPUT':all_core_rasters})

# Convert each raster to a unique point data set
os.chdir(pop_raster_dir)

for lyr in glob.glob("*.tif"):
    
    rstr_to_points = processing.run("native:pixelstopoints", 
                   {'INPUT_RASTER':lyr,
                    'RASTER_BAND':1,
                    'FIELD_NAME':'VALUE',
                    'OUTPUT':'memory:'})
    
    rstr_to_points = rstr_to_points['OUTPUT']
    
    meancoords = meancoords_dir + '/' + 'mc_' + lyr[:-11] + '.shp'

    meancoords_mmry = processing.run("native:meancoordinates", 
                   {'INPUT':rstr_to_points,
                    'WEIGHT':'VALUE',
                    'UID':'fid',
                    'OUTPUT':meancoords})

# Merge memory layers containing mean coordinates for each individual feature to
# create a single layer with all mean coordinates

meancoords_dp = os.path.join(meancoords_dir, '*.shp')
meancoordfiles = glob.glob(meancoords_dp)
print(meancoordfiles)

meancoords_process = processing.run("native:mergevectorlayers", 
               {'LAYERS':meancoordfiles,
                'CRS':None,
                'OUTPUT':'memory:'})

meancoords_output = meancoords_process['OUTPUT']

# Join core information to population weighted centroids. This information will
# be of use when creating isochrones for each core. Use field calculator to add
# a field with the Unique ID which is stored in the layer  field.

#Add field with unique ID to mean coord points
add_field = processing.run("native:fieldcalculator", 
               {'INPUT':meancoords_output,
                'FIELD_NAME':'Unique_ID',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'right(  "layer",  length( "layer" ) - 6)',
                'OUTPUT':'memory:'})

add_field = add_field['OUTPUT']

#Join attributes
pop_weighted_ctrds = os.path.join(dir_01_outputs, "pop_weighted_ctrds.shp")

processing.run("native:joinattributestable",
               {'INPUT':add_field,
                'FIELD':'Unique_ID',
                'INPUT_2':all_cores,
                'FIELD_2':'Unique_ID',
                'FIELDS_TO_COPY':[],
                'METHOD':0,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'OUTPUT':pop_weighted_ctrds})



#%% Prepare Data For Travel Distance Matrix
# Create convex hull for core polygons

cvx_hull_cores = processing.run("native:convexhull", 
                                {'INPUT':all_cores,
                                 'OUTPUT':'memory:'})

cvx_hull_cores = cvx_hull_cores['OUTPUT']

# Extract vertices from convex hulls of core polygons
core_vrtcs = processing.run("native:extractvertices", 
               {'INPUT': cvx_hull_cores,
                'OUTPUT':'memory:'})

core_vrtcs = core_vrtcs['OUTPUT']

# Reproject centroids and vertices to WGS4
reproj_ctrds = processing.run("native:reprojectlayer", 
                              {'INPUT':pop_weighted_ctrds,
                               'TARGET_CRS':QgsCoordinateReferenceSystem('EPSG:4326'),
                               'OPERATION':'+proj=pipeline +step +inv +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84 +step +proj=unitconvert +xy_in=rad +xy_out=deg',
                               'OUTPUT':'memory:'})

reproj_ctrds = reproj_ctrds['OUTPUT']

reproj_vrtcs = processing.run("native:reprojectlayer", 
                              {'INPUT':core_vrtcs,
                               'TARGET_CRS':QgsCoordinateReferenceSystem('EPSG:4326'),
                               'OPERATION':'+proj=pipeline +step +inv +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84 +step +proj=unitconvert +xy_in=rad +xy_out=deg',
                               'OUTPUT':'memory:'})

reproj_vrtcs = reproj_vrtcs['OUTPUT']

# Add geometry
ctrds_geo = processing.run("qgis:exportaddgeometrycolumns", 
               {'INPUT':reproj_ctrds,
                'CALC_METHOD':0,
                'OUTPUT':'memory:'})

ctrds_geo = ctrds_geo['OUTPUT']

vrtcs_geo = processing.run("qgis:exportaddgeometrycolumns", 
               {'INPUT':reproj_vrtcs,
                'CALC_METHOD':0,
                'OUTPUT':'memory:'})

vrtcs_geo = vrtcs_geo['OUTPUT']

# Export attribute table as CSV
ctrds_csv = os.path.join(dir_01_outputs, "ctrds_csv.xlsx")

vrtcs_csv = os.path.join(dir_01_outputs, "vrtcs_csv.xlsx")

processing.run("native:exporttospreadsheet", 
               {'LAYERS':ctrds_geo,
                'USE_ALIAS':False,
                'FORMATTED_VALUES':False,
                'OUTPUT':ctrds_csv,
                'OVERWRITE':True})

processing.run("native:exporttospreadsheet", 
               {'LAYERS':vrtcs_geo,
                'USE_ALIAS':False,
                'FORMATTED_VALUES':False,
                'OUTPUT':vrtcs_csv,
                'OVERWRITE':True})

# Convert CSVs to DataFrames
ctrds_df = pd.read_excel(ctrds_csv)
vrtcs_df = pd.read_excel(vrtcs_csv)

# Create location column containing a list of x and y coords
merge_coords_vc = vrtcs_df.assign(location_vc=vrtcs_df[['xcoord', 'ycoord']].values.tolist())

#Create a nested list containing all x and y coords for the associated ID
vc_coords = merge_coords_vc.groupby('Id')['location_vc'].apply(list)

#Join dataframes
merge_dfs = pd.merge(ctrds_df, vc_coords, on='Id', how='left')

#Merge x and y columns in centroids dataframe
merge_coords_cd = merge_dfs.assign(location_cd=merge_dfs[['xcoord', 'ycoord']].values.tolist())

#Extract necessary columns from larger dataframe
vc_cd_df = merge_coords_cd.filter(['Id', 'location_cd', 'location_vc'])
    
#Concatenate coordinate lists. Centroid coordinates should be first in the new list
vc_cd_df['locations'] = vc_cd_df.apply(lambda x: [x['location_cd']] + x['location_vc'], axis=1)

vc_cd_df = vc_cd_df.filter(['Id', 'locations'])

vc_cd_csv = os.path.join(dir_01_outputs, "vc_cd_csv.csv")

vc_cd_df.to_csv(vc_cd_csv)
#%%Use ORS API to create travel distance matrix for each polygon
column_names = ['Id', 'avg_dist']

dm_df = pd.DataFrame(columns = column_names)

for index, row in vc_cd_df.iterrows():
    
    length = (len(row.locations))
# =============================================================================
#     print(length)
#     print(list(range(1,length)))
#     print(row.locations[1:length])
# =============================================================================
    body = {"locations":row.locations,
            "destinations":list(range(1,length)),
            "id": row.Id,
            "metrics":["distance"],
            "sources":[0],
            "units":"m"}

    headers = {
        'Accept': 'application/json, application/geo+json, application/gpx+xml, img/png; charset=utf-8',
        'Authorization': ors_key,
        'Content-Type': 'application/json; charset=utf-8'
        }

    call = requests.post('https://api.openrouteservice.org/v2/matrix/driving-car', json=body, headers=headers)

    call_stat_and_reason = (call.status_code, call.reason)

    geojson = call.text
    
    js = json.loads(geojson)
    
    #Append average distance to dataframe
    dist = js['distances']
    dist_list = reduce(lambda x,y: x+y, dist)
    cleaned_list = []
    for val in dist_list:
        if val != None :
            cleaned_list.append(val)
    
    sum_list = sum(cleaned_list)
    len_list = len(cleaned_list)
    if sum_list == 0 or len_list == 0:
        continue
    avg_dist = sum_list / len_list
    pt_id = js['metadata']['id']
    appnd_df = pd.DataFrame({'Id':[pt_id],
                    'avg_dist':[avg_dist]})
   
    dm_df = pd.concat([appnd_df, dm_df], ignore_index=True)
    
    sleep(sleep_count)
    
# Save dataframe as CSV    
dm_csv = os.path.join(dir_01_outputs, 'distance_matrix.csv')
dm_df.to_csv(dm_csv)

# Convert id to int
#id_to_int = os.path.join(dir_01_outputs, "id_to_int.dbf")
id_to_int = processing.run("native:fieldcalculator", 
               {'INPUT':dm_csv,
                'FIELD_NAME':'id_int',
                'FIELD_TYPE':1,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':' to_int("Id") ',
                'OUTPUT':'memory:'})

id_to_int = id_to_int['OUTPUT']

# join CSV containing distance values to population weighted centroids layer
dist_output = processing.run("native:joinattributestable", 
               {'INPUT':pop_weighted_ctrds,
                'FIELD':'Id',
                'INPUT_2':id_to_int,
                'FIELD_2':'id_int',
                'FIELDS_TO_COPY':['avg_dist'],
                'METHOD':1,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'OUTPUT':'memory:'})

dist_output = dist_output['OUTPUT']

#%%
# Reproject population weighted centroids to get lat and long
pop_weighted_ctrds_WGS84 = processing.run("native:reprojectlayer", 
               {'INPUT':dist_output,
                'TARGET_CRS':QgsCoordinateReferenceSystem('EPSG:4326'),
                'OPERATION':'+proj=pipeline +step +inv +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84 +step +proj=unitconvert +xy_in=rad +xy_out=deg',
                'OUTPUT':'memory:'})

pop_weighted_ctrds_WGS84 = pop_weighted_ctrds_WGS84['OUTPUT']

# Add lat and long to population weighted centroids
add_geo = processing.run("qgis:exportaddgeometrycolumns", 
               {'INPUT':pop_weighted_ctrds_WGS84,
                'CALC_METHOD':0,
                'OUTPUT':'memory:'})

add_geo = add_geo['OUTPUT']

# Reproject centroids layer back to appropriate coordinate system
repr_ae = processing.run("native:reprojectlayer", 
                            {'INPUT':add_geo,
                             'TARGET_CRS':, #Fill value
                             'OPERATION':'+proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad +step +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84',
                             'OUTPUT':'memory:'})

repr_ae = repr_ae['OUTPUT']

# Create final core
rm_null = processing.run("native:extractbyexpression", 
               {'INPUT':repr_ae,
                'EXPRESSION':' "avg_dist" IS NOT NULL',
                'OUTPUT':'memory:'})

rm_null = rm_null['OUTPUT']

core_pw_centroids = os.path.join(dir_01_outputs, "core_pw_centroids.shp")

processing.run("native:fieldcalculator", 
               {'INPUT':rm_null,
                'FIELD_NAME':'iso_dist',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':' "avg_dist" + ' + str(polycore_threshold),
                'OUTPUT':core_pw_centroids})

core_pw_ctrds_table = os.path.join(dir_01_outputs, "core_pw_ctrds_table.xlsx")

processing.run("native:exporttospreadsheet", 
               {'LAYERS':core_pw_centroids,
                'USE_ALIAS':False,
                'FORMATTED_VALUES':False,
                'OUTPUT':core_pw_ctrds_table,
                'OVERWRITE':True})

#%% Isochrone maps to associate cores
core_pw_centroids_df = pd.read_excel(core_pw_ctrds_table)
core_pw_centroids_df = core_pw_centroids_df.reset_index()
core_pw_centroids_df = core_pw_centroids_df[core_pw_centroids_df.iso_dist <= 120000] #Isochrone distance cannot exceed 120,000 meters

for f in os.listdir(dir_iso_outputs):
       os.remove(os.path.join(dir_iso_outputs, f))

for index, row in core_pw_centroids_df.iterrows():
    
    body = {"locations":[[(row["xcoord"]),
                          (row["ycoord"])]],
            "range":[row["iso_dist"]],
            "id":str((row["Id"])),
            "range_type":"distance",
            "units":"m"}

    headers = {
        'Accept': 'application/json, application/geo+json, application/gpx+xml, img/png; charset=utf-8',
       'Authorization': ors_key,
       'Content-Type': 'application/json; charset=utf-8'
       }
    
    call = requests.post('https://api.openrouteservice.org/v2/isochrones/driving-car', json=body, headers=headers)

    call_stat_and_reason = (call.status_code, call.reason)

    geojson = call.text

    gdf = gpd.read_file(geojson)

    isochrone_x = os.path.join(dir_iso_outputs, "iso_" + str(row["Id"]) + ".shp")
    
    gdf.to_file(isochrone_x)
    
    sleep(sleep_count)
    
# Merge vector layers
# Obtain only shapefiles from our isochrone directory folder
shp_iso = os.path.join(dir_iso_outputs, '*.shp')
shp_files = glob.glob(shp_iso)
print(shp_files)

#iso_merge = os.path.join(dir_01_outputs, "isochrones_merged.shp")
iso_merge = processing.run("native:mergevectorlayers",
               {'LAYERS': shp_files,
               'CRS':, # Fill value
               'OUTPUT':'memory:'})
    
iso_merge = iso_merge['OUTPUT']
    
# Add iso_id field
isochrones = os.path.join(dir_01_outputs, "core_isochrones.shp")

processing.run("native:fieldcalculator", 
               {'INPUT':iso_merge,
                'FIELD_NAME':'iso_id',
                'FIELD_TYPE':2,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'substr("layer",regexp_match("layer",\'iso_\')+4)',
                'OUTPUT':isochrones})
#%% Determine polycores

# To connect cores, join attributes by location. This would
# give us the core attributes for potentially multiple polygons.
#core_with_iso = os.path.join(dir_01_outputs, "core_test.shp")
core_with_iso = processing.run("native:joinattributesbylocation", 
               {'INPUT':all_cores,
                'JOIN':isochrones,
                'PREDICATE':[0,1,2,3,4,5,6],
                'JOIN_FIELDS':['iso_id'],
                'METHOD':0,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'OUTPUT':'memory:'})

core_with_iso = core_with_iso['OUTPUT']

# Dissolve overlapping features
#dissolve_iso_id = os.path.join(dir_01_outputs, "disslv_iso_id.shp")
dissolve_iso_id = processing.run("native:dissolve", 
               {'INPUT':core_with_iso,
                'FIELD':['iso_id'],
                'OUTPUT':'memory:'})

dissolve_iso_id = dissolve_iso_id['OUTPUT']

# Dissolve by Id
#dissolve_id = os.path.join(dir_01_outputs, "disslv_id.shp")
dissolve_id = processing.run("native:dissolve", 
               {'INPUT':dissolve_iso_id,
                'FIELD':['Id'],
                'OUTPUT':'memory:'})

dissolve_id = dissolve_id['OUTPUT']

# Join attributes of multipart polygons with overlapping parts
#join_overlap = os.path.join(dir_01_outputs, "join_overlap.shp")
join_overlap = processing.run("native:joinattributesbylocation", 
               {'INPUT':dissolve_id,
                'JOIN':dissolve_id,
                'PREDICATE':[0,1,2,3,4,5,6],
                'JOIN_FIELDS':[],
                'METHOD':0,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'OUTPUT':'memory:'})

join_overlap = join_overlap['OUTPUT']

# Dissolve by iso_Id_2
#iso_id_2= os.path.join(dir_01_outputs, "iso_id_2.shp")
iso_id_2 = processing.run("native:dissolve", 
               {'INPUT':join_overlap,
                'FIELD':['iso_id_2'],
                'OUTPUT':'memory:'})

iso_id_2 = iso_id_2['OUTPUT']

# Dissolve by Id once more
id_dissolve = processing.run("native:dissolve", 
               {'INPUT':iso_id_2,
                'FIELD':['Id'],
                'OUTPUT':'memory:'})

id_dissolve = id_dissolve['OUTPUT']

# Give new unique ID
new_id = processing.run("native:fieldcalculator", 
               {'INPUT':id_dissolve,
                'FIELD_NAME':'new_id',
                'FIELD_TYPE':1,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'$Id',
                'OUTPUT':'memory:'})

new_id = new_id['OUTPUT']

#Only retain new_id field
new_cores = os.path.join(dir_01_outputs, "cores.shp")

processing.run("native:retainfields", {'INPUT':new_id,
                                       'FIELDS':['new_id'],
                                       'OUTPUT':new_cores})
#%%Join non-core clusters to core clusters
# get centroids for non core clusters
# Extract by mask and split population raster for each polygon
non_core_clusters = QgsVectorLayer(non_core_clusters, '', 'ogr')

# Iterate over features in shapefile

# Split non-core cluster layer by ID attribute
split_ncc = processing.run("native:splitvectorlayer", 
               {'INPUT':non_core_clusters,
                'FIELD':'Id',
                'FILE_TYPE':0,
                'OUTPUT':dir_ncc})
    
# Get list of filepaths for each shapefile in directory and clip raster by each
for filename in os.listdir(dir_ncc):
    # Index filename to pull ID to assign to raster
    # print(filename[3:-5])
    
    f = os.path.join(dir_ncc, filename)

    ncc_rasters = os.path.join(nc_pr_dir, str(filename[3:-5]) + "_raster.tif")

    processing.run("gdal:cliprasterbymasklayer",
                       {'INPUT':pg_t1_fp,
                        'MASK':f,
                        'SOURCE_CRS':None,
                        'TARGET_CRS':None,
                        'NODATA':None,
                        'ALPHA_BAND':False,
                        'CROP_TO_CUTLINE':True,
                        'KEEP_RESOLUTION':True,
                        'SET_RESOLUTION':False,
                        'X_RESOLUTION':None,
                        'Y_RESOLUTION':None,
                        'MULTITHREADING':False,
                        'OPTIONS':'',
                        'DATA_TYPE':0,
                        'EXTRA':'',
                        'OUTPUT':ncc_rasters})

# Convert each raster to a unique point data set
os.chdir(nc_pr_dir)

for lyr in glob.glob("*.tif"):
    print(lyr)
    
    rstr_to_points = processing.run("native:pixelstopoints", 
                   {'INPUT_RASTER':lyr,
                    'RASTER_BAND':1,
                    'FIELD_NAME':'VALUE',
                    'OUTPUT':'memory:'})
    
    rstr_to_points = rstr_to_points['OUTPUT']
    
    meancoords = nc_mc_dir + '/' + 'mc_' + lyr[:-11] + '.shp'

    processing.run("native:meancoordinates", 
                   {'INPUT':rstr_to_points,
                    'WEIGHT':'VALUE',
                    'UID':'fid',
                    'OUTPUT':meancoords})    
    
# Merge layers containing mean coordinates for each individual feature to
# create a single layer with all mean coordinates

meancoords_dp = os.path.join(nc_mc_dir, '*.shp')
meancoordfiles = glob.glob(meancoords_dp)
print(meancoordfiles)

#meancoords_output = os.path.join(dir_01_outputs, 'meancoords.shp')
meancoords_process = processing.run("native:mergevectorlayers", 
               {'LAYERS':meancoordfiles,
                'CRS':None,
                'OUTPUT':'memory:'})

meancoords_output = meancoords_process['OUTPUT']

# Join core information to population weighted centroids. This information will
# be of use when creating isochrones for each core. Use field calculator to add
# a field with the Unique ID which is stored in the layer  field.

#Add field with unique ID to mean coord points
#add_field = os.path.join(dir_01_outputs, "add_field.shp")

add_field = processing.run("native:fieldcalculator", 
               {'INPUT':meancoords_output,
                'FIELD_NAME':'Unique_ID',
                'FIELD_TYPE':1,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'right(  "layer",  length( "layer" ) -3)',
                'OUTPUT':'memory:'})

add_field = add_field['OUTPUT']

#Join attributes
ncc_pwc = os.path.join(dir_01_outputs, "ncc_pwc.shp")

processing.run("native:joinattributestable",
               {'INPUT':add_field,
                'FIELD':'Unique_ID',
                'INPUT_2':non_core_clusters,
                'FIELD_2':'Id',
                'FIELDS_TO_COPY':[],
                'METHOD':0,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'OUTPUT':ncc_pwc})

#%% Create travel distance matrix
# Create convex hull for core polygons

cvx_hull_nc = processing.run("native:convexhull", 
                                {'INPUT':non_core_clusters,
                                 'OUTPUT':'memory:'})

cvx_hull_nc = cvx_hull_nc['OUTPUT']

# Extract vertices from ncc polygons
core_vrtcs = processing.run("native:extractvertices", 
               {'INPUT': cvx_hull_nc,
                'OUTPUT':'memory:'})

core_vrtcs = core_vrtcs['OUTPUT']

# Reproject centroids and vertices to WGS4
reproj_ctrds = processing.run("native:reprojectlayer", 
                              {'INPUT':ncc_pwc,
                               'TARGET_CRS':QgsCoordinateReferenceSystem('EPSG:4326'),
                               'OPERATION':'+proj=pipeline +step +inv +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84 +step +proj=unitconvert +xy_in=rad +xy_out=deg',
                               'OUTPUT':'memory:'})

reproj_ctrds = reproj_ctrds['OUTPUT']

reproj_vrtcs = processing.run("native:reprojectlayer", 
                              {'INPUT':core_vrtcs,
                               'TARGET_CRS':QgsCoordinateReferenceSystem('EPSG:4326'),
                               'OPERATION':'+proj=pipeline +step +inv +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84 +step +proj=unitconvert +xy_in=rad +xy_out=deg',
                               'OUTPUT':'memory:'})

reproj_vrtcs = reproj_vrtcs['OUTPUT']

# Add geometry
ctrds_geo = processing.run("qgis:exportaddgeometrycolumns", 
               {'INPUT':reproj_ctrds,
                'CALC_METHOD':0,
                'OUTPUT':'memory:'})

ctrds_geo = ctrds_geo['OUTPUT']

vrtcs_geo = processing.run("qgis:exportaddgeometrycolumns", 
               {'INPUT':reproj_vrtcs,
                'CALC_METHOD':0,
                'OUTPUT':'memory:'})

vrtcs_geo = vrtcs_geo['OUTPUT']

# Export attribute table as CSV
ctrds_csv = os.path.join(dir_01_outputs, "ncc_ctrds_csv.xlsx")

vrtcs_csv = os.path.join(dir_01_outputs, "ncc_vrtcs_csv.xlsx")

processing.run("native:exporttospreadsheet", 
               {'LAYERS':ctrds_geo,
                'USE_ALIAS':False,
                'FORMATTED_VALUES':False,
                'OUTPUT':ctrds_csv,
                'OVERWRITE':True})

processing.run("native:exporttospreadsheet", 
               {'LAYERS':vrtcs_geo,
                'USE_ALIAS':False,
                'FORMATTED_VALUES':False,
                'OUTPUT':vrtcs_csv,
                'OVERWRITE':True})

# Convert CSV to DataFrame
ctrds_df = pd.read_excel(ctrds_csv)

vrtcs_df = pd.read_excel(vrtcs_csv)

# Create location column containing a list of x and y coords
merge_coords_vc = vrtcs_df.assign(location_vc=vrtcs_df[['xcoord', 'ycoord']].values.tolist())

#Create a nested list containing all x and y coords for the associated ID
vc_coords = merge_coords_vc.groupby('Id')['location_vc'].apply(list)

#Join dataframes
merge_dfs = pd.merge(ctrds_df, vc_coords, on='Id', how='left')

#Merge x and y columns in centroids dataframe
merge_coords_cd = merge_dfs.assign(location_cd=merge_dfs[['xcoord', 'ycoord']].values.tolist())

#Extract necessary columns from larger dataframe
vc_cd_df = merge_coords_cd.filter(['Id', 'location_cd', 'location_vc'])

# Remove nan if applicable
vc_cd_df = pd.DataFrame.dropna(vc_cd_df)
    
#Concatenate coordinate lists. Centroid coordinates should be first in the new list
vc_cd_df['locations'] = vc_cd_df.apply(lambda x: [x['location_cd']] + x['location_vc'], axis=1)

vc_cd_df = vc_cd_df.filter(['Id', 'locations'])
#%% Use ORS API to create travel distance matrix for each polygon
column_names = ['Id', 'avg_dist']

dm_df = pd.DataFrame(columns = column_names)

for index, row in vc_cd_df.iterrows():
    
    length = (len(row.locations))
# =============================================================================
#     print(length)
#     print(list(range(1,length)))
#     print(row.locations[1:length])
# =============================================================================
    body = {"locations":row.locations,
            "destinations":list(range(1,length)),
            "id": row.Id,
            "metrics":["distance"],
            "sources":[0],
            "units":"m"}

    headers = {
        'Accept': 'application/json, application/geo+json, application/gpx+xml, img/png; charset=utf-8',
        'Authorization': ors_key,
        'Content-Type': 'application/json; charset=utf-8'
        }

    call = requests.post('https://api.openrouteservice.org/v2/matrix/driving-car', json=body, headers=headers)

    call_stat_and_reason = (call.status_code, call.reason)

    geojson = call.text
    
    js = json.loads(geojson)
    
    #Append average distance to dataframe
    dist = js['distances']
    dist_list = reduce(lambda x,y: x+y, dist)
    cleaned_list = []
    for val in dist_list:
        if val != None :
            cleaned_list.append(val)
    sum_list = sum(cleaned_list)
    len_list = len(cleaned_list)
    if sum_list == 0 or len_list == 0:
        continue
    avg_dist = sum_list / len_list
    pt_id = js['metadata']['id']
    appnd_df = pd.DataFrame({'Id':[pt_id],
                    'avg_dist':[avg_dist]})
   
    dm_df = pd.concat([appnd_df, dm_df], ignore_index=True)
    
    sleep(sleep_count)
#%%    
# Save dataframe as CSV    
dm_csv = os.path.join(dir_01_outputs, 'ncc_distance_matrix.csv')
    
dm_df.to_csv(dm_csv)

# Convert idf field from string to int
id_to_int = processing.run("native:fieldcalculator", 
               {'INPUT':dm_csv,
                'FIELD_NAME':'id_int',
                'FIELD_TYPE':1,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'to_int(to_real("Id"))',
                'OUTPUT':'memory:'})

id_to_int = id_to_int['OUTPUT']

# join CSV containing distance values to population weighted centroids layer
dist_output = processing.run("native:joinattributestable", 
               {'INPUT':ncc_pwc,
                'FIELD':'Id',
                'INPUT_2':id_to_int,
                'FIELD_2':'id_int',
                'FIELDS_TO_COPY':['avg_dist'],
                'METHOD':1,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'OUTPUT':'memory:'})

dist_output = dist_output['OUTPUT']

#%%
# Reproject population weighted centroids to get lat and long
nc_pwc_ctrds_WGS84 = processing.run("native:reprojectlayer", 
               {'INPUT':dist_output,
                'TARGET_CRS':QgsCoordinateReferenceSystem('EPSG:4326'),
                'OPERATION':'+proj=pipeline +step +inv +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84 +step +proj=unitconvert +xy_in=rad +xy_out=deg',
                'OUTPUT':'memory:'})

nc_pwc_ctrds_WGS84 = nc_pwc_ctrds_WGS84['OUTPUT']

# Add lat and long to population weighted centroids
add_geo = processing.run("qgis:exportaddgeometrycolumns", 
               {'INPUT':nc_pwc_ctrds_WGS84,
                'CALC_METHOD':0,
                'OUTPUT':'memory:'})

add_geo = add_geo['OUTPUT']

# Reproject centroids layer back to appropriate coordinate system
repr_ae = processing.run("native:reprojectlayer", 
                            {'INPUT':add_geo,
                             'TARGET_CRS':, # Fill value
                             'OPERATION':'+proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad +step +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84',
                             'OUTPUT':'memory:'})

repr_ae = repr_ae['OUTPUT']

# Remove nulls
# We may not be able to calculate a distance matrix for some clusters as they do
# not have travel data thereby creating nulls
rm_null = processing.run("native:extractbyexpression",
               {'INPUT':repr_ae,
                'EXPRESSION':' "avg_dist" IS NOT NULL',
                'OUTPUT':'memory:'})

rm_null = rm_null['OUTPUT']

# Create final core
nc_pw_centroids = os.path.join(dir_01_outputs, "nc_pw_centroids.shp")

processing.run("native:fieldcalculator", 
               {'INPUT':rm_null,
                'FIELD_NAME':'iso_dist',
                'FIELD_TYPE':0,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':' "avg_dist" + ' + str(nc_to_core_threshold),
                'OUTPUT':nc_pw_centroids})

nc_pw_ctrds_table = os.path.join(dir_01_outputs, "nc_pw_ctrds_table.xlsx")

processing.run("native:exporttospreadsheet", 
               {'LAYERS':nc_pw_centroids,
                'USE_ALIAS':False,
                'FORMATTED_VALUES':False,
                'OUTPUT':nc_pw_ctrds_table,
                'OVERWRITE':True})

#%% Isochrone maps to associate cores

nc_pw_centroids_df = pd.read_excel(nc_pw_ctrds_table)
nc_pw_centroids_df = nc_pw_centroids_df.reset_index()
nc_pw_centroids_df = nc_pw_centroids_df[nc_pw_centroids_df.iso_dist <= 120000]

for index, row in nc_pw_centroids_df.iterrows():
    
    body = {"locations":[[(row["xcoord"]),
                          (row["ycoord"])]],
            "range":[row["iso_dist"]],
            "id":str((row["Id"])),
            "range_type":"distance",
            "units":"m"}

    headers = {
        'Accept': 'application/json, application/geo+json, application/gpx+xml, img/png; charset=utf-8',
       'Authorization': ors_key,
       'Content-Type': 'application/json; charset=utf-8'
       }
    
    call = requests.post('https://api.openrouteservice.org/v2/isochrones/driving-car', json=body, headers=headers)

    call_stat_and_reason = (call.status_code, call.reason)

    geojson = call.text

    gdf = gpd.read_file(geojson)

    isochrone_x = os.path.join(nc_iso_dir, "iso_" + str(row["Id"]) + ".shp")
    
    gdf.to_file(isochrone_x)
    
    sleep(sleep_count)
    
# Merge vector layers
# Obtain only shapefiles from our isochrone directory folder
shp_iso = os.path.join(nc_iso_dir, '*.shp')
shp_files = glob.glob(shp_iso)
print(shp_files)

iso_merge = processing.run("native:mergevectorlayers",
               {'LAYERS': shp_files,
               'CRS':, # Fill value
               'OUTPUT':'memory:'})
    
iso_merge =iso_merge['OUTPUT']
    
# Add nci_id field 
isochrones = os.path.join(dir_01_outputs, "non_core_isochrones.shp")

processing.run("native:fieldcalculator", 
               {'INPUT':iso_merge,
                'FIELD_NAME':'nci_id',
                'FIELD_TYPE':1,
                'FIELD_LENGTH':0,
                'FIELD_PRECISION':0,
                'FORMULA':'substr("layer",regexp_match("layer",\'iso_\')+4)',
                'OUTPUT':isochrones})
#%%
#associate_iso_w_core = os.path.join(dir_01_outputs, "join_iso_to_nc.shp")
associate_iso_w_core = processing.run("native:joinattributesbylocation", 
               {'INPUT':isochrones,
                'JOIN':new_cores,
                'PREDICATE':[0,1,2,3,4,5,6],
                'JOIN_FIELDS':['new_id'],
                'METHOD':2,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'OUTPUT':'memory:'})

associate_iso_w_core = associate_iso_w_core['OUTPUT']

#Join this table to non core clusters
#join_rslt_to_nc = os.path.join(dir_01_outputs, "join_rslt_to_nc.shp")
join_rslt_to_nc = processing.run("native:joinattributestable", 
               {'INPUT':non_core_clusters,
                'FIELD':'Id',
                'INPUT_2':associate_iso_w_core,
                'FIELD_2':'nci_id',
                'FIELDS_TO_COPY':['nci_id','new_id'],
                'METHOD':1,
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'OUTPUT':'memory:'})

join_rslt_to_nc = join_rslt_to_nc['OUTPUT']

#Join non core clusters to cores by field
#join_nc_to_core = os.path.join(dir_01_outputs, "join_nc_to_core.shp")
join_nc_to_core = processing.run("native:mergevectorlayers", 
               {'LAYERS':[new_cores,
                          join_rslt_to_nc],
                'CRS':None,
                'OUTPUT':'memory:'})

join_nc_to_core = join_nc_to_core['OUTPUT']

#Dissolve by field
#dissolve_id = os.path.join(dir_01_outputs, "dissolve_test.shp")
dissolve_id = processing.run("native:dissolve", 
               {'INPUT':join_nc_to_core,
                'FIELD':['new_id'],
                'OUTPUT':'memory:'})

dissolve_id = dissolve_id['OUTPUT']

#Remove null agglomerations
#This creates our final agglomerations which will be used as input for 02_script

rm_null = processing.run("native:extractbyexpression",
               {'INPUT':dissolve_id,
                'EXPRESSION':' "new_id" IS NOT NULL',
                'OUTPUT':'memory:'})

rm_null = rm_null['OUTPUT']

rm_fields = processing.run("native:retainfields", 
               {'INPUT':rm_null,
                'FIELDS':['new_id'],
                'OUTPUT':'memory:'})

rm_fields = rm_fields['OUTPUT']

unconfigured_aggloms = os.path.join(dir_01_outputs, "unconfigured_aggloms.shp")

processing.run("native:renametablefield", 
               {'INPUT':rm_fields,
                'FIELD':'new_id',
                'NEW_NAME':'Id',
                'OUTPUT':unconfigured_aggloms})

#%% Remove large intermediary files if wanted

# =============================================================================
# shutil.rmtree(dir_iso_outputs)
# shutil.rmtree(meancoords_dir)
# shutil.rmtree(nc_iso_dir)
# shutil.rmtree(nc_mc_dir)
# shutil.rmtree(nc_pr_dir)
# shutil.rmtree(dir_ncc)
# shutil.rmtree(pop_raster_dir)
# shutil.rmtree(split_vctr_dir)
# =============================================================================
