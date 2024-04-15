# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 10:56:56 2022

@author: orioncr

Script for calculating LCR, PGR, LCRPGR, and supporting metrics

Replace country abbreviation line 68
Replace country boundary line 122

"""
# %% Import libraries and setup QGIS

from processing.core.Processing import Processing
import processing
from qgis.analysis import QgsNativeAlgorithms
import os
import numpy as np
import rasterio
import skimage
import skimage.filters.rank as rank
from skimage.morphology import diamond 
import math
import plotly
import imageio
from matplotlib import pyplot as plt
from pathlib import Path
from osgeo import gdal
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
from qgis.PyQt.QtCore import QVariant
from qgis.utils import iface
from qgis.analysis import QgsRasterCalculatorEntry
# Initialize QGIS
# Supply path to qgis install location
# You can get this using QgsApplication.prefixPath() in QGIS Python console
QgsApplication.setPrefixPath('C:/PROGRA~1/QGIS32~1.4/apps/qgis-ltr', True)
# Create a reference to the QgsApplication.  Setting the
# second argument to False disables the GUI.
qgs = QgsApplication([], False)
# Load providers
qgs.initQgis()

# Add provider to import native algorithms and initialize processing
Processing.initialize()

# QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())

# %% Set General Variables
# Enter the abbreviation of your country or region of interest here (as a string)
country_abbrv = 'XX'

# Insert appropriate projection
# Current example: Africa Albers Equal Area
projection = 'ESRI:102022'

# Threshold value for initial hotspots
LCRPGR_threshold = '1'

# Set agglomeration density threshold to remove low density agglomerations if needed
agglom_density_threshold = ''

# %% Directories
dir_temp = r"Path to folder" + country_abbrv +"_temp"

if not os.path.exists(dir_temp):
    os.makedirs(dir_temp)

dir_00_outputs_t1 = r"Path to 00_outputs_t1"

dir_00_outputs_t2 = r"Path to 00_outputs_t2"

dir_01_outputs_t1 = r"Path to 01_outputs_t1"

dir_01_outputs_t2 = r"Path to 01_outputs_t2"

dir_02_outputs = r"Path to 02_outputs_"

dir_03_outputs = r"Path to 03_outputs"

if not os.path.exists(dir_03_outputs):
    os.makedirs(dir_03_outputs)
    
dir_country_boundaries = r"Path to country boundary data"
  
#for f in os.listdir(dir_03_outputs):
#    os.remove(os.path.join(dir_03_outputs, f))

# %% Set Variables
# Delineated agglomeration vectors
aggloms_t1 = os.path.join(
    dir_02_outputs, country_abbrv + "_t1_Agglomerations.shp")

aggloms_t2 = os.path.join(
    dir_02_outputs, country_abbrv + "_t2_Agglomerations.shp")

# Developed rasters
developed_t1 = os.path.join(
    dir_00_outputs_t1, "Developed_" + country_abbrv + "_t1.tif")

developed_t2 = os.path.join(
    dir_00_outputs_t2, "Developed_" + country_abbrv + "_t2.tif")

# Urbanized open space rasters
uos_t1 = os.path.join(
    dir_00_outputs_t1, "Urbanized_Open_Space_" + country_abbrv + "_t1.tif")

uos_t2 = os.path.join(
    dir_00_outputs_t2, "Urbanized_Open_Space_" + country_abbrv + "_t2.tif")

# Country Boundary
#Insert country boundary shapefile name
ctry_bndry_fp = os.path.join(dir_country_boundaries, 'Country boundary shapefile name')

# Raster characteristics
raster = rasterio.open(developed_t1)
#raster_t2 = rasterio.open(developed_t2)

cell_size = raster.res[0]

cell_area = cell_size**2

# %% Pre-Processing

country_boundary = os.path.join(dir_03_outputs, "Country_Boundary_" + country_abbrv + ".shp")

processing.run("native:reprojectlayer", {'INPUT':ctry_bndry_fp,
                                         'TARGET_CRS':QgsCoordinateReferenceSystem(projection),
                                         'OPERATION':'+proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad +step +proj=aea +lat_0=0 +lon_0=25 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0 +ellps=WGS84',
                                         'OUTPUT':country_boundary})

# Merge agglomerations from t1 and t2 so we can access when needed
t1_and_t2_merged = processing.run("native:joinattributestable",
                                  {'INPUT': aggloms_t2,
                                   'FIELD': 'ID',
                                   'INPUT_2': aggloms_t1,
                                   'FIELD_2': 'ID',
                                   'FIELDS_TO_COPY': [],
                                   'METHOD': 1,
                                   'DISCARD_NONMATCHING': False,
                                   'PREFIX': '',
                                   'OUTPUT': 'memory:'})

t1_and_t2_merged = t1_and_t2_merged['OUTPUT']

# %% SDG 11.3.1 Calculations
#Secondary indicators
#Built up area per capita
BUA_per_capita_t1 = processing.run("native:fieldcalculator", 
                                {'INPUT':t1_and_t2_merged,
                                 'FIELD_NAME':'BUA_per_t1',
                                 'FIELD_TYPE':0,
                                 'FIELD_LENGTH':0,
                                 'FIELD_PRECISION':0,
                                 'FORMULA':' round(("t1_BUA(m²)"  /  "t1_Pop"), 2) ',
                                 'OUTPUT':'memory:'})

BUA_per_capita_output_t1 = BUA_per_capita_t1['OUTPUT']

BUA_per_capita_t2 = processing.run("native:fieldcalculator", 
                                {'INPUT':BUA_per_capita_output_t1,
                                 'FIELD_NAME':'BUA_per_t2',
                                 'FIELD_TYPE':0,
                                 'FIELD_LENGTH':0,
                                 'FIELD_PRECISION':0,
                                 'FORMULA':' round(("t2_BUA(m²)"  /  "t2_Pop"), 2) ',
                                 'OUTPUT':'memory:'})

BUA_per_capita_output_t2 = BUA_per_capita_t2['OUTPUT']

#Total change in built up area
total_change_in_BUA = processing.run("native:fieldcalculator", 
                                     {'INPUT':BUA_per_capita_output_t2,
                                      'FIELD_NAME':'%_chg_BUA',
                                      'FIELD_TYPE':1,
                                      'FIELD_LENGTH':0,
                                      'FIELD_PRECISION':0,
                                      'FORMULA':'(( "t2_BUA(m²)" -  "t1_BUA(m²)" ) /  ("t1_BUA(m²)" )) * 100',
                                      'OUTPUT':'memory:'})

total_change_in_BUA_output = total_change_in_BUA['OUTPUT']

# Compute LCR
LCR = processing.run("native:fieldcalculator",
                     {'INPUT': total_change_in_BUA_output,
                      'FIELD_NAME': 'LCR',
                      'FIELD_TYPE': 0,
                      'FIELD_LENGTH': 0,
                      'FIELD_PRECISION': 0,
                      'FORMULA': ' round((("t2_BUA(m²)" -  "t1_BUA(m²)")/ "t1_BUA(m²)" ) * (1/(5)), 6) ',
                      'OUTPUT': 'memory:'})

LCR_output = LCR['OUTPUT']

# Compute PGR
PGR = processing.run("native:fieldcalculator",
                     {'INPUT': LCR_output,
                      'FIELD_NAME': 'PGR',
                      'FIELD_TYPE': 0,
                      'FIELD_LENGTH': 0,
                      'FIELD_PRECISION': 0,
                      'FORMULA': 'round((ln(  "t2_Pop" /  "t1_Pop" ))/ (5), 6)',
                      'OUTPUT': 'memory:'})

PGR_output = PGR['OUTPUT']

#LCRPGR or LUE
LCRPGR_output = os.path.join(dir_03_outputs, country_abbrv + "_LCRPGR.shp")

LCRPGR = processing.run("native:fieldcalculator",
                        {'INPUT':PGR_output,
                         'FIELD_NAME': 'LCRPGR',
                         'FIELD_TYPE': 0,
                         'FIELD_LENGTH': 0,
                         'FIELD_PRECISION': 0,
                         'FORMULA': ' round(("LCR" /  "PGR"), 6) ',
                         'OUTPUT': LCRPGR_output})
# %% Development Characteristics
# In the Atlas of Urban Expansion, infill is defined as the new development that
# occurs in time 2 that falls within the urbanized open space of time 1

infill_calc = os.path.join(dir_temp, "infill_calc_temp.tif")

processing.run("qgis:rastercalculator",
               # Change expression values to reflect raster and band 
               {'EXPRESSION': '"Developed_NG_t2@1" + "Urbanized_Open_Space_NG_t1@1"',
                'LAYERS': [developed_t2, uos_t1],
                'CELLSIZE': 0,
                'EXTENT': None,
                'CRS': None,
                'OUTPUT': infill_calc})

infill_rstr = os.path.join(dir_03_outputs, "Infill.tif")

processing.run("native:reclassifybytable",
               {'INPUT_RASTER': infill_calc,
                'RASTER_BAND': 1,
                'TABLE': ['0', '0', '0', '1', '1', '0', '2', '2', '1', '3', '1000000', '0'],
                'NO_DATA': 0,
                'RANGE_BOUNDARIES': 2,
                'NODATA_FOR_MISSING': False,
                'DATA_TYPE': 0,
                'OUTPUT': infill_rstr})

infill_zonal = processing.run("native:zonalstatisticsfb",
                              {'INPUT': LCRPGR_output,
                               'INPUT_RASTER': infill_rstr,
                               'RASTER_BAND': 1,
                               'COLUMN_PREFIX': 'if_',
                               'STATISTICS': [1],
                               'OUTPUT': 'memory:'})

infill_zonal = infill_zonal['OUTPUT']


infill_zonal_fp = os.path.join(dir_03_outputs, "infill_zonal.shp")

processing.run("native:fieldcalculator",
               {'INPUT': infill_zonal,
                'FIELD_NAME': 'infill',
                'FIELD_TYPE': 1,
                'FIELD_LENGTH': 0,
                'FIELD_PRECISION': 0,
                'FORMULA': ' "if_sum" *' + str(cell_area),
                'OUTPUT': infill_zonal_fp})
# %% Extension

# Subtract the time 1 developed raster from the time 2 developed raster
ext_subtract = os.path.join(dir_temp, "ext_substract.tif")

processing.run("qgis:rastercalculator",
               # Change expression values to reflect raster and band 
               {'EXPRESSION': '"Developed_NG_t2@1" - "Developed_NG_t1@1"',
                'LAYERS': [developed_t2, developed_t1],
                'CELLSIZE': 0,
                'EXTENT': None,
                'CRS': None,
                'OUTPUT': ext_subtract})

# Reclassify
ext_reclass = os.path.join(dir_temp, "ext_reclass.tif")

processing.run("native:reclassifybytable",
               {'INPUT_RASTER': ext_subtract,
                'RASTER_BAND': 1,
                'TABLE': ['-1', '-1', '0','0', '0', '0', '1', '1', '1'],
                'NO_DATA': -9999,
                'RANGE_BOUNDARIES': 2,
                'NODATA_FOR_MISSING': False,
                'DATA_TYPE': 5 ,
                'OUTPUT': ext_reclass})

ext_translate = os.path.join(dir_temp, "ext_translate.tif")

processing.run("gdal:translate",            
               {'INPUT':ext_reclass,
                'TARGET_CRS':None,
                'NODATA':None,
                'COPY_SUBDATASETS':False,
                'OPTIONS':'',
                'EXTRA':'',
                'DATA_TYPE':1,
                'OUTPUT':ext_translate})

# Find difference between country boundaries and agglomerations
diff = os.path.join(dir_temp, "diff.shp")

processing.run("native:difference",
               {'INPUT':country_boundary,
                'OVERLAY':aggloms_t1,
                'OUTPUT':diff})

# Clip this reclassifed layer to agglomeration boundaries of time 1
clip_reclass = os.path.join(dir_temp, "clip_reclass.tif")

processing.run("gdal:cliprasterbymasklayer",
               {'INPUT': ext_translate,
                'MASK': diff,
                'SOURCE_CRS': None,
                'TARGET_CRS': None,
                'NODATA': None,
                'ALPHA_BAND': False,
                'CROP_TO_CUTLINE': True,
                'KEEP_RESOLUTION': False,
                'SET_RESOLUTION': False,
                'X_RESOLUTION': None,
                'Y_RESOLUTION': None,
                'MULTITHREADING': False,
                'OPTIONS': '',
                'DATA_TYPE': 0,
                'EXTRA': '',
                'OUTPUT': clip_reclass})


#%%
# Find contiguous developed pixel clusters
# Start with image filtering
# Use rank sum from skimage to find pixels with neighbors of a value greater than 0

with rasterio.open(clip_reclass, 'r') as ds:
    arr = ds.read(1)  # read all raster values

print(arr)
print(arr.shape)

rook = diamond(1, dtype=np.uint8)

rank_sum = rank.sum(arr, rook)

cont_dev = os.path.join(dir_temp, 'contiguous_developed.tif')
x_pixels = len(rank_sum[1])
y_pixels = len(rank_sum)
driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(cont_dev,x_pixels, y_pixels, 1,gdal.GDT_Byte)
dataset.GetRasterBand(1).WriteArray(rank_sum)

file = clip_reclass
data0 = gdal.Open(file)
geotrans=data0.GetGeoTransform()  #get GeoTranform from existed 'data0'
proj=data0.GetProjection() #you can get from an existing tif
dataset.SetGeoTransform(geotrans)
dataset.SetProjection(proj)
dataset.FlushCache()
dataset=None
#%%
# Reclassify and remove pixels that have less than a certain amount of neighbors
cont_pixels = os.path.join(dir_temp, "cont_pixels.tif")

processing.run("native:reclassifybytable",
               {'INPUT_RASTER':cont_dev,
                'RASTER_BAND':1,
                'TABLE':['0','1','0','2','256','1'],
                'NO_DATA':-9999,
                'RANGE_BOUNDARIES':2,
                'NODATA_FOR_MISSING':False,
                'DATA_TYPE':0,
                'OUTPUT':cont_pixels})

# Raster calculator to determine where clip_reclass and cont_pixels overlap

ext_potential = os.path.join(dir_temp, "ext_potential.tif")

processing.run("qgis:rastercalculator",
               # Change expression values to reflect raster and band 
               {'EXPRESSION': '"clip_reclass@1" + "cont_pixels@1"',
                'LAYERS': [clip_reclass, cont_pixels],
                'CELLSIZE': 0,
                'EXTENT': None,
                'CRS': None,
                'OUTPUT': ext_potential})

pot_translate = os.path.join(dir_temp, "pot_translate.tif")

processing.run("gdal:translate",            
               {'INPUT':ext_potential,
                'TARGET_CRS':None,
                'NODATA':None,
                'COPY_SUBDATASETS':False,
                'OPTIONS':'',
                'EXTRA':'',
                'DATA_TYPE':1,
                'OUTPUT':pot_translate})

### Not working because of missing gdal.bat file
# =============================================================================
# vectorize_rast = processing.run("gdal:polygonize", 
#                 {'INPUT':pot_translate,
#                 'BAND':1,
#                 'FIELD':'COUNT',
#                 'EIGHT_CONNECTEDNESS':False,
#                 'EXTRA':'',
#                 'OUTPUT':'memory:'})
# =============================================================================

vectorize_rast = processing.run("native:pixelstopolygons", 
               {'INPUT_RASTER':pot_translate,
                'RASTER_BAND':1,
                'FIELD_NAME':'COUNT',
                'OUTPUT':'memory:'})

vectorize_rast = vectorize_rast['OUTPUT']

# Had to manually create the vectorized raster because gdal isnt working with python
#vectorize_rast = os.path.join(dir_temp, "vectorize_rast.shp")

# Extract by expression
count_2_vctr = processing.run("native:extractbyexpression", 
               {'INPUT':vectorize_rast,
                'EXPRESSION':'"COUNT" = 2',
                'OUTPUT':'memory:'})

count_2_vctr = count_2_vctr['OUTPUT']
# Find which of these polygons overlap with the boundaries of time 1 agglomerations
# and extract
ext_polygons = processing.run("native:extractbylocation", 
               {'INPUT':count_2_vctr,
                'PREDICATE':[0,1,3,4,5,6,7],
                'INTERSECT':aggloms_t1,
                'OUTPUT':'memory:'})

ext_polygons = ext_polygons['OUTPUT']
# Extract by mask: Find the pixels that met the neighborhood criteria that fall within the overlapping polygons

# These are the extension pixels
extension = os.path.join(dir_temp, "extension.tif")

processing.run("gdal:cliprasterbymasklayer", 
               {'INPUT':pot_translate,
                'MASK':ext_polygons,
                'SOURCE_CRS':None,
                'TARGET_CRS':None,
                'NODATA':None,
                'ALPHA_BAND':False,
                'CROP_TO_CUTLINE':True,
                'KEEP_RESOLUTION':False,
                'SET_RESOLUTION':False,
                'X_RESOLUTION':None,
                'Y_RESOLUTION':None,
                'MULTITHREADING':False,
                'OPTIONS':'',
                'DATA_TYPE':0,
                'EXTRA':'',
                'OUTPUT':extension})
#Calculate zonal statistics for these pixels within the T2 agglomeration

ext_zonal = processing.run("native:zonalstatisticsfb", 
               {'INPUT':infill_zonal_fp,
                'INPUT_RASTER':extension,
                'RASTER_BAND':1,
                'COLUMN_PREFIX':'ext_',
                'STATISTICS':[0,1,2],
                'OUTPUT':'memory:'})

ext_zonal = ext_zonal['OUTPUT']

#Calculate extension area for each agglomeration
ext_area = processing.run("native:fieldcalculator",
               {'INPUT': ext_zonal,
                'FIELD_NAME': 'extension',
                'FIELD_TYPE': 1,
                'FIELD_LENGTH': 0,
                'FIELD_PRECISION': 0,
                'FORMULA': ' "ext_count" *' + str(cell_area),
                'OUTPUT': 'memory:'})

ext_area = ext_area['OUTPUT']

field_names_aggloms = ext_area.fields().names()
print(field_names_aggloms)

#%%
# Clean shapefile to obtain only fields wanted
agglom_w_calcs = os.path.join(dir_03_outputs, "aggloms_w_calcs.shp")

retain_fields = processing.run("native:retainfields", 
               {'INPUT': ext_area,
               'FIELDS':['ID', 'Area(m²)', 'Area(m²)_2', 't1_Pop', 't2_Pop', 't1_BUA(m²)', 't2_BUA(m²)', 'BUA_per_t1', 'BUA_per_t2',
                         '%_chg_BUA', 'LCR', 'PGR', 'LCRPGR', 'infill', 'extension'],
               'OUTPUT':'memory:'})

retain_fields_output = retain_fields['OUTPUT']

processing.run("native:refactorfields", 
               {'INPUT':retain_fields_output,
                'FIELDS_MAPPING':[{'expression': '"ID"','length': 10,'name': 'ID','precision': 0,'type': 4},
                                  {'expression': '"Area(m²)_2"','length': 20,'name': 'Area(m²)_1','precision': 0,'type': 2},
                                  {'expression': '"Area(m²)"','length': 20,'name': 'Area(m²)_2','precision': 0,'type': 2},
                                  {'expression': '"t1_Pop"','length': 20,'name': 't1_Pop','precision': 0,'type': 2},
                                  {'expression': '"t2_Pop"','length': 20,'name': 't2_Pop','precision': 0,'type': 2},
                                  {'expression': '"t1_BUA(m²)"','length': 20,'name': 't1_BUA(m²)','precision': 0,'type': 2},
                                  {'expression': '"t2_BUA(m²)"','length': 20,'name': 't2_BUA(m²)','precision': 0,'type': 2},
                                  {'expression': '"BUA_per_t1"','length': 23,'name': 'BUA_per_t1','precision': 15,'type': 2},
                                  {'expression': '"BUA_per_t2"','length': 23,'name': 'BUA_per_t2','precision': 15,'type': 2},
                                  {'expression': '"%_chg_BUA"','length': 10,'name': '%_chg_BUA','precision': 0,'type': 4},
                                  {'expression': '"LCR"','length': 23,'name': 'LCR','precision': 6,'type': 6},
                                  {'expression': '"PGR"','length': 23,'name': 'PGR','precision': 6,'type': 6},
                                  {'expression': '"LCRPGR"','length': 23,'name': 'LCRPGR','precision': 6,'type': 6},
                                  {'expression': '"infill"','length': 10,'name': 'infill','precision': 0,'type': 4},
                                  {'expression': '"extension"','length': 10,'name': 'extension','precision': 0,'type': 4}],
                'OUTPUT':agglom_w_calcs})
# =============================================================================
# for f in os.listdir(dir_temp):
#      os.remove(os.path.join(dir_temp, f))
# =============================================================================
     
#%% Associate Places With Agglomerations
city_data = os.path.join(dir_01_outputs_t1,country_abbrv + "_OSM_cities.shp")
 
town_data = os.path.join(dir_01_outputs_t1,country_abbrv + "_OSM_towns.shp")

# Cities will be joined first

city_within_core = processing.run("native:joinbynearest", 
               {'INPUT':agglom_w_calcs,
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

# Now we join towns

town_within_core = processing.run("native:joinbynearest", 
               {'INPUT':ncc_output,
                'INPUT_2':town_data,
                'FIELDS_TO_COPY':[],
                'DISCARD_NONMATCHING':False,
                'PREFIX':'',
                'NEIGHBORS':1,
                'MAX_DISTANCE':5000,
                'OUTPUT':'memory:'})

town_within_core_output = town_within_core['OUTPUT']

#town_cores = os.path.join(dir_01_outputs, 't1_Town_Cores.shp')

town_cores = processing.run("native:extractbyexpression", 
                            {'INPUT':town_within_core_output,
                             'EXPRESSION':' "full_id" IS NOT NULL',
                             'OUTPUT':'memory:'})

town_cores = town_cores['OUTPUT']

ntc = processing.run("native:extractbyexpression",
                     {'INPUT':town_within_core_output,
                      'EXPRESSION':' "full_id" IS NULL',
                      'OUTPUT':'memory:'})

ntc_output = ntc['OUTPUT']

combine_cores = processing.run("native:mergevectorlayers", 
                       {'LAYERS':[core_with_locality, ntc_output],
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

add_field = processing.run("native:addautoincrementalfield", {'INPUT':rm_duplicates,
                                                  'FIELD_NAME':'Unique_ID',
                                                  'START':0,
                                                  'MODULUS':0,
                                                  'GROUP_FIELDS':[],                 
                                                  'SORT_EXPRESSION':'',
                                                  'SORT_ASCENDING':True,
                                                  'SORT_NULLS_FIRST':False,
                                                  'OUTPUT':'memory:'})

add_field_op = add_field['OUTPUT']

# we will just use year two for density calculations
# remove aggloms with a density less than set threshold to remove outliers and filter
# misclassified areas

calc_density = processing.run("native:fieldcalculator", 
                              {'INPUT':add_field_op,
                               'FIELD_NAME':'density',
                               'FIELD_TYPE':1,
                               'FIELD_LENGTH':0,
                               'FIELD_PRECISION':0,
                               'FORMULA':' "t2_Pop"  /  ("Area(m²)_2" / 1000000) ',
                               'OUTPUT':'memory:'})

cd_output = calc_density ['OUTPUT']

rm_outliers = processing.run("native:extractbyexpression", {'INPUT':cd_output,
                                                            'EXPRESSION':' "density" >= ' + agglom_density_threshold,
                                                            'OUTPUT':'memory:'})

rm_output = rm_outliers ['OUTPUT']

final_aggloms = os.path.join(dir_03_outputs, country_abbrv + "_final_agglomerations.shp")

drop_fields = processing.run("native:deletecolumn", 
               {'INPUT':rm_output,
                'COLUMN':['full_id','osm_id','distance','feature_x','feature_y',
                          'nearest_x','nearest_y', 'n','full_id_2','osm_id_2','full_id_3',
                          'osm_id_3','full_id_4','osm_id_4','layer','path','Unique_ID'],
                'OUTPUT':final_aggloms})
#%% Extract agglomerations with LCRPGR greater than threshold value
LCRPGR_hotspot = os.path.join(dir_03_outputs, country_abbrv + "_LCRPGR_hotspots.shp")

processing.run("native:extractbyexpression", 
               {'INPUT':final_aggloms,
                'EXPRESSION':' "LCRPGR" >= ' + LCRPGR_threshold,
                'OUTPUT':LCRPGR_hotspot})

LCRPGR_hotspot_csv = os.path.join(dir_03_outputs, country_abbrv + "_LCRPGR_hotspots.csv")

processing.run("native:exporttospreadsheet", 
               {'LAYERS':LCRPGR_hotspot,
                'USE_ALIAS':False,
                'FORMATTED_VALUES':False,
                'OUTPUT':LCRPGR_hotspot_csv,
                'OVERWRITE':True})

#%%Categorize agglomerations by population sizes if desired

hotspots_pop_50k_or_less = os.path.join(dir_03_outputs, country_abbrv + "_hotspots_pop_50k_or_less.shp")

processing.run("native:extractbyexpression", 
               {'INPUT':LCRPGR_hotspot,
                'EXPRESSION':' "t2_Pop" <= 49999',
                'OUTPUT':hotspots_pop_50k_or_less})

hotspots_pop_50k_to_100k = os.path.join(dir_03_outputs, country_abbrv + "_hotspots_pop_50k_to_100k.shp")

processing.run("native:extractbyexpression", 
               {'INPUT':LCRPGR_hotspot,
                'EXPRESSION':' "t2_Pop" >= 50000 AND  "t2_Pop" <= 99999',
                'OUTPUT':hotspots_pop_50k_to_100k})

hotspots_pop_greater_than_100k = os.path.join(dir_03_outputs, country_abbrv + "_hotspots_pop_greater_than_100k.shp")

processing.run("native:extractbyexpression", 
               {'INPUT':LCRPGR_hotspot,
                'EXPRESSION':' "t2_Pop" >= 100000',
                'OUTPUT':hotspots_pop_greater_than_100k})


