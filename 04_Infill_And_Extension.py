# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 10:56:56 2022

@author: orioncr

Script for calculating LCR, PGR, LCRPGR or LUE, Supporting Metrics,
Growth Patterns and Green Space

This calculates the correct infill and extension values

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

#QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())

# %% Set General Variables
# Enter the abbreviation of your country of interest here (as a string)
country_abbrv = 'ET_90m'

projection = 'ESRI:102022'

LCRPGR_threshold = '1'

agglom_density_threshold = '1500'

# %% Directories
dir_temp = r"D:\/" + country_abbrv +"_temp"

#for f in os.listdir(dir_temp):
#    os.remove(os.path.join(dir_temp, f))
dir_00_outputs_t1 = r"N:/RStor/jvogeler/Lab/users/orioncr/scripts/outputs/00_derive_urban_extent/" + \
    country_abbrv + "/t1"

dir_00_outputs_t2 = r"N:/RStor/jvogeler/Lab/users/orioncr/scripts/outputs/00_derive_urban_extent/" + \
    country_abbrv + "/t2"

dir_01_outputs_t1 = r"N:/RStor/jvogeler/Lab/users/orioncr/scripts/outputs/01_configure_agglomerations/" + \
    country_abbrv + "/t1/outputs"

dir_01_outputs_t2 = r"N:/RStor/jvogeler/Lab/users/orioncr/scripts/outputs/01_configure_agglomerations/" + \
    country_abbrv + "/t2/outputs"

dir_02_outputs = r"N:/RStor/jvogeler/Lab/users/orioncr/scripts/outputs/02_final_agglomerations/" + country_abbrv

dir_03_outputs = r"N:/RStor/jvogeler/Lab/users/orioncr/scripts/outputs/03_calculations/" + country_abbrv

#if not os.path.exists(dir_03_outputs):
#    os.makedirs(dir_03_outputs)
    
dir_country_boundaries = r"N:/RStor/jvogeler/Lab/users/orioncr/country_data/countryboundaries/" + country_abbrv
  
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

final_aggloms = os.path.join(dir_03_outputs, country_abbrv + "_final_agglomerations.shp")

# Country Boundary
ctry_bndry_fp = os.path.join(dir_country_boundaries, 'eth_admbnda_adm0_csa_bofedb_itos_2021.shp')

# Raster characteristics
raster = rasterio.open(developed_t1)
#raster_t2 = rasterio.open(developed_t2)

cell_size = raster.res[0]

cell_area = cell_size**2

# %% Hotspots
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
               
#%%Categorize agglomerations by population sizes

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


# %% Infill

dv_change = os.path.join(dir_temp, "change_t1_to_t2.tif")

processing.run("qgis:rastercalculator",
               {'EXPRESSION':'"Developed_ET_t2@1" - "Developed_ET_t1@1"',
                'LAYERS': [developed_t2, developed_t1],
                'CELLSIZE': 0,
                'EXTENT': None,
                'CRS': None,
                'OUTPUT': dv_change})

new_development = os.path.join(dir_03_outputs, "new_development.tif")

processing.run("native:reclassifybytable",
               {'INPUT_RASTER': dv_change,
                'RASTER_BAND': 1,
                'TABLE': ['-1', '-1', '0', '1', '1', '1'],
                'NO_DATA': 0,
                'RANGE_BOUNDARIES': 2,
                'NODATA_FOR_MISSING': False,
                'DATA_TYPE': 1,
                'OUTPUT': new_development})

clip_infill = os.path.join(dir_03_outputs, "infill_clip.tif")

processing.run("gdal:cliprasterbymasklayer",
               {'INPUT':new_development,
                'MASK':aggloms_t1,
                'SOURCE_CRS':None,
                'TARGET_CRS':None,
                'NODATA':0,
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
                'OUTPUT':clip_infill})

infill_raster = os.path.join(dir_03_outputs, "infill_raster.tif")

processing.run("native:fillnodata", 
               {'INPUT':clip_infill,
                'BAND':1,
                'FILL_VALUE':0,
                'OUTPUT':infill_raster})

infill_zonal = processing.run("native:zonalstatisticsfb",
                              {'INPUT': LCRPGR_hotspot,
                               'INPUT_RASTER': infill_raster,
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
                'FORMULA': ' "if_sum" * ' + str(cell_area),
                'OUTPUT': infill_zonal_fp})

infill_display_rstr = os.path.join(dir_03_outputs, "infill_display_rstr.tif")

processing.run("gdal:cliprasterbymasklayer",
               {'INPUT':infill_raster,
                'MASK':LCRPGR_hotspot,
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
                'OUTPUT':infill_display_rstr})

# %% Extension and leapfrog
new_development_fill = os.path.join(dir_03_outputs, "new_development_fill.tif")

processing.run("native:fillnodata", 
               {'INPUT':new_development,
                'BAND':1,
                'FILL_VALUE':0,
                'OUTPUT':new_development_fill})

ext_lf_new = os.path.join(dir_temp, "ext_lf_new.tif")

processing.run("qgis:rastercalculator", 
               {'EXPRESSION':'"new_development_fill@1" - "infill_raster@1"',
                'LAYERS':[new_development_fill, infill_raster],
                'CELLSIZE':0,
                'EXTENT':None,
                'CRS':None,
                'OUTPUT':ext_lf_new})

ext_lf_zonal = processing.run("native:zonalstatisticsfb",
                              {'INPUT': infill_zonal_fp,
                               'INPUT_RASTER': ext_lf_new,
                               'RASTER_BAND': 1,
                               'COLUMN_PREFIX': 'ext_',
                               'STATISTICS': [1],
                               'OUTPUT': 'memory:'})

ext_lf_zonal = ext_lf_zonal['OUTPUT']

final_with_dev_calcs = os.path.join(dir_03_outputs, country_abbrv + "_final_with_dev_calcs.shp")

processing.run("native:fieldcalculator",
               {'INPUT': ext_lf_zonal,
                'FIELD_NAME': 'extension',
                'FIELD_TYPE': 1,
                'FIELD_LENGTH': 0,
                'FIELD_PRECISION': 0,
                'FORMULA': ' "ext_sum" * ' + str(cell_area),
                'OUTPUT': final_with_dev_calcs})

ext_lf_display_rstr = os.path.join(dir_03_outputs, "ext_lf_display_rstr.tif")

processing.run("gdal:cliprasterbymasklayer",
               {'INPUT':ext_lf_new,
                'MASK':ext_lf_zonal,
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
                'OUTPUT':ext_lf_display_rstr})

new_development_clip = os.path.join(dir_03_outputs, "new_development_clip.tif")

processing.run("gdal:cliprasterbymasklayer",
               {'INPUT':new_development_fill,
                'MASK':final_with_dev_calcs,
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
                'OUTPUT':new_development_clip})

new_da = os.path.join(dir_03_outputs, "new_da.html")

processing.run("native:rasterlayeruniquevaluesreport",
               {'INPUT':new_development_clip,
                'BAND':1,
                'OUTPUT_HTML_FILE':new_da})







