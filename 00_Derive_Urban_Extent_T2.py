# -*- coding: utf-8 -*-
"""
Created on June 17th 2022

@author: orioncr
00_Derive_Urban_Extent_T2

NOTES
Use the Urban Extent method to extract Urban, Suburban, Captured
Open Space, and Fringe Open Space layers from built-up land cover or developed land use data.
This script spatially analyzes data using arcpy.
"""

#%% Import libraries
import os
import arcpy
import math
import rasterio

# All outputs will be in the following coordinate system
# Substitute appropriate coordinate system. Africa Albers Equal Area Conic is 
# used as an example
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("Africa Albers Equal Area Conic")

# Overwrite outputs
arcpy.env.overwriteOutput = True

#%% Set General Variables
#Enter the abbreviation of your country or region of interest here (as a string)
country_abbrv = "XX"

#Enter the file name of the land cover or land use data
lclu_tile = "XX.tif"

#Assign projection spatial reference to object
projection = arcpy.env.outputCoordinateSystem

#Define variables
#The cpatured open space variable determines the maximum size of a captured open space.
#In the Atlas of Urban Expansion, captured open space only considers areas less than 200 hectares in area
#In meters, this value is 2,000,000
captured_open_space_area = '2000000'

#The developed lclu value variable is the value that represents the developed 
#land use class or built-up land cover class in the input map. 
developed_lclu_value = 0

#%% Set Directories
#Set directory to store outputs
dir_outputs = r""

if not os.path.exists(dir_outputs):
    os.makedirs(dir_outputs)

#Set directory for accessing land use or land cover maps
dir_landcoverdata = r""

#Set directory for accessing gridded population data
dir_popdata = r""

#Set directory for accessing country or region boundary data
dir_boundarydata = r"‪"

#Utilize this code chunk if wanting to remove all outputs
# =============================================================================
# for f in os.listdir(dir_outputs):
#     os.remove(os.path.join(dir_outputs, f))
# =============================================================================
    
#%% Pre-processing
### Create raster with built up or developed pixels only

### Load land cover or land use tile/s
year_2 = os.path.join(dir_landcoverdata, lclu_tile)

# These values will need to be changed as input changes
# We use rasterio to obtain information about our input raster such as cell
# size and cell area
raster = rasterio.open(year_2)

cell_size = raster.res[0]

cell_area = cell_size**2

# =============================================================================
# Mosaic to new raster to merge raster tiles
# This is only necessary if raster tile does not cover extent of study region
# arcpy.management.MosaicToNewRaster(input_rasters, 
# output_location, 
# raster_dataset_name_with_extension, 
# {coordinate_system_for_the_raster}, 
# {pixel_type}, 
# {cellsize}, 
# number_of_bands, 
# {mosaic_method},
#  {mosaic_colormap_mode})
# =============================================================================

### Reclass for developed pixels only
### Reclass value for built up class to 1, remaining classes are reclassified as 0
Developed_List = [[0,2,0], [developed_lclu_value,developed_lclu_value,1], [4,7,0]]

Developed_Remap = arcpy.sa.RemapValue(Developed_List)

#fp_Developed_Reclass = os.path.join(dir_outputs, "Developed_Reclass" + country_abbrv + "_t2.tif")

Developed_Reclass = arcpy.sa.Reclassify(year_2,
                                 "Value", 
                                 Developed_Remap,
                                 "DATA"); 
#Developed_Reclass.save(fp_Developed_Reclass)

### Project raster
fp_Developed = os.path.join(dir_outputs, "Developed_" + country_abbrv + "_t2.tif")

Developed = arcpy.management.ProjectRaster(Developed_Reclass,
                                              fp_Developed,
                                              projection)
#%%Processing

### Compute focal statistics to determine density of built-up area or developed
### land use area per square kilometer
#fp_Developed_FocalStat = os.path.join(dir_outputs, "Developed_FocalStat" + country_abbrv + "_t2.tif")

Developed_FocalStat = arcpy.sa.FocalStatistics(fp_Developed,
                                        "Circle 564 MAP",
                                        "SUM",
                                        "DATA",
                                        90); 
#Developed_FocalStat.save(fp_Developed_FocalStat)

### Reclassify the focal statistics output to determine urbaneness

## In the below tool, we have to know how many pixels fall within the walking distance circle
# We use the maximum pixels within the 1 km2 buffer 
pixels_within_1km_circle = Developed_FocalStat.maximum

## Reclassification

fp_USR_Reclass = os.path.join(dir_outputs, "USR_Reclass_" + country_abbrv + "_t2.tif")
    
urban_bua = pixels_within_1km_circle * 0.5

suburban_bua = pixels_within_1km_circle * 0.25
                          
USR_List = [[1,suburban_bua,1], [suburban_bua,urban_bua,2], [urban_bua,pixels_within_1km_circle,3]]

USR_Remap = arcpy.sa.RemapValue(USR_List)

USR_Reclass = arcpy.sa.Reclassify(Developed_FocalStat,
                                  "Value",
                                  USR_Remap,
                                  "DATA");
USR_Reclass.save(fp_USR_Reclass)

# Use Raster Calculator to combine BU layer with USR layer

#fp_BU_USR_RasterCalc = os.path.join(dir_outputs, "BU_USR_" + country_abbrv + "_t2.tif")

BU_USR_RasterCalc =  arcpy.sa.RasterCalculator([USR_Reclass,fp_Developed],
                                             ["x","y"],
                                             "x+y");
#BU_USR_RasterCalc.save(fp_BU_USR_RasterCalc)

### Reclassify output from raster calculation to remove rural pixels
### Step 10
#fp_BU_US_Reclass = os.path.join(dir_outputs, "BU_US_" + country_abbrv + "_t2.tif")

BU_US_List = [[0,"NODATA"],
              [1, "NODATA"],
              [2, "NODATA"],
              [3,2],
              [4,1],
              ["NODATA","NODATA"]]

BU_US_Remap = arcpy.sa.RemapValue(BU_US_List)

BU_US_Reclass = arcpy.sa.Reclassify(BU_USR_RasterCalc,
                                    "Value",
                                    BU_US_Remap,
                                    "DATA");
#BU_US_Reclass.save(fp_BU_US_Reclass)

### Buffer resulting urban and suburban pixels to 100m to define fringe open spaces
### Step 11
#fp_FringeOpenSpace_UC = os.path.join(dir_outputs, "FOS_UC_" + country_abbrv + "_t2.tif")

max_distance = (math.sqrt(((cell_size/2)**2)*(2))) + (100)

FringeOpenSpace_UC = arcpy.sa.EucDistance(BU_US_Reclass,
                                          max_distance,
                                          cell_size,
                                          None,
                                          "PLANAR",
                                          None,
                                          None);
#FringeOpenSpace_UC.save(fp_FringeOpenSpace_UC)

### Reclassify buffer results to a single value
### Step 12
#fp_FringeOpenSpace_Reclass = os.path.join(dir_outputs, "FOS_Reclass_" + country_abbrv + "_t2.tif")

Fringe_Open_Space_List = [[0,10000000,1]]

Fringe_Open_Space_Remap = arcpy.sa.RemapValue(Fringe_Open_Space_List)

FringeOpenSpace_Reclass = arcpy.sa.Reclassify(FringeOpenSpace_UC,
                                              "VALUE",
                                              Fringe_Open_Space_Remap,
                                              "NODATA")

#FringeOpenSpace_Reclass.save(fp_FringeOpenSpace_Reclass)

### Merge urban and suburban pixels with reclassified buffer from previous step
### Step 13

### Define null values for two raster outputs
#fp_DefineNullforBU = os.path.join(dir_outputs, "BU_US_Null_" + country_abbrv + "_t2.tif")

DefineNullforBU = arcpy.sa.RasterCalculator([BU_US_Reclass,BU_US_Reclass],
                                            ["x","y"],
                                            'Con(IsNull(x),0,y)');
#DefineNullforBU.save(fp_DefineNullforBU)

#fp_DefineNullforFOS = os.path.join(dir_outputs, "FOS_Null_" + country_abbrv + "_t2.tif")

DefineNullforFOS = arcpy.sa.RasterCalculator([FringeOpenSpace_Reclass,FringeOpenSpace_Reclass],
                                             ["x","y"],
                                             'Con(IsNull(x),0,y)');
#DefineNullforFOS.save(fp_DefineNullforFOS)

### Merge the two outputs from previous step through Raster Calculator addition
### Step 14
#fp_MergeNullOutputs = os.path.join(dir_outputs, "Merge_Null_BU_FOS_" + country_abbrv + "_t2.tif")

MergeNullOutputs = arcpy.sa.RasterCalculator([DefineNullforBU,DefineNullforFOS],
                                             ["x","y"],
                                             "x+y");
#MergeNullOutputs.save(fp_MergeNullOutputs)

### Set Null value to zero so you are left with urban, suburban and fringe open space pixels
### Step 15
#Null_USFOS = os.path.join(dir_outputs, "Null_USFOS_" + country_abbrv + "_t2.tif")

Null_USFOS = arcpy.sa.RasterCalculator([MergeNullOutputs,MergeNullOutputs],
                                             ["x","y"],
                                             'SetNull(x ==0,y)');
#Null_USFOS.save(Null_USFOS)

### Merge urban, suburban and FOS classes into one using reclassify
### 17
#fp_USFOS_Reclass = os.path.join(dir_outputs, "USFOS_Reclass_" + country_abbrv + "_t2.tif")

USFOSList = [[0, "NODATA"],[1, 1], [2,1], [3,1]]

USFOSRemap = arcpy.sa.RemapValue(USFOSList)

USFOS_Reclass = arcpy.sa.Reclassify(Null_USFOS,
                                    "Value",
                                    USFOSRemap,
                                    "DATA");
#USFOS_Reclass.save(fp_USFOS_Reclass)

### Convert output raster to polygon

fp_USFOS_Polygon = os.path.join(dir_outputs, "USFOS_Polygon_" + country_abbrv + "_t2.shp")

with arcpy.EnvManager(outputZFlag="Disabled", outputMFlag="Disabled"): 
    USFOS_Polygon = arcpy.conversion.RasterToPolygon(USFOS_Reclass, 
                                                     fp_USFOS_Polygon,
                                     "NO_SIMPLIFY",
                                     "Value", 
                                     "SINGLE_OUTER_PART",
                                     None)

#%% Extract Captured Open Space
### Create polygons for open space gaps within USFOS Polygon
#fp_ID_COS_Polygons = os.path.join(dir_outputs, "Identify_COS_Polygons_" + country_abbrv + "_t2.shp")

ID_COS_Polygons = arcpy.analysis.Union(fp_USFOS_Polygon,
                     r"memory/ID_COS_Polygons",
                     "ALL",
                     None,
                     "NO_GAPS")

### Get a list of field names to be used in selection
field_names = [f.name for f in arcpy.ListFields(ID_COS_Polygons)]

print(field_names)

### Obtain all COS Polygons by selecting features with FID of -1
#fp_COS_Polygons = os.path.join(dir_outputs, "COS_Polygons_" + country_abbrv + "_t2.shp")

COS_Polygons = arcpy.analysis.Select(ID_COS_Polygons, 
                      r"memory/COS_Polygons",
                      field_names[2] + " = -1")

### Calculate the area and a value for each polygon
Calc_Area_For_COS_Polygons = arcpy.management.CalculateField(COS_Polygons, 
                                "Area", 
                                "!shape.area!", 
                                "PYTHON3", 
                                '',
                                "LONG",
                                "NO_ENFORCE_DOMAINS")

Calc_Value_For_COS_Polygons = arcpy.management.CalculateField(COS_Polygons,
                                                             "Value",
                                                             "1",
                                                             "PYTHON3",
                                                             '',
                                                             "SHORT",
                                                             "NO_ENFORCE_DOMAINS")

### Select polygons with an area less than threshold
### 200 hectares or 2000000 square meters is the default threshold from Angel et al.
COS_Polygons_Less_Than_Threshold = os.path.join(dir_outputs, "COS_Polygons_Less_Than_Threshold_" + country_abbrv + "_t2.shp")

Select_COS_Polygons_Less_Than = arcpy.analysis.Select(COS_Polygons,
                      COS_Polygons_Less_Than_Threshold,
                      "Area < " + captured_open_space_area)

### Fill holes that are less than 2000000 sq meters in area - this is the captured open space
#fp_fill_COS = os.path.join(dir_outputs, "Fill_COS_" + country_abbrv + "_t2.shp")

fill_COS = arcpy.management.EliminatePolygonPart(fp_USFOS_Polygon, 
                                      r"memory/fill_COS",
                                      "AREA",
                                      captured_open_space_area + " SquareMeters",
                                      0,
                                      "CONTAINED_ONLY")

field_names = [f.name for f in arcpy.ListFields(fill_COS)]

print(field_names)

### Convert polygon to raster
#fp_convert_p_to_r = os.path.join(dir_outputs, "Polygon_To_Raster_" + country_abbrv + "_t2.tif")

convert_p_to_r = arcpy.conversion.PolygonToRaster(fill_COS,
                                 "ORIG_FID",
                                 r"memory\convert_p_to_r",
                                 "CELL_CENTER",
                                 "NONE",
                                 cell_size,
                                 "BUILD")

### Convert raster to polygon
### At this step, it is very important to set max vertices as this will impact
### our use of the OpenRouteService API in the next script. The API currently 
### only allows 3500 routes maximum
fp_urban_clusters = os.path.join(dir_outputs, "Urban_Clusters_" + country_abbrv + "_t2.shp")

arcpy.conversion.RasterToPolygon(convert_p_to_r,
                                 fp_urban_clusters,
                                 "SIMPLIFY",
                                 "Value",
                                 "SINGLE_OUTER_PART",
                                 None)
### Calculate area of each of the resultant polygons using geometry tool

with arcpy.EnvManager(outputCoordinateSystem=None):
    CalcGeo_USFOS_Polygon = arcpy.management.CalculateGeometryAttributes(fp_urban_clusters,
                                                 "Area AREA_GEODESIC",
                                                 '',
                                                 "SQUARE_KILOMETERS",
                                                 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]',
                                                 "SAME_AS_INPUT")
    
#%% Extract individual layers for use in LUE and pattern calculations

#Urban
fp_Urban_Reclass = os.path.join(dir_outputs, "Urban_" + country_abbrv + "_t2.tif")

Urban_List = [[0,"NODATA"],
              [1, "NODATA"],
              [2, "NODATA"],
              [3,"NODATA"],
              [4,1],
              ["NODATA","NODATA"]]

Urban_Remap = arcpy.sa.RemapValue(Urban_List)

Urban_Reclass = arcpy.sa.Reclassify(BU_USR_RasterCalc,
                                    "Value",
                                    Urban_Remap,
                                    "DATA");
Urban_Reclass.save(fp_Urban_Reclass)

#Suburban
fp_Suburban_Reclass = os.path.join(dir_outputs, "Suburban_" + country_abbrv + "_t2.tif")

Suburban_List = [[0,"NODATA"],
              [1, "NODATA"],
              [2, "NODATA"],
              [3,1],
              [4,"NODATA"],
              ["NODATA","NODATA"]]

Suburban_Remap = arcpy.sa.RemapValue(Suburban_List)

Suburban_Reclass = arcpy.sa.Reclassify(BU_USR_RasterCalc,
                                    "Value",
                                    Suburban_Remap,
                                    "DATA");
Suburban_Reclass.save(fp_Suburban_Reclass)

#Rural
fp_Rural_Reclass = os.path.join(dir_outputs, "Rural_" + country_abbrv + "_t2.tif")

Rural_List = [[0,"NODATA"],
              [1, "NODATA"],
              [2, 1],
              [3,"NODATA"],
              [4,"NODATA"],
              ["NODATA","NODATA"]]

Rural_Remap = arcpy.sa.RemapValue(Rural_List)

Rural_Reclass = arcpy.sa.Reclassify(BU_USR_RasterCalc,
                                    "Value",
                                    Rural_Remap,
                                    "DATA");
Rural_Reclass.save(fp_Rural_Reclass)

#Captured Open Space

fp_COS = os.path.join(dir_outputs, "Captured_Open_Space_" + country_abbrv + "_t2.tif")

arcpy.conversion.PolygonToRaster(COS_Polygons_Less_Than_Threshold,
                                 "Value",
                                 fp_COS,
                                 "CELL_CENTER",
                                 "NONE",
                                 cell_size,
                                 "BUILD")

#Fringe Open Space

arcpy.conversion.RasterToPolygon(FringeOpenSpace_Reclass,
                                 r"memory\FOS_Polygon",
                                 "NO_SIMPLIFY",
                                 "Value",
                                 "SINGLE_OUTER_PART",
                                 None)

arcpy.conversion.RasterToPolygon(BU_US_Reclass,
                                 r"memory\US_Polygon",
                                 "NO_SIMPLIFY",
                                 "Value",
                                 "SINGLE_OUTER_PART",
                                 None)

arcpy.analysis.Erase(r"memory\FOS_Polygon",
                     r"memory\US_Polygon",
                     r"memory\FOS_Mask",
                     None)

fp_FOS = os.path.join(dir_outputs, "Fringe_Open_Space_" + country_abbrv + "_t2.tif")

FOS = arcpy.sa.ExtractByMask(FringeOpenSpace_Reclass,
                       r"memory\FOS_Mask");
FOS.save(fp_FOS)

# Urbanized Open Space

arcpy.management.MosaicToNewRaster(fp_FOS + ';' + fp_COS,
                                   dir_outputs,
                                   "Urbanized_Open_Space_" + country_abbrv + "_t2.tif",
                                   'PROJCS["Africa_Albers_Equal_Area_Conic",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",25.0],PARAMETER["Standard_Parallel_1",20.0],PARAMETER["Standard_Parallel_2",-23.0],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]', "8_BIT_UNSIGNED",
                                   None,
                                   1,
                                   "LAST",
                                   "FIRST")

fp_UOS = os.path.join(dir_outputs, "Urbanized_Open_Space_" + country_abbrv + "_t2.tif")

