#----------------------------------------------------------------------------------------------------------------------------
# Copyright 2020 Paul Gosling
#
#   Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the
#   License. You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an
#   "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific
#   language governing permissions and limitations under the License.
#
#----------------------------------------------------------------------------------------------------------------------------
# Title: Automated Projection Selection for ArcGIS (ArcMap Python Add-in)
#
# Purpose: Script provides functionality to run a projection selection process based on user inputs of the purpose of the GIS
#          output/project and a dataset which defines the geographical footprint of interest. Each purpose has a defined set 
#          of candidate projections based on the footprint characteristics. The selection process includes assessment of 
#          distortion in terms of distance, area and shape. Finally, the user can select each of the candidate projection 
#          files and set the data frame coordinate system to graphically view the results of applying that projection.
#
# Author: Paul Gosling
# Version: 1.0
# Published: 3 Jan 2020
#----------------------------------------------------------------------------------------------------------------------------

import arcpy
import pythonaddins
import math
import random
import sys
import os
from arcpy import env
env.overwriteOutput = True

class ProjectedCS(object):

    # Function to initialise parameters for PCS objects created by selection process.
    def __init__(self, pcs_type, proj_name, in_sr, pcs_para1, pcs_para2, pcs_para3, pcs_para4, pcs_para5, pcs_para6):
        self.pcs_type = pcs_type
        self.proj_name = proj_name
        self.in_sr = in_sr
        self.pcs_para1 = pcs_para1
        self.pcs_para2 = pcs_para2
        self.pcs_para3 = pcs_para3
        self.pcs_para4 = pcs_para4
        self.pcs_para5 = pcs_para5
        self.pcs_para6 = pcs_para6

    # Function to create WKT string and SpatialReference object for PCS objects initialised during the selection process.
    def create_sr(self):

    # For each PCS type (categorised by the parameters required to define the PCS) create parameter string using inputs.
        if self.pcs_type == "A":
            parameter_string = 'PARAMETER["Central_Meridian",{0}]'.format(self.pcs_para1)
        elif self.pcs_type == "B":
            parameter_string = 'PARAMETER["Central_Meridian",{0}],PARAMETER["Latitude_Of_Origin",{1}]'.format(self.pcs_para1, self.pcs_para2)
        elif self.pcs_type == "C":
            parameter_string = 'PARAMETER["Central_Meridian",{0}],PARAMETER["Standard_Parallel_1",{1}]'.format(self.pcs_para1, self.pcs_para2)
        elif self.pcs_type == "D":
            parameter_string = 'PARAMETER["Central_Meridian",{0}],PARAMETER["Standard_Parallel_1",{1}],PARAMETER["Standard_Parallel_2",{2}],'\
                               'PARAMETER["Latitude_Of_Origin",{3}]'.format(self.pcs_para1, self.pcs_para2, self.pcs_para3, self.pcs_para4)
        elif self.pcs_type == "E":
            parameter_string = 'PARAMETER["Longitude_Of_Center",{0}],PARAMETER["Latitude_Of_Center",{1}]'.format(self.pcs_para1, self.pcs_para2)
        elif self.pcs_type == "F":
            parameter_string = 'PARAMETER["Central_Meridian",{0}],PARAMETER["Scale_Factor",{1}],PARAMETER["Latitude_Of_Origin",{2}]'\
                               .format(self.pcs_para1, self.pcs_para2, self.pcs_para3)
        elif self.pcs_type == "G":
            parameter_string = 'PARAMETER["Central_Meridian",{0}],PARAMETER["Option",{1}]'.format(self.pcs_para1, self.pcs_para2)
        elif self.pcs_type == "H":
            parameter_string = 'PARAMETER["Latitude_Of_1st_Point",{0}],PARAMETER["Latitude_Of_2nd_Point",{1}],PARAMETER["Scale_Factor",{2}],'\
                               'PARAMETER["Longitude_Of_1st_Point",{3}],PARAMETER["Longitude_Of_2nd_Point",{4}],PARAMETER["Latitude_Of_Center",{5}]'\
                               .format(self.pcs_para1, self.pcs_para2, self.pcs_para3, self.pcs_para4, self.pcs_para5, self.pcs_para6)
        elif self.pcs_type == "I":
            parameter_string = 'PARAMETER["Central_Meridian",{0}],PARAMETER["Standard_Parallel_1",{1}],PARAMETER["Standard_Parallel_2",{2}],'\
                               'PARAMETER["Scale_Factor",{3}],PARAMETER["Latitude_Of_Origin",{4}]'\
                               .format(self.pcs_para1, self.pcs_para2, self.pcs_para3, self.pcs_para4, self.pcs_para5)
        elif self.pcs_type == "J":
            parameter_string = 'PARAMETER["Central_Meridian",{0}],PARAMETER["Central_Parallel",{1}]'.format(self.pcs_para1, self.pcs_para2)
        elif self.pcs_type == "K":
            parameter_string = 'PARAMETER["Latitude_Of_1st_Point",{0}],PARAMETER["Latitude_Of_2nd_Point",{1}],PARAMETER["Longitude_Of_1st_Point",{2}],'\
                               'PARAMETER["Longitude_Of_2nd_Point",{3}]'.format(self.pcs_para1, self.pcs_para2, self.pcs_para3, self.pcs_para4)
        elif self.pcs_type == "L":
            parameter_string = 'PARAMETER["Longitude_Of_Center",{0}],PARAMETER["Latitude_Of_Center",{1}],PARAMETER["Height",{2}]'\
                               .format(self.pcs_para1, self.pcs_para2, self.pcs_para3)

    # Create WKT string using projection name, GCS string derived from input dataset spatial reference, and parameter string.
    # Create SpatialReference object by loading from WKT string.
        input_sr_string = self.in_sr.exportToString()
        input_GCS_string = input_sr_string.split(";")[0]
        input_GCS_string = input_GCS_string.replace("'", '"')
        pcs_name = '"Custom_' + self.proj_name + '"'
        wkt_string = 'PROJCS[{0},{1},PROJECTION["{2}"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],{3},UNIT["Meter",1.0]]'\
                     .format(pcs_name, input_GCS_string, self.proj_name, parameter_string)
        spatial_ref = arcpy.SpatialReference()
        spatial_ref.loadFromString(wkt_string)

        return spatial_ref, wkt_string

# Function to determine input type grouping where function uses polygon to define footprint of interest, and returns footprint characteristics to be used in PCS definitions.
def footprint_type(in_lyr, out_path, in_sr, map, df):

    # Read input feature geometry values including GCS footprint centroid lat/long as initial centroid approximation.
    cursor = arcpy.da.SearchCursor(in_lyr, ["SHAPE@"])
    for row in cursor:
        lat_max, lat_min = row[0].extent.YMax, row[0].extent.YMin
        lon_max, lon_min = row[0].extent.XMax, row[0].extent.XMin
        centre = row[0].trueCentroid
        lon_centre_gcs, lat_centre_gcs  = centre.X, centre.Y
    del row, cursor

    # Use approximate centroid to project input feature using Azimuthal Equidistant projection based on that location.
    env.workspace = out_path
    equidistant_fc = os.path.join(out_path, "footprint_centroid")
    pcs_equidistant = ProjectedCS("B", "Azimuthal_Equidistant", in_sr, lon_centre_gcs, lat_centre_gcs, "", "", "", "")
    centroid_sr, equidistant_string = pcs_equidistant.create_sr()
    arcpy.Project_management(in_lyr, equidistant_fc, centroid_sr)

    # Determine centroid of projected feature and create new feature class containing that point using Azimuthal Equidistant projection.
    cursor = arcpy.da.SearchCursor(equidistant_fc, ["SHAPE@TRUECENTROID"])
    for row in cursor:
        x_centre, y_centre = row[0]
    del row, cursor
    arcpy.CreateFeatureclass_management(out_path, "pt_centroid", "POINT", spatial_reference = centroid_sr)
    centroid_fc = os.path.join(out_path, "pt_centroid")
    centroid = (x_centre, y_centre)
    cursor = arcpy.da.InsertCursor(centroid_fc, ["SHAPE@XY"])
    cursor.insertRow([centroid])
    del cursor

    # Project to GCS and read lat/long as better approximation of geographic centre of input feature.
    centroid_fc_gcs = os.path.join(out_path, "pt_centroid_gcs")
    arcpy.Project_management(centroid_fc, centroid_fc_gcs, in_sr)
    cursor = arcpy.da.SearchCursor(centroid_fc_gcs, ["SHAPE@XY"])
    for row in cursor:
        lon_centre, lat_centre = row[0]
    del row, cursor

    # Determine input footprint size by adding geodesic area attribute to input feature.
    arcpy.AddGeometryAttributes_management(in_lyr, "AREA_GEODESIC", Area_Unit = "SQUARE_KILOMETERS")
    cursor = arcpy.da.SearchCursor(in_lyr, ["AREA_GEO"])
    for row in cursor:
        in_size = row[0]
    del row, cursor

    # Calculate total area size of input ellipsoid. Designed to work for any ellipsoid including other celestial bodies, assuming use of the same lat/long system.
    coordinates = [(-180.0,-90.0),(-180.0,90.0),(180.0,90.0),(180.0,-90.0)]
    arcpy.CreateFeatureclass_management(out_path, "ellipsoid_polygon", "POLYGON", spatial_reference = in_sr)
    ellipsoid_fc = os.path.join(out_path, "ellipsoid_polygon")
    with arcpy.da.InsertCursor(ellipsoid_fc, ['SHAPE@']) as cursor:
        cursor.insertRow([coordinates])
    del cursor, coordinates
    arcpy.AddGeometryAttributes_management(ellipsoid_fc, "AREA_GEODESIC", Area_Unit = "SQUARE_KILOMETERS")
    cursor = arcpy.da.SearchCursor(ellipsoid_fc, ["AREA_GEO"])
    for row in cursor:
        ellipsoid_size = row[0]
    del row, cursor
    half_size = ellipsoid_size / 2
    quarter_size = ellipsoid_size / 4

    # Use area size to define World input footprint type. Lower limit is half the total area size calculated for input ellipsoid.
    # Reset lon/lat of centre to avoid rounding errors in centroid calculation.
    if in_size >= half_size:
        in_type = "World"
        lon_centre, lat_centre = 0.0, 0.0

    # Use area size to define Hemisphere input footprint type, between quarter and half the total area size calculated for input ellipsoid.
    elif in_size > quarter_size and in_size < half_size:
        in_type = "Hemisphere"

    # If footprint geometry encompasses single Pole or centroid in Polar region then set appropriate type, reset latitude to Pole and longitude to Greenwich Meridian
    # if footprint encompasses the entire polar region.
    elif lat_max > 89 or lat_centre > 70:
        in_type = "North Pole"
        lat_centre = 90.0
        if lon_max > 179 and lon_min < -179:
            lon_centre = 0.0
    elif lat_min < -89 or lat_centre < -70:
        in_type = "South Pole"
        lat_centre = -90.0
        if lon_max > 179 and lon_min < -179:
            lon_centre = 0.0

    # If footprint geometry crosses Equator or centroid is within certain latitude then set appropriate type and reset latitude centre to Equator.
    elif (lat_max > 0 and lat_min < 0) or abs(lat_centre) < 15:
        in_type = "Equatorial"
        lat_centre = 0.0

    # By process of elimination all remaining footprints are at Middle latitudes so set appropriate type.
    else:
        in_type = "Middle"

    # Remove all layers generated by function from data frame.
    layer_list = arcpy.mapping.ListLayers(map, "", df)
    for layer in layer_list:
        if layer.name == "footprint_centroid" or layer.name == "pt_centroid" or layer.name == "pt_centroid_gcs" or layer.name == "ellipsoid_polygon":
            arcpy.mapping.RemoveLayer(df, layer)

    return in_type, in_size, lon_centre, lat_centre, lat_max, lat_min

# Function to define footprint extent characteristics grouping.
def footprint_extent(in_lyr, out_path, in_sr, in_lon, in_lat, map, df):

    # Create new temp file and project input features using custom Azimuthal Equidistant projection at centroid of input features.
    env.workspace = out_path
    extent_fc = os.path.join(out_path,"extent_fc")
    pcs_extent = ProjectedCS("B", "Azimuthal_Equidistant", in_sr, in_lon, in_lat, "", "", "", "")
    extent_sr, extent_string = pcs_extent.create_sr()
    arcpy.Project_management(in_lyr, extent_fc, extent_sr)

    # Determine minimum and maximum extents of projected feature. FE/FN zero at centre of footprint so min_x, min_y will be negative.
    desc_extent = arcpy.Describe(extent_fc)
    eqd_extent = desc_extent.extent
    min_x, min_y, max_x, max_y = eqd_extent.XMin, eqd_extent.YMin, eqd_extent.XMax, eqd_extent.YMax

    # Calculate total XY extent containing all features and determine extent group from ratio calculation.
    delta_x = max_x - min_x
    delta_y = max_y - min_y
    extent_ratio = delta_x / delta_y
    if extent_ratio > 1.25:
        extent_group = "EastWest"
    elif extent_ratio < 0.8:
        extent_group = "NorthSouth"
    else:
        extent_group = "Equal"

    # Remove all layers generated by function from data frame.
    layer_list = arcpy.mapping.ListLayers(map, "", df)
    for layer in layer_list:
        if layer.name == "extent_fc":
            arcpy.mapping.RemoveLayer(df, layer)

    return extent_group, extent_ratio

# Function to determine input type grouping where function uses set of individual point/line/polygon features to define footprint of interest.
# Returns footprint characteristics to be used in PCS definitions.
def location_type(in_purpose, in_lyr, out_path, in_sr, map, df):

    # Make feature layer from input feature class for Select Layer functions below.
    env.workspace = out_path
    arcpy.MakeFeatureLayer_management(in_lyr, "lyr")

    # Determine approximate mean centre from GCS input features.
    mc_approx_gcs = os.path.join(out_path, "mc_approx_gcs")
    arcpy.MeanCenter_stats(in_lyr, mc_approx_gcs)
    cursor = arcpy.da.SearchCursor(mc_approx_gcs, ["SHAPE@XY"])
    for row in cursor:
        mc_lon, mc_lat = row[0]
    del row, cursor

    # Use approximate mean centre to project input features using Azimuthal Equidistant projection based on that location.
    pcs_eqd_approx = ProjectedCS("B", "Azimuthal_Equidistant", in_sr, mc_lon, mc_lat, "", "", "", "")
    mc_approx_sr, eqd_approx_string = pcs_eqd_approx.create_sr()
    eqd_approx_fc = os.path.join(out_path, "eqd_approx_fc")
    arcpy.Project_management(in_lyr, eqd_approx_fc, mc_approx_sr)

    # Recalculate mean centre using projected features and determine GCS lat/lon as better approximation of geographic centre of input features.
    mc_eqd = os.path.join(out_path, "mean_center_eqd")
    arcpy.MeanCenter_stats(eqd_approx_fc, mc_eqd)
    mc_gcs = os.path.join(out_path, "mean_center_gcs")
    arcpy.Project_management(mc_eqd, mc_gcs, in_sr)
    cursor = arcpy.da.SearchCursor(mc_gcs, ["SHAPE@XY"])
    for row in cursor:
        lon_centre, lat_centre = row[0]
    del row, cursor

    # Determine maximum latitude extents for GCS input features for Pole type definition and conic projection standard parallel determination.
    arcpy.CreateFeatureclass_management(out_path, "input_extent", "POLYGON", spatial_reference = in_sr)
    extent_fc = os.path.join(out_path, "input_extent")
    cursor = arcpy.da.InsertCursor(extent_fc, ["SHAPE@"])
    desc_layer = arcpy.Describe(in_lyr)
    layer_extent = desc_layer.extent
    extent_polygon = layer_extent.polygon
    cursor.insertRow([extent_polygon])
    del cursor
    cursor = arcpy.da.SearchCursor(extent_fc, ["SHAPE@"])
    for row in cursor:
        lat_max, lat_min = row[0].extent.YMax, row[0].extent.YMin
    del row, cursor

    # Use Select by Location to define maximum limit for Geospatial Analysis Distance purpose. Buffer limit at 500Km is gut feeling and requires further study.
    # Equatorial type not necessary for this function.
    if in_purpose == "Geospatial Analysis (Distance)":
        buffer_analysis = os.path.join(out_path, "buffer_analysis")
        arcpy.Buffer_analysis(mc_gcs, buffer_analysis, "500000 Meters", method = "GEODESIC")
        arcpy.SelectLayerByLocation_management("lyr", "COMPLETELY_WITHIN", buffer_analysis, invert_spatial_relationship = "INVERT")
        count_analysis = int(arcpy.GetCount_management("lyr")[0])
        if count_analysis == 0:   # Selection inverted so if all features within limit then count is zero.
            in_type = "Analysis"
        else:
            in_type = "World"

    # Create line dataset to measure geodesic distance on the Equator from zero to 180E. Calculate limits based on this distance to enable determination of whether all
    # input features fall within buffers at these sizes.
    else:
        coordinates = [(0.0,0.0),(180.0,0.0)]
        arcpy.CreateFeatureclass_management(out_path, "eq_distance", "POLYLINE", spatial_reference = in_sr)
        eqdist_fc = os.path.join(out_path, "eq_distance")
        with arcpy.da.InsertCursor(eqdist_fc, ['SHAPE@']) as cursor:
            cursor.insertRow([coordinates])
        del cursor, coordinates
        arcpy.AddGeometryAttributes_management(eqdist_fc, "LENGTH_GEODESIC", Length_Unit = "KILOMETERS")
        cursor = arcpy.da.SearchCursor(eqdist_fc, ["LENGTH_GEO"])
        for row in cursor:
            ellipsoid_dist = row[0]
        del row, cursor
        half_dist = ellipsoid_dist / 2
        quarter_dist = ellipsoid_dist / 4

    # Use Select by Location to define footprint types for various functions which use input features to define footprint of interest. Set Middle when all features fall within buffer at 1/4 the
    # equator distance, Hemisphere within 1/2 distance, World otherwise. Reset central lon/lat for World type.
        buffer_middle = os.path.join(out_path,"buffer_middle")
        middle_dist = str(quarter_dist) + " Kilometers"
        arcpy.Buffer_analysis(mc_gcs, buffer_middle, middle_dist, method = "GEODESIC")
        arcpy.SelectLayerByLocation_management("lyr", "COMPLETELY_WITHIN", buffer_middle, invert_spatial_relationship = "INVERT")
        count_middle = int(arcpy.GetCount_management("lyr")[0])
        buffer_hemisphere = os.path.join(out_path,"buffer_hemisphere")
        hemisphere_dist = str(half_dist) + " Kilometers"
        arcpy.Buffer_analysis(mc_gcs, buffer_hemisphere, hemisphere_dist, method = "GEODESIC")
        arcpy.SelectLayerByLocation_management("lyr", "COMPLETELY_WITHIN", buffer_hemisphere, invert_spatial_relationship = "INVERT")
        count_hemisphere = int(arcpy.GetCount_management("lyr")[0])
        if count_middle == 0:   # Selection inverted so if all features within limit then count is zero.
            in_type = "Middle"
        elif count_hemisphere == 0:   # Selection inverted so if all features within limit then count is zero.
            in_type = "Hemisphere"
        else:
            in_type = "World"
            lon_centre, lat_centre = 0.0, 0.0

    # For Middle size footprints, if footprint extent geometry crosses Equator or centroid is within certain latitude then reset type and set latitude centre to Equator.
        if in_type == "Middle":
            if (lat_max > 0 and lat_min < 0) or abs(lat_centre) < 15:
                in_type = "Equatorial"
                lat_centre = 0.0

    # For non-World size footprints, if footprint extends to Pole or centroid in Polar region then reset type, set latitude centre to Pole and longitude to Greenwich Meridian
    # if footprint encompasses the entire polar region.
    if in_type != "World":
        if lat_max > 89 or lat_centre > 70:
            in_type = "North Pole"
            lat_centre = 90.0
            if lon_max > 179 and lon_min < -179:
                lon_centre = 0.0
        elif lat_min < -89 or lat_centre < -70:
            in_type = "South Pole"
            lat_centre = -90.0
            if lon_max > 179 and lon_min < -179:
                lon_centre = 0.0

    # Remove all layers generated by function from data frame.
    layer_list = arcpy.mapping.ListLayers(map, "", df)
    for layer in layer_list:
        if layer.name == "lyr" or layer.name == "input_extent" or layer.name == "mc_approx_gcs" or layer.name == "eqd_approx_fc"\
           or layer.name == "mean_center_eqd" or layer.name == "mean_center_gcs" or layer.name == "buffer_analysis"\
           or layer.name == "eq_distance" or layer.name == "buffer_middle" or layer.name == "buffer_hemisphere":
            arcpy.mapping.RemoveLayer(df, layer)

    return in_type, lon_centre, lat_centre, lat_max, lat_min

# Function to calculate standard parallels for conic projections
def conic_sp(in_latmax, in_latmin, in_k):

    # Determine standard parallels for input footprint using min/max latitudes read from geometry in footprint_type function, and k constant value defining extent type.
    # Final values are rounded to nearest 0.25 (15 mins) to ensure sensible standard parallel values for quoting on map publication
    delta_lat = in_latmax - in_latmin
    adjust_lat = delta_lat / in_k
    conic_sp1 = (round((in_latmax - adjust_lat) / 0.25)) * 0.25
    conic_sp2 = (round((in_latmin + adjust_lat) / 0.25)) * 0.25

    return conic_sp1, conic_sp2

# Function to create extent polygon for footprint defined by raster dataset or raster catalog.
# In this case the data fills the footprint described by the extent object boundaries. This function would not be appropriate for feature classes.
def raster_polygon(in_lyr, out_path, in_sr, map, df):

    # Make copy of input raster dataset in GDB to ensure extent polygon function below uses the correct spatial reference.
    env.workspace = out_path
    in_copy = os.path.join(out_path, "raster_copy")
    arcpy.CopyRaster_management(in_lyr, in_copy)

    # Create extent feature class and insert polygon based on input layer extent object.
    arcpy.CreateFeatureclass_management(out_path, "raster_extent", "POLYGON", spatial_reference = in_sr)
    extent_fc = os.path.join(out_path, "raster_extent")
    cursor = arcpy.da.InsertCursor(extent_fc,["SHAPE@"])
    desc_layer = arcpy.Describe(in_copy)
    layer_extent = desc_layer.extent
    extent_polygon = layer_extent.polygon
    cursor.insertRow([extent_polygon])
    del cursor

    # Remove all layers generated by function from data frame.
    layer_list = arcpy.mapping.ListLayers(map, "", df)
    for layer in layer_list:
        if layer.name == "raster_copy" or layer.name == "raster_extent":
            arcpy.mapping.RemoveLayer(df, layer)

    return extent_fc

# Function to read input features and determine whether all start points identical.
def nav_condition(in_lyr, feat_count):

    cursor = arcpy.da.SearchCursor(in_lyr, ["SHAPE@"])
    i = 0
    for row in cursor:
        lon, lat = row[0].firstPoint.X, row[0].firstPoint.Y
        if i == 0:
            start_lon, start_lat = lon, lat
            i = i = 1
        else:
            if lon == start_lon and lat == start_lat:
                i = i + 1
    del row, cursor
    if i == feat_count:
        nav_single_pt = "TRUE"
    else:
        nav_single_pt = "FALSE"

    return nav_single_pt, start_lon, start_lat

# Function to carry out projection selection based on user inputs. For each map purpose run functions to determine footprint characteristics.
# Create ProjectedCS class objects for selections based on conceptual diagram for that purpose.
def projection_selection(out_gdb, in_purpose, in_layer, in_type, in_count, in_sr, in_map, in_df):

# Create list to hold ProjectedCS class objects and initialise variables ensuring these have value when not defined by map purpose requirements.
    projections, notes = [], []
    geo_size_str, geo_sp_str, geo_extent_str = "N/A", "N/A", "N/A"
    geo_lon, geo_lat, geo_lat_max, geo_lat_min, geo_k = 0.0, 0.0, 0.0, 0.0, 0
    line_details, point_details = (), ()

# Create notes to be added to list when certain criteria met, later adding to results text file to provide assistance to user.
    conic_note = "  - Conic (Albers, Lambert Conformal, Equidistant) - if distortion assessment is employed and one of these projections is the number \n"\
                 "    one rank for the initial standard parallel and K value settings, then all remaining K values between 3 and 7 will be used to \n"\
                 "    create and assess additional standard parallel options.\n"
    sf_note = "  - Projection using Scale Factor parameter (Transverse Mercator or Stereographic) - if distortion assessment is employed and one of \n"\
              "    these projections is the number one rank for the initial SF setting of 1, then SF values between 0.9996 and 0.9999 (step 0.0001) \n"\
              "    will be used to create and assess additional options.\n"
    eqd_cyl_note = "  - Equidistant Cylindrical - standard parallel is set to 45.0 to provide an alternative to Plate Carree. The user is advised \n"\
                   "    that this parameter can be modified using the Coordinate System properties to provide different looking outputs as required.\n"
    winkel_tripel_note = "  - Winkel Tripel - standard parallel is set to 40.0 as per the Times Atlas version. The user is advised that this \n"\
                         "    parameter can be modified using the Coordinate System properties to provide different looking outputs as required.\n"
    vnsp_note = "  - Vertical Near-Side Perspective - height parameter is set to 5,000 or 10,000 Km depending on the footprint category. The user is advised \n"\
                "    that this parameter can be modified using the Coordinate System properties to provide different looking outputs as required.\n"
    world_note = "  - World projections - for the World footprint type note that the central longitude has been set to the Greenwich Meridian. The user is advised \n"\
                 "    that this parameter can be modified using the Coordinate System properties to re-centre the map as necessary. For example, 10E is often \n"\
                 "    used to ensure that Russia is shown as a single entity rather than being split by the 180 meridian.\n"
    grid_note = "  - Footprints smaller than World or Hemisphere size - the user is advised that in many situations an additional sensible choice can be \n"\
                "    the reference system used by the local or national mapping agency, or in global use such as the UTM system. This is especially true \n"\
                "    if certain datasets to be used in the project are already referenced to that system. Note that these systems mostly use conformal \n"\
                "    map projections such as Transverse Mercator and Lambert Conformal Conic.\n"
    geo_analysis_area_note = "  - When selecting Geospatial Analysis (Area) purpose the projection selection is directed to use either Thematic Vector or Thematic Raster \n"\
                             "    selection diagrams depending on input data type. If necessary, the Minimum Bounding Geometry tool is used to create a single polygon \n"\
                             "    input for Thematic Vector, so please note this tool can lead to inappropriate results. Please review the footprint_mbg feature class and if \n"\
                             "    necessary create your own single polygon which encompasses the polygon features then select Thematic Vector purpose.\n"
    equal_area_note = "  - When selecting Thematic Vector / Raster or Geospatial Analysis (Area) purposes the distortion assessment is based only on area distortion. \n"\
                      "    In all cases the combined index result is almost exactly zero and the ranking order is effectively arbitrary. Therefore, the user should view each \n"\
                      "    option graphically before making their selection decision using qualitative reasons.\n"

# For Geo Analysis Area purpose, reset to relevant Thematic purpose depending on data type.
    if in_purpose == "Geospatial Analysis (Area)":
        if in_type == "RasterDataset":
            in_purpose = "Thematic Raster"
        else:
            in_purpose = "Thematic Vector"
            if in_count != 1:  # If input has multiple polygons create single input using Minimum Bounding Geometry tool.
                out_layer = os.path.join(out_gdb, "footprint_mbg")
                arcpy.MinimumBoundingGeometry_management(in_layer, out_layer, "CONVEX_HULL", "ALL")
                in_layer = out_layer
            notes.append(geo_analysis_area_note)

# General Reference: standard mapping type which focusses on the location of geographical features, e.g. atlas maps, topographical maps etc.
    if in_purpose == "General Reference":
        geo_type, geo_size, geo_lon, geo_lat, geo_lat_max, geo_lat_min = footprint_type(in_layer, out_gdb, in_sr, in_map, in_df)
        geo_size_str = str(round(geo_size, 2))
        if geo_type == "World":
            projections.append(ProjectedCS("C", "Winkel_Tripel", in_sr, geo_lon, 40.0, "", "", "", ""))
            notes.append(winkel_tripel_note)
            projections.append(ProjectedCS("A", "Robinson", in_sr, geo_lon, "", "", "", "", ""))
            projections.append(ProjectedCS("A", "NaturalEarth", in_sr, geo_lon, "", "", "", "", ""))
            notes.append(world_note)
        elif geo_type == "Hemisphere":
            projections.append(ProjectedCS("E", "Orthographic", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("B", "Lambert_Azimuthal_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("L", "Vertical_Near_Side_Perspective", in_sr, geo_lon, geo_lat, 10000000.0, "", "", ""))
            notes.append(vnsp_note)
        elif geo_type == "North Pole" or geo_type == "South Pole":
            projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
            notes.append(sf_note)
            projections.append(ProjectedCS("B", "Lambert_Azimuthal_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
        elif geo_type == "Equatorial":
            projections.append(ProjectedCS("C", "Mercator", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("A", "Sinusoidal", in_sr, geo_lon, "", "", "", "", ""))
        elif geo_type == "Middle":
            geo_extent_group, geo_extent_ratio = footprint_extent(in_layer, out_gdb, in_sr, geo_lon, geo_lat, in_map, in_df)
            geo_extent_str = "{0} (Extent X/Y ratio = {1:.3f})".format(geo_extent_group, geo_extent_ratio)
            if geo_extent_group == "EastWest":
                geo_k = 7  # k = 7 for EW extent as per Maling(1992)
                geo_sp1, geo_sp2 = conic_sp(geo_lat_max, geo_lat_min, geo_k)
                geo_sp_str = "\n   Standard parallels: SP1 = {0}, SP2 = {1} (K = {2})".format(str(geo_sp1), str(geo_sp2), str(geo_k))
                projections.append(ProjectedCS("D", "Equidistant_Conic", in_sr, geo_lon, geo_sp1, geo_sp2, geo_lat, "", ""))
                projections.append(ProjectedCS("I", "Lambert_Conformal_Conic", in_sr, geo_lon, geo_sp1, geo_sp2, 1.0, geo_lat, ""))
                notes.append(conic_note)
            elif geo_extent_group == "NorthSouth":
                projections.append(ProjectedCS("F", "Transverse_Mercator", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
                notes.append(sf_note)
                projections.append(ProjectedCS("C", "Bonne", in_sr, geo_lon, geo_lat, "", "", "", ""))
            elif geo_extent_group == "Equal":
                projections.append(ProjectedCS("B", "Lambert_Azimuthal_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
                projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))

# Thematic Vector: standard mapping type which focusses on the geographical distribution of a phenomena rather than the location, e.g. statistical maps.
# Note: code also used for Geospatial Analysis (Area) purpose. Some analysis techniques are based on an accurate measurement of the area of an object,
# e.g. population density.
    elif in_purpose == "Thematic Vector":
        notes.append(equal_area_note)
        geo_type, geo_size, geo_lon, geo_lat, geo_lat_max, geo_lat_min = footprint_type(in_layer, out_gdb, in_sr, in_map, in_df)
        geo_size_str = str(round(geo_size, 2))
        if geo_type == "World":
            projections.append(ProjectedCS("A", "Mollweide", in_sr, geo_lon, "", "", "", "", ""))
            projections.append(ProjectedCS("G", "Goode_Homolosine", in_sr, geo_lon, 0.0, "", "", "", ""))
            projections.append(ProjectedCS("G", "Goode_Homolosine", in_sr, geo_lon, 1.0, "", "", "", ""))
            projections.append(ProjectedCS("G", "Goode_Homolosine", in_sr, geo_lon, 2.0, "", "", "", ""))
            projections.append(ProjectedCS("A", "Hammer_Aitoff", in_sr, geo_lon, "", "", "", "", ""))
            projections.append(ProjectedCS("A", "Sinusoidal", in_sr, geo_lon, "", "", "", "", ""))
            notes.append(world_note)
        elif geo_type == "Hemisphere":
            projections.append(ProjectedCS("B", "Lambert_Azimuthal_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
        elif geo_type == "North Pole" or geo_type == "South Pole":
            projections.append(ProjectedCS("B", "Lambert_Azimuthal_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
        elif geo_type == "Equatorial":
            projections.append(ProjectedCS("A", "Sinusoidal", in_sr, geo_lon, "", "", "", "", ""))
            projections.append(ProjectedCS("C", "Cylindrical_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
        elif geo_type == "Middle":
            geo_extent_group, geo_extent_ratio = footprint_extent(in_layer, out_gdb, in_sr, geo_lon, geo_lat, in_map, in_df)
            geo_extent_str = "{0} (Extent X/Y ratio = {1:.3f})".format(geo_extent_group, geo_extent_ratio)
            if geo_extent_group == "EastWest":
                geo_k = 7  # k = 7 for EW extent as per Maling(1992)
                geo_sp1, geo_sp2 = conic_sp(geo_lat_max, geo_lat_min, geo_k)
                geo_sp_str = "\n   Standard parallels: SP1 = {0}, SP2 = {1} (K = {2})".format(str(geo_sp1), str(geo_sp2), str(geo_k))
                projections.append(ProjectedCS("D", "Albers", in_sr, geo_lon, geo_sp1, geo_sp2, geo_lat, "", ""))
                notes.append(conic_note)
                projections.append(ProjectedCS("C", "Bonne", in_sr, geo_lon, geo_lat, "", "", "", ""))
            elif geo_extent_group == "NorthSouth":
                projections.append(ProjectedCS("F", "Transverse_Cylindrical_Equal_Area", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
                notes.append(sf_note)
                projections.append(ProjectedCS("C", "Bonne", in_sr, geo_lon, geo_lat, "", "", "", ""))
            elif geo_extent_group == "Equal":
                projections.append(ProjectedCS("B", "Lambert_Azimuthal_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
                projections.append(ProjectedCS("C", "Bonne", in_sr, geo_lon, geo_lat, "", "", "", ""))

# Thematic Raster: standard mapping type which focusses on the geographical distribution of a phenomena rather than the location, e.g. land use / land cover.
# Note: code also used for Geospatial Analysis (Area) purpose. Some analysis techniques are based on an accurate measurement of the area of an object,
# e.g. population density.
    elif in_purpose == "Thematic Raster":
        notes.append(equal_area_note)
        geo_type, geo_size, geo_lon, geo_lat, geo_lat_max, geo_lat_min = footprint_type(in_layer, out_gdb, in_sr, in_map, in_df)
        geo_size_str = str(round(geo_size, 2))
        if geo_type == "World":
            projections.append(ProjectedCS("G", "Goode_Homolosine", in_sr, geo_lon, 0.0, "", "", "", ""))
            projections.append(ProjectedCS("G", "Goode_Homolosine", in_sr, geo_lon, 1.0, "", "", "", ""))
            projections.append(ProjectedCS("G", "Goode_Homolosine", in_sr, geo_lon, 2.0, "", "", "", ""))
            projections.append(ProjectedCS("A", "Sinusoidal", in_sr, geo_lon, "", "", "", "", ""))
            projections.append(ProjectedCS("A", "Mollweide", in_sr, geo_lon, "", "", "", "", ""))
            projections.append(ProjectedCS("B", "Wagner_IV", in_sr, 0.0, 0.0, "", "", "", ""))
            projections.append(ProjectedCS("B", "Wagner_VII", in_sr, 0.0, 0.0, "", "", "", ""))
            notes.append(world_note)
        elif geo_type == "Hemisphere":
            projections.append(ProjectedCS("B", "Lambert_Azimuthal_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("A", "Sinusoidal", in_sr, geo_lon, "", "", "", "", ""))
        elif geo_type == "North Pole" or geo_type == "South Pole":
            projections.append(ProjectedCS("B", "Lambert_Azimuthal_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
        elif geo_type == "Equatorial" or geo_type == "Middle":
            projections.append(ProjectedCS("B", "Lambert_Azimuthal_Equal_Area", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("A", "Sinusoidal", in_sr, geo_lon, "", "", "", "", ""))
            projections.append(ProjectedCS("C", "Bonne", in_sr, geo_lon, geo_lat, "", "", "", ""))

# Geospatial Analysis (Distance): numerous analysis techniques and GIS tools require Cartesian XY inputs and are based on the calculation of Euclidean
# distance between objects, e.g. Spatial Distribution, Nearest Neighbour, Spatial Autocorrelation, Clustering and geostatistical techniques such as
# Inverse Distance Weighting (IDW).
    elif in_purpose == "Geospatial Analysis (Distance)":
        geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr, in_map, in_df)
        if geo_type == "World":
            pythonaddins.MessageBox("ArcGIS analysis techniques are not applicable for an area of this size. Please select a different input dataset.",\
                                    "Warning", 0)
            sys.exit()
        elif geo_type == "North Pole" or geo_type == "South Pole":
            projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
            notes.append(sf_note)
        else:
            geo_extent_group, geo_extent_ratio = footprint_extent(in_layer, out_gdb, in_sr, geo_lon, geo_lat, in_map, in_df)
            geo_extent_str = "{0} (Extent X/Y ratio = {1:.3f})".format(geo_extent_group, geo_extent_ratio)
            if geo_extent_group == "EastWest":
                geo_k = 7  # k = 7 for EW extent as per Maling(1992)
                geo_sp1, geo_sp2 = conic_sp(geo_lat_max, geo_lat_min, geo_k)
                geo_sp_str = "\n   Standard parallels: SP1 = {0}, SP2 = {1} (K = {2})".format(str(geo_sp1), str(geo_sp2), str(geo_k))
                projections.append(ProjectedCS("D", "Equidistant_Conic", in_sr, geo_lon, geo_sp1, geo_sp2, geo_lat, "", ""))
                projections.append(ProjectedCS("I", "Lambert_Conformal_Conic", in_sr, geo_lon, geo_sp1, geo_sp2, 1.0, geo_lat, ""))
                notes.append(conic_note)
            elif geo_extent_group == "NorthSouth":
                projections.append(ProjectedCS("F", "Transverse_Mercator", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
                notes.append(sf_note)
                projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
            elif geo_extent_group == "Equal":
                projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
                projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
                notes.append(sf_note)

# Navigation Routes (Geodesic): Navigation based on travelling the shortest distance between two locations, also referred to as Great Circle navigation.
    elif in_purpose == "Navigation Routes (Geodesic)":
        single_start, geo_lon, geo_lat = nav_condition(in_layer, in_count)
        if in_count == 1:
            geo_type = "1 Line"
            # AddGeomAttributes used as line feature geometry object does not include mid value required for projection definition.
            arcpy.AddGeometryAttributes_management(in_layer, "LENGTH_GEODESIC; LINE_START_MID_END")
            cursor = arcpy.da.SearchCursor(in_layer,["LENGTH_GEO", "START_X", "START_Y", "END_X", "END_Y", "MID_X", "MID_Y"])
            for row in cursor:
                length_geo = row[0]
                lon1, lat1, lon2, lat2 = row[1], row[2], row[3], row[4]
                if lon2 > 180:   # Reset longitude for geodesic line features created using the XY to Line function which cross 180.
                    lon2 = lon2 - 360
                centre_lon, centre_lat = row[5], row[6]
            del row, cursor
            projections.append(ProjectedCS("H", "Hotine_Oblique_Mercator_Two_Point_Center", in_sr, lat1, lat2, 1.0, lon1, lon2, centre_lat))
            projections.append(ProjectedCS("K", "Two_Point_Equidistant", in_sr, lat1, lat2, lon1, lon2, "", ""))
            projections.append(ProjectedCS("E", "Gnomonic", in_sr, centre_lon, centre_lat, "", "", "", ""))
            projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, centre_lon, centre_lat, "", "", "", ""))
            line_details = (lon1, lat1, centre_lon, centre_lat, lon2, lat2, length_geo)
        elif single_start == "TRUE":   # Checks whether all lines have same start point and therefore azimuthal projection from that point more logical selection.
            geo_type = "Single Start Point"
            projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("E", "Gnomonic", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("E", "Orthographic", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
            notes.append(sf_note)
        else:
            geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr, in_map, in_df)
            if geo_type == "World":
                projections.append(ProjectedCS("C", "Winkel_Tripel", in_sr, geo_lon, 40.0, "", "", "", ""))
                notes.append(winkel_tripel_note)
                projections.append(ProjectedCS("A", "Robinson", in_sr, geo_lon, "", "", "", "", ""))
                projections.append(ProjectedCS("A", "NaturalEarth", in_sr, geo_lon, "", "", "", "", ""))
                notes.append(world_note)
            elif geo_type == "Hemisphere":
                projections.append(ProjectedCS("E", "Orthographic", in_sr, geo_lon, geo_lat, "", "", "", ""))
                projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
                notes.append(sf_note)
                projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
            elif geo_type == "North Pole" or geo_type == "South Pole":
                projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
                notes.append(sf_note)
                projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
                projections.append(ProjectedCS("E", "Gnomonic", in_sr, geo_lon, geo_lat, "", "", "", ""))
            elif geo_type == "Equatorial":
                projections.append(ProjectedCS("E", "Gnomonic", in_sr, geo_lon, geo_lat, "", "", "", ""))
            elif geo_type == "Middle":
                geo_k = 6  # k = 6 for standard extent as per several refs
                geo_sp1, geo_sp2 = conic_sp(geo_lat_max, geo_lat_min, geo_k)
                geo_sp_str = "\n   Standard parallels: SP1 = {0}, SP2 = {1} (K = {2})".format(str(geo_sp1), str(geo_sp2), str(geo_k))
                projections.append(ProjectedCS("E", "Gnomonic", in_sr, geo_lon, geo_lat, "", "", "", ""))
                projections.append(ProjectedCS("I", "Lambert_Conformal_Conic", in_sr, geo_lon, geo_sp1, geo_sp2, 1.0, geo_lat, ""))
                notes.append(conic_note)

# Navigation Routes (Loxodrome): Navigation primarily for maritime use based on sailing a constant bearing between two locations, also referred to as
# Rhumb-Line navigation.
    elif in_purpose == "Navigation Routes (Loxodrome)":
        single_start, geo_lon, geo_lat = nav_condition(in_layer, in_count)
        if in_count == 1:
            geo_type = "1 Line"
            # AddGeomAttributes used as line feature geometry object does not include mid value required for projection definition.
            arcpy.AddGeometryAttributes_management(in_layer, "LENGTH_GEODESIC; LINE_START_MID_END")
            cursor = arcpy.da.SearchCursor(in_layer,["LENGTH_GEO", "START_X", "START_Y", "MID_X", "MID_Y", "END_X", "END_Y"])
            for row in cursor:
                length_geo = row[0]
                lon1, lat1, lon2, lat2 = row[1], row[2], row[5], row[6]
                if lon2 > 180:   # Reset longitude for geodesic line features created using the XY to Line function which cross 180.
                    lon2 = lon2 - 360
                centre_lon, centre_lat = row[3], row[4]
            del row, cursor
            projections.append(ProjectedCS("C", "Mercator", in_sr, centre_lon, 0.0, "", "", "", ""))
            projections.append(ProjectedCS("J", "Loximuthal", in_sr, centre_lon, centre_lat, "", "", "", ""))
            line_details = (lon1, lat1, centre_lon, centre_lat, lon2, lat2, length_geo)
        elif single_start == "TRUE":
            geo_type = "Single Start Point"
            projections.append(ProjectedCS("C", "Mercator", in_sr, geo_lon, 0.0, "", "", "", ""))
            projections.append(ProjectedCS("J", "Loximuthal", in_sr, centre_lon, centre_lat, "", "", "", ""))
        else:
            geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr, in_map, in_df)
            if geo_type == "World":
                projections.append(ProjectedCS("C", "Winkel_Tripel", in_sr, geo_lon, 40.0, "", "", "", ""))
                notes.append(winkel_tripel_note)
                projections.append(ProjectedCS("A", "Robinson", in_sr, geo_lon, "", "", "", "", ""))
                projections.append(ProjectedCS("A", "NaturalEarth", in_sr, geo_lon, "", "", "", "", ""))
                notes.append(world_note)
            elif geo_type == "North Pole" or geo_type == "South Pole":
                projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
                notes.append(sf_note)
                projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
            else:
                projections.append(ProjectedCS("C", "Mercator", in_sr, geo_lon, 0.0, "", "", "", ""))
                projections.append(ProjectedCS("J", "Loximuthal", in_sr, centre_lon, centre_lat, "", "", "", ""))
                geo_k = 6  # k = 6 for standard extent as per several refs
                geo_sp1, geo_sp2 = conic_sp(geo_lat_max, geo_lat_min, geo_k)
                geo_sp_str = "\n   Standard parallels: SP1 = {0}, SP2 = {1} (K = {2})".format(str(geo_sp1), str(geo_sp2), str(geo_k))
                projections.append(ProjectedCS("I", "Lambert_Conformal_Conic", in_sr, geo_lon, geo_sp1, geo_sp2, 1.0, geo_lat, ""))
                projections.append(ProjectedCS("D", "Equidistant_Conic", in_sr, geo_lon, geo_sp1, geo_sp2, geo_lat, "", ""))
                notes.append(conic_note)

# Ranges of Activity: Visualisation of phenomena distance from single/multiple source locations, e.g. radio / telephone signals,
# or tsunami distance from earthquake epicenter.
    elif in_purpose == "Ranges of Activity":
        if in_count == 1:
            geo_type = "1 Feature"
            cursor = arcpy.da.SearchCursor(in_layer, ["SHAPE@XY"])
            for row in cursor:
                geo_lon, geo_lat = row[0]
            del row, cursor
            projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("E", "Orthographic", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
            notes.append(sf_note)
        else:
            geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr, in_map, in_df)
            projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
            notes.append(sf_note)

# Flow Patterns: Maps which generally display a symbolised arrow to highlight the movement of objects between locations, e.g. migration,
# commercial distribution, airline routes.
    elif in_purpose == "Flow Patterns":
        if in_count == 1:
            geo_type = "1 Feature"
            cursor = arcpy.da.SearchCursor(in_layer, ["SHAPE@TRUECENTROID"])
            for row in cursor:
                geo_lon, geo_lat  = row[0]
            del row, cursor
            projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("E", "Orthographic", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
            projections.append(ProjectedCS("L", "Vertical_Near_Side_Perspective", in_sr, geo_lon, geo_lat, 10000000.0, "", "", ""))
            notes.append(sf_note)
            notes.append(vnsp_note)
        elif in_count == 2:
            geo_type = "2 Features"
            array = []
            cursor = arcpy.da.SearchCursor(in_layer, ["SHAPE@TRUECENTROID"])
            for row in cursor:
                point = row[0]
                array.append(point)
            del row, cursor
            lon1, lat1 = array[0]
            lon2, lat2 = array[1]
            projections.append(ProjectedCS("K", "Two_Point_Equidistant", in_sr, lat1, lat2, lon1, lon2, "", ""))
            point_details = (lon1, lat1, lon2, lat2)
        else:
            geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr, in_map, in_df)
            if geo_type == "World":
                projections.append(ProjectedCS("C", "Winkel_Tripel", in_sr, geo_lon, 40.0, "", "", "", ""))
                notes.append(winkel_tripel_note)
                projections.append(ProjectedCS("A", "Robinson", in_sr, geo_lon, "", "", "", "", ""))
                projections.append(ProjectedCS("A", "NaturalEarth", in_sr, geo_lon, "", "", "", "", ""))
                notes.append(world_note)
            else:
                projections.append(ProjectedCS("E", "Orthographic", in_sr, geo_lon, geo_lat, "", "", "", ""))
                if geo_type == "Middle":
                    projections.append(ProjectedCS("L", "Vertical_Near_Side_Perspective", in_sr, geo_lon, geo_lat, 5000000.0, "", "", ""))
                else:
                    projections.append(ProjectedCS("L", "Vertical_Near_Side_Perspective", in_sr, geo_lon, geo_lat, 10000000.0, "", "", ""))
                notes.append(vnsp_note)
                projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
                notes.append(sf_note)
                projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))

# World Index: Maps of the entire globe which overlay information on basic geographical outlines, e.g. time zones, climate zones.
    elif in_purpose == "World Index":
        geo_type = "World"
        geo_lat_max, geo_lat_min = 90.0, -90.0
        projections.append(ProjectedCS("A", "Miller_Cylindrical", in_sr, geo_lon, "", "", "", "", ""))
        projections.append(ProjectedCS("A", "Plate_Carree", in_sr, geo_lon, "", "", "", "", ""))
        projections.append(ProjectedCS("C", "Equidistant_Cylindrical", in_sr, geo_lon, 45.0, "", "", "", ""))
        notes.append(eqd_cyl_note)
        notes.append(world_note)

# Where the footprint size is smaller than the largest settings, add note regarding additional choice of projections used by local or national map grids.
    if geo_type != "World" and geo_type != "Hemisphere":
        notes.append(grid_note)

    # Remove all layers generated by function from data frame.
    layer_list = arcpy.mapping.ListLayers(in_map, "", in_df)
    for layer in layer_list:
        if layer.name == "footprint_mbg":
            arcpy.mapping.RemoveLayer(in_df, layer)

    return projections, notes, geo_type, geo_size_str, geo_lon, geo_lat, geo_lat_max, geo_lat_min, geo_k,\
           geo_sp_str, geo_extent_str, line_details, point_details

# Function to create feature class of fibonacci lattice points required by distortion assessment measures.
def distortion_features(geo_type, in_lyr, in_type, out_gdb, in_sr, map, df):

    # Assign subjective value for natural number N to control the total number of total lattice points created.
    if geo_type == "World":
        N = 100
    elif geo_type == "Hemisphere":
        N = 250
    elif geo_type == "Equatorial" or geo_type == "North Pole" or geo_type == "South Pole":
        N = 500
    else:
        N = 2500

    # Create feature class of fibonaaci lattice spiral lat/longs, see Gonzalez (2010) and Baselga (2018).
    dist_name = "distortion_pts"
    dist_pts = os.path.join(out_gdb, dist_name)
    arcpy.CreateFeatureclass_management(out_gdb, dist_name, "Point", spatial_reference = in_sr)
    P = (2 * N) + 1
    phi = (1 + math.sqrt(5)) / 2   # Calculates the Golden Ratio.
    i = 0 - N
    cursor = arcpy.da.InsertCursor(dist_pts, ["SHAPE@XY"])
    while i <= N:
        lat_rad = math.asin((2 * i) / float(P))
        lat_deg = math.degrees(lat_rad)
        lon_rad = (2 * math.pi * math.fmod(i, phi)) / phi
        lon_deg = math.degrees(lon_rad)
        if lon_deg > 180.0:
            lon_deg = lon_deg - 360.0
        elif lon_deg < -180.0:
            lon_deg = lon_deg + 360.0
        dist_pt = arcpy.Point(lon_deg, lat_deg)
        cursor.insertRow([dist_pt])
        i = i + 1
    del cursor

    # Select features from total distortion feature class within limits based on input layer data type.
    arcpy.MakeFeatureLayer_management(dist_pts, "dist_lyr")
    if in_type == "Point":  # Subjective choice of 1,000Km from point features.
        arcpy.SelectLayerByLocation_management("dist_lyr", "WITHIN_A_DISTANCE_GEODESIC", in_lyr, "1000000 Meters")
    elif in_type == "Polyline":  # Subjective choice of 100Km from line features.
        arcpy.SelectLayerByLocation_management("dist_lyr", "WITHIN_A_DISTANCE_GEODESIC", in_lyr, "100000 Meters")
    else:  # Handles polygons and also raster datasets where an input layer polygon has been generated.
        arcpy.SelectLayerByLocation_management("dist_lyr", "COMPLETELY_WITHIN", in_lyr)
    dist_features = "distortion_features"
    dist_fc = os.path.join(out_gdb, dist_features)
    arcpy.CopyFeatures_management("dist_lyr", dist_fc)
    pt_count = int(arcpy.GetCount_management(dist_fc).getOutput(0))

    # Add field and populate with randomly generated buffer distance. Unit "Meters" specified in field as input spatial ref is GCS so would by default
    # use angular units.
    arcpy.AddField_management(dist_fc, "Buffer_Dist", "TEXT")
    buffer_limit = 1000000  # May be better to consider different values for different sized footprints.
    cursor = arcpy.da.UpdateCursor(dist_fc, ["Buffer_Dist"])
    for row in cursor:
        random_dist = random.randint(0, buffer_limit)
        row[0] = str(random_dist) + " Meters"
        cursor.updateRow(row)
    del row, cursor

    # Remove all layers generated by function from data frame.
    layer_list = arcpy.mapping.ListLayers(map, "", df)
    for layer in layer_list:
        if layer.name == "distortion_pts" or layer.name == "dist_lyr" or layer.name == "distortion_features":
            arcpy.mapping.RemoveLayer(df, layer)

    return dist_fc, pt_count

# Function to calculate distance distortion measures for input points.
def distance_distortion(in_lyr, proj_in_lyr, out_path, in_count, in_sr, map, df):

    # Create two lists of point geometry objects for original and projected input points.
    geo_points = []
    cursor = arcpy.da.SearchCursor(in_lyr, ["SHAPE@XY"])
    for row in cursor:
        geo_x, geo_y = row[0]
        geo_point = arcpy.Point(geo_x, geo_y)
        geo_points.append(geo_point)
    del row, cursor
    proj_points = []
    cursor = arcpy.da.SearchCursor(proj_in_lyr, ["SHAPE@XY"])
    for row in cursor:
        proj_x, proj_y = row[0]
        proj_point = arcpy.Point(proj_x, proj_y)
        proj_points.append(proj_point)
    del row, cursor

    # Create gdb table and add fields for start and end points,
    env.workspace = out_path
    arcpy.CreateTable_management(out_path, "distance_table")
    arcpy.AddField_management("distance_table", "Start_X", "DOUBLE")
    arcpy.AddField_management("distance_table", "Start_Y", "DOUBLE")
    arcpy.AddField_management("distance_table", "End_X", "DOUBLE")
    arcpy.AddField_management("distance_table", "End_Y", "DOUBLE")

    # Randomly generate start and end points for geo distances from input points and add values to gdb table. If statement included to ensure
    # that start and end are different points. Create list of Euclidean distances calculated for equivalent projected points.
    m = 1
    m_limit = in_count * 5
    euclidean_distances = []
    cursor = arcpy.da.InsertCursor("distance_table", ["Start_X", "Start_Y", "End_X", "End_Y"])
    while m <= m_limit:
        start_pt = random.randint(0, (in_count - 1))
        end_pt = random.randint(0, (in_count - 1))
        if end_pt != start_pt:
            start_x, start_y = geo_points[start_pt].X, geo_points[start_pt].Y
            end_x, end_y = geo_points[end_pt].X, geo_points[end_pt].Y
            cursor.insertRow([start_x, start_y, end_x, end_y])
            delta_X = proj_points[end_pt].X - proj_points[start_pt].X
            delta_Y = proj_points[end_pt].Y - proj_points[start_pt].Y
            euclidean_dist = math.sqrt((math.pow(delta_X, 2) + math.pow(delta_Y, 2)))
            euclidean_distances.append(euclidean_dist)
            m = m + 1
    del cursor

    # Create geodesic lines from geo points. Add field to fc and add Euclidean distances from list.
    distance_fc = os.path.join(out_path, "Geodesic_Distances")
    arcpy.XYToLine_management("distance_table", distance_fc, "Start_X", "Start_Y", "End_X", "End_Y", "GEODESIC", spatial_reference = in_sr)
    arcpy.AddGeometryAttributes_management(distance_fc, "LENGTH_GEODESIC", Length_Unit = "METERS")
    arcpy.AddField_management(distance_fc, "Euclidean_Distance", "DOUBLE")
    cursor = arcpy.da.UpdateCursor(distance_fc, ["Euclidean_Distance"])
    a = 0
    for row in cursor:
        row[0] = euclidean_distances[a]
        cursor.updateRow(row)
        a = a + 1
    del row, cursor

    # Use geo and Euclidean distances to calculate distortion measure, using formulas by Peters (1975, cited in Canters, 2002).
    cursor = arcpy.da.SearchCursor(distance_fc, ["LENGTH_GEO", "Euclidean_Distance"])
    sum = 0
    for row in cursor:
        row_index = abs(row[0] - row[1]) / abs(row[0] + row[1])
        sum = sum + row_index
    del row, cursor
    EP = sum / m_limit
    K1 = (1 + EP) / (1 - EP)

    # Remove layer and table generated by function from data frame.
    layer_list = arcpy.mapping.ListLayers(map, "", df)
    for layer in layer_list:
        if layer.name == "Geodesic_Distances":
            arcpy.mapping.RemoveLayer(df, layer)
    table_list = arcpy.mapping.ListTableViews(map, "", df)
    for table in table_list:
        if table.name == "distance_table":
            arcpy.mapping.RemoveTableView(df, table)

    return EP

# Function to calculate area distortion measures for input buffers, using adaptation of formulas by Peters (Canters, 2002).
def area_distortion(in_lyr, in_count):

    cursor = arcpy.da.SearchCursor(in_lyr, ["AREA_GEO", "SHAPE@AREA"])
    sum = 0
    for row in cursor:
        row_index = abs(row[0] - row[1]) / abs(row[0] + row[1])
        sum = sum + row_index
    del row, cursor
    EA = sum / in_count
    KA = (1 + EA) / (1 - EA)

    return EA

# Function to calculate shape distortion measure for input points.
def shape_distortion(geo_fc, circ_fc, in_count, map, df):

    # Make layers to enable OBJECTID selection for geodesic and circular buffers. Determine intersect and union area and calculate shape index by Lee and Sallee (1970).
    arcpy.MakeFeatureLayer_management(geo_fc, "geodesic_lyr")
    arcpy.MakeFeatureLayer_management(circ_fc, "circle_lyr")
    shp_count = 200   # Shape distortion calculations take significant resources so limit required to reduce calculation time.
    if in_count < shp_count:
        shp_count = in_count   # If point selection number is already less than limit then reset maximum to input feature count.
    id = 1
    sum_r = 0
    while id <= shp_count:
        arcpy.SelectLayerByAttribute_management ("geodesic_lyr", "NEW_SELECTION", '"OBJECTID" = ' + str(id))
        arcpy.SelectLayerByAttribute_management ("circle_lyr", "NEW_SELECTION", '"OBJECTID" = ' + str(id))
        intersect_geom = arcpy.Intersect_analysis(["geodesic_lyr", "circle_lyr"], arcpy.Geometry())
        union_geom = arcpy.Union_analysis(["geodesic_lyr", "circle_lyr"], arcpy.Geometry())
        intersect_area = 0
        for i in intersect_geom:
            intersect_area = intersect_area + i.area
        union_area = 0
        for u in union_geom:
            union_area = union_area + abs(u.area)   # On testing some union areas calculated as negative and same effect seen in ArcGIS so absolute values used.
        del i, u
        r = 1 - (intersect_area / union_area)
        sum_r = sum_r + r
        id = id + 1

    # Calculate mean shape index by dividing total of individual indices by number of points.
    shp_index = sum_r / in_count

    # Remove all layers generated by function from data frame.
    layer_list = arcpy.mapping.ListLayers(map, "", df)
    for layer in layer_list:
        if layer.name == "geodesic_lyr" or layer.name == "circle_lyr":
            arcpy.mapping.RemoveLayer(df, layer)

    return shp_index, shp_count

# Function to run distortion assessment measures for each SpatialReference object in candidate list.
def distortion_measures(in_lyr, out_path, sr, in_id, pt_count, in_weight, in_sr, in_map, in_df):

    # Create 4 new feature classes for each spatial reference object: project input points to spatial ref; create geodesic buffers from
    # original points; project geodesic buffers to spatial ref; create Euclidean buffers (circles) using same buffer distance for projected points.
    distortion_proj_fc = os.path.join(out_path, "projected_points_" + str(in_id))
    arcpy.Project_management(in_lyr, distortion_proj_fc, sr)
    geo_buffer_fc = os.path.join(out_path, "geo_buffers_" + str(in_id))
    arcpy.Buffer_analysis(in_lyr, geo_buffer_fc, "Buffer_Dist", method = "GEODESIC")
    arcpy.AddGeometryAttributes_management(geo_buffer_fc, "AREA_GEODESIC", Area_Unit = "SQUARE_METERS")
    proj_buffer_fc = os.path.join(out_path, "projected_buffers_" + str(in_id))
    arcpy.Project_management(geo_buffer_fc, proj_buffer_fc, sr, preserve_shape = "PRESERVE_SHAPE")
    circle_fc = os.path.join(out_path, "buffer_projected_points_" + str(in_id))
    arcpy.Buffer_analysis(distortion_proj_fc, circle_fc, "Buffer_Dist", method = "PLANAR")

    # Run distortion measure functions and write results to list of tuples.
    if in_weight[0] == 1:
        out_EP = distance_distortion(in_lyr, distortion_proj_fc, out_path, pt_count, in_sr, in_map, in_df)
    else:
        out_EP = 0.0
    if in_weight[1] == 1:
        out_EA = area_distortion(proj_buffer_fc, pt_count)
    else:
        out_EA = 0.0
    if in_weight[2] == 1:
        shape_index, shape_count = shape_distortion(proj_buffer_fc, circle_fc, pt_count, in_map, in_df)
    else:
        shape_index, shape_count = 0.0, 0
    combined_index = out_EP + out_EA + shape_index
    dist_measures = (out_EP, out_EA, shape_index, shape_count, combined_index)

    # Remove all layers generated by function from data frame.
    layer_list = arcpy.mapping.ListLayers(in_map, "projected_points_*", in_df)
    for layer in layer_list:
        arcpy.mapping.RemoveLayer(in_df, layer)
    del layer
    layer_list = arcpy.mapping.ListLayers(in_map, "geo_buffers_*", in_df)
    for layer in layer_list:
        arcpy.mapping.RemoveLayer(in_df, layer)
    del layer
    layer_list = arcpy.mapping.ListLayers(in_map, "projected_buffers_*", in_df)
    for layer in layer_list:
        arcpy.mapping.RemoveLayer(in_df, layer)
    del layer
    layer_list = arcpy.mapping.ListLayers(in_map, "buffer_projected_points_*", in_df)
    for layer in layer_list:
        arcpy.mapping.RemoveLayer(in_df, layer)
    del layer

    return dist_measures

# Function to determine optimal conic projection based on the value of k which minimises the Combined Index of distortion.
def optimal_conic(in_sr, in_proj, id, pcs_name, in_k, in_lon, in_lat, in_lat_max, in_lat_min, dist_fc, out_path, gdb_path, pt_count, in_weights, map, df):

# For each potential k value (except the original setting) create list of candidate ProjectedCS objects.
    k_list = [3, 4, 5, 6, 7]
    k_index = k_list.index(in_k)
    k_list.pop(k_index)
    conic_list = []
    for k in k_list:
        in_sp1, in_sp2 = conic_sp(in_lat_max, in_lat_min, k)
        sp_str = "\n   Standard Parallels = SP1: {0}, SP2: {1} (K = {2})".format(str(in_sp1), str(in_sp2), str(k))
        if pcs_name == "Custom_Lambert_Conformal_Conic":
            pcs = ProjectedCS("I", "Lambert_Conformal_Conic", in_sr, in_lon, in_sp1, in_sp2, 1.0, in_lat, "")
        elif pcs_name == "Custom_Equidistant_Conic":
            pcs = ProjectedCS("D", "Equidistant_Conic", in_sr, in_lon, in_sp1, in_sp2, in_lat, "", "")
        elif pcs_name == "Custom_Albers":
            pcs = ProjectedCS("D", "Albers", in_sr, in_lon, in_sp1, in_sp2, in_lat, "", "")

# Create spatial reference object, .prj file and list of details.
        pcs_sr, pcs_string = pcs.create_sr()
        pcs_prj = "{0}_PCS_{1}.prj".format(in_proj, str(id))
        pcs_prj_filename = os.path.join(out_path, pcs_prj)
        pcs_prj_file = open(pcs_prj_filename, "w")
        pcs_prj_file.write(pcs_string)
        pcs_prj_file.close()

# Run distortion assessment measures for each candidate, compare and determine
        dist = distortion_measures(dist_fc, gdb_path, pcs_sr, id, pt_count, in_weights, in_sr, map, df)
        output_details = (pcs_sr, id, pcs_string, pcs_prj, dist[0], dist[1], dist[2], dist[3], dist[4], sp_str)
        conic_list.append(output_details)
        id = id + 1
    del k

    return conic_list

# Function to determine optimal scale factor which minimises Combined Index of distortion for Transverse Mercator and Stereographic projections.
def optimal_sf(in_sr, in_proj, id, pcs_name, in_lon, in_lat, dist_fc, out_path, gdb_path, pt_count, in_weights, map, df):

# For each potential scale factor value (except the original setting of 1.0) create list of candidate ProjectedCS objects.
    sf_list = []
    sf = 0.9996
    while sf <= 0.9999:
        if pcs_name == "Custom_Transverse_Mercator":
            pcs = ProjectedCS("F", "Transverse_Mercator", in_sr, in_lon, sf, in_lat, "", "", "")
        elif pcs_name == "Custom_Stereographic":
            pcs = ProjectedCS("F", "Stereographic", in_sr, in_lon, sf, in_lat, "", "", "")

# Create spatial reference object, .prj file and list of details.
        pcs_sr, pcs_string = pcs.create_sr()
        pcs_prj = "{0}_PCS_{1}.prj".format(in_proj, str(id))
        pcs_prj_filename = os.path.join(out_path, pcs_prj)
        pcs_prj_file = open(pcs_prj_filename, "w")
        pcs_prj_file.write(pcs_string)
        pcs_prj_file.close()

# Run distortion assessment measures for each candidate and create list of pcs and distortion details.
        dist = distortion_measures(dist_fc, gdb_path, pcs_sr, id, pt_count, in_weights, in_sr, map, df)
        output_details = (pcs_sr, id, pcs_string, pcs_prj, dist[0], dist[1], dist[2], dist[3], dist[4], "")
        sf_list.append(output_details)
        sf = sf + 0.0001
        id = id + 1
    del k

    return sf_list

# Class to handle user input of a Project Name using free text (max 10 characters).
class Project_Name(object):
    """Implementation for AutomatedProjectionSelection_addin.project_name (ComboBox)"""
    def __init__(self):
        self.editable = True
        self.enabled = True
        self.width = 'WWWWWWWW'
        self.value = ""

    def onEditChange(self, text):
    # Ensure current Map Document is a saved MXD file. Set this_map, this_df and output_path variables for later use.
        try:
            self.this_map = arcpy.mapping.MapDocument("CURRENT")
            self.this_df = self.this_map.activeDataFrame
            if self.this_map.filePath != "":
                desc_map = arcpy.Describe(self.this_map.filePath)
                self.output_path = desc_map.path

    # Set project name variable from user text and enable other user inputs.
                self.value = text
                if len(self.value) > 10:
                    pythonaddins.MessageBox("Project name has maximum of 10 characters!", "Warning", 0)
                    self.value = self.value[:10]
                    self.refresh()
                self.name = self.value
                input_purpose.enabled = True
                reset_inputs.enabled = True

    # Warn user if map document needs to be defined.
            else:
                pythonaddins.MessageBox("Current map document undefined. Please open MXD file or save blank document to enable tools to run.",\
                                        "Warning", 0)
                self.value = ""
                self.refresh()

    # Use Try/Except statement to catch any exceptions with selecting project name.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem selecting project name. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)
            self.value = ""
            self.refresh()

    def refresh(self):
        pass

# Class to handle user input using the Purpose drop-down list.
class Input_Purpose(object):
    """Implementation for AutomatedProjectionSelection_addin.input_purpose (ComboBox)"""
    def __init__(self):
        self.items = ["General Reference", "Thematic Vector", "Thematic Raster", "Geospatial Analysis (Distance)", "Geospatial Analysis (Area)",\
                        "Navigation Routes (Geodesic)", "Navigation Routes (Loxodrome)", "Ranges of Activity", "Flow Patterns", "World Index"]
        self.editable = True
        self.enabled = False
        self.dropdownWidth = 'WWWWWWWWWWWWWWWW'
        self.width = 'WWWWWWWWWWWWWWWW'
        self.value = ""

    def onSelChange(self, selection):

    # If purpose changed then reset various other inputs, including input layer selection as each purpose has specific data type requirements.
        try:
            run_tool.checked = False
            if input_layer.value != "":
                input_layer.value = ""
                input_layer.refresh()
            if spatial_reference.value != "":
                spatial_reference.value = ""
                spatial_reference.refresh()
            spatial_reference.enabled = False
            apply_projection.enabled = False

    # Set purpose variable from user selection and ensure footprint definition user inputs are enabled.
            self.purpose = selection
            input_layer.enabled = True

    # Create list of appropriate data types and set distortion weights (distance / area / shape) for each purpose.
            if self.purpose == "General Reference":
                self.data_types = ["Polygon"]
                self.dist_weights = (0, 1, 1)
            elif self.purpose == "Thematic Vector":
                self.data_types = ["Polygon"]
                self.dist_weights = (0, 1, 0)
            elif self.purpose == "Thematic Raster":
                self.data_types = ["RasterDataset", "RasterCatalog"]
                self.dist_weights = (0, 1, 0)
            elif self.purpose == "Geospatial Analysis (Distance)":
                self.data_types = ["Point", "Polygon", "Polyline", "RasterDataset"]
                self.dist_weights = (1, 0, 0)
            elif self.purpose == "Geospatial Analysis (Area)":
                self.data_types = ["Polygon", "RasterDataset"]
                self.dist_weights = (0, 1, 0)
            elif self.purpose == "Navigation Routes (Geodesic)":
                self.data_types = ["Polyline"]
                self.dist_weights = (0, 1, 1)
            elif self.purpose == "Navigation Routes (Loxodrome)":
                self.data_types = ["Polyline"]
                self.dist_weights = (0, 1, 1)
            elif self.purpose == "Ranges of Activity":
                self.data_types = ["Point"]
                self.dist_weights = (1, 0, 0)
            elif self.purpose == "Flow Patterns":
                self.data_types = ["Point", "Polygon"]
                self.dist_weights = (0, 1, 1)
            elif self.purpose == "World Index":
                self.data_types = ["Polygon"]
                self.dist_weights = (0, 1, 1)

    # Use Try/Except statement to catch any exceptions with selecting input purpose.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem selecting input purpose. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)

    def refresh(self):
        pass

# Class to handle user input using the Input Layer drop-down list.
class Input_Layer(object):
    """Implementation for AutomatedProjectionSelection_addin.input_layer (ComboBox)"""
    def __init__(self):
        self.editable = True
        self.enabled = False
        self.dropdownWidth = 'WWWWWWWWWWWW'
        self.width = 'WWWWWWWWWWWW'
        self.value = ""

    def onSelChange(self, selection):

    # If input layer changed then reset various other inputs.
        try:
            run_tool.checked = False
            if spatial_reference.value != "":
                spatial_reference.value = ""
                spatial_reference.refresh()
            spatial_reference.enabled = False
            apply_projection.enabled = False

    # Set input layer variable from user selection and zoom to extent of selected layer. Enable Run Tool button.
            self.select_layer = selection
            project_name.this_df.extent = self.select_layer.getExtent()
            run_tool.enabled = True

    # Use Try/Except statement to catch any exceptions with selecting input layer.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem selecting input layer. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)

    def onFocus(self, focused):

    # On selecting drop-down generate list of layers in current data frame of appropriate data type and geographic coordinate system.
        try:
            if focused:
                self.items = []
                layer_list = arcpy.mapping.ListLayers(project_name.this_map, "", project_name.this_df)
                for layer in layer_list:
                    desc_layer = arcpy.Describe(layer)
                    if desc_layer.dataType != "GroupLayer":
                        layer_sr = desc_layer.spatialReference
                        if layer_sr.type == "Geographic":
                            if desc_layer.dataType == "FeatureLayer":
                                desc_featurelayer = desc_layer.featureClass
                                data_type = desc_featurelayer.shapeType
                            elif desc_layer.dataType == "RasterLayer" or desc_layer.dataType == "RasterCatalogLayer":
                                desc_rasterlayer = desc_layer.dataElement
                                data_type = desc_rasterlayer.datasetType
                            else:
                                data_type = "Other"
                            if data_type in input_purpose.data_types:
                                self.items.append(layer)

    # Use Try/Except statement to catch any exceptions with selecting input layer.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem selecting input layer. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)

    def refresh(self):
        pass

# Class to handle running of selected tool.
class Run_Tool(object):
    """Implementation for AutomatedProjectionSelection_addin.run_tool (Button)"""
    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):

    # Set variables from input fields, then disable all input fields to ensure Reset clicked before running again..
        try:
            proj_name = project_name.name
            purpose = input_purpose.purpose
            selected_layer = input_layer.select_layer
            layer_name = selected_layer.name
            project_name.enabled = False
            input_purpose.enabled = False
            input_layer.enabled = False

    # Set name for results text file and create File Geodatabase, unless Project Name has already been used then delete all previous feature classes in GDB.
    # Delete any .prj files with the Project Name in output folder.
            self.result_name = project_name.name + "_results.txt"
            gdb_name = project_name.name + "_temp.gdb"
            self.gdb_path = os.path.join(project_name.output_path, gdb_name)
            env.workspace = self.gdb_path
            if arcpy.Exists(self.gdb_path) == False:
                arcpy.CreateFileGDB_management(project_name.output_path, gdb_name)
            else:
                fc_list = arcpy.ListFeatureClasses()
                for fc in fc_list:
                    arcpy.Delete_management(fc)
            env.workspace = project_name.output_path
            prj_name = project_name.name + "*.prj"
            prj_list = arcpy.ListFiles(prj_name)
            for prj in prj_list:
                arcpy.Delete_management(prj)

    # Message for user to verify selections. Continue if satisfied, exit if not.
            tool_message = pythonaddins.MessageBox("Output File Path: {0}\nResults Text File: {1}\nGDB Name: {2}\n\nThe output path will contain a results text file, "\
                                                   "a GDB containing temporary files and candidate projection .prj files.\n\nIf a tool has already been run using the "\
                                                   "same project name then previous GDB feature classes and .prj files will be deleted. Please rename now if you wish "\
                                                   "to retain these previous items.\n\nPlease be aware that the following process can take many minutes to finish."
                                                   .format(project_name.output_path, self.result_name, gdb_name), "Message", 1)

    # Use Try/Except statement to catch any exceptions with setting the input variables.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem setting the input variables. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)

    # For selected input layer, either make copy of feature class in GDB or use raster_polygon function to create extent layer. Set spatial reference, data type
    # and feature count variables required by selection tools. Check for zero features if layer is a feature class.
        try:
            if tool_message == "OK":
                desc_layer = arcpy.Describe(selected_layer)
                input_sr = desc_layer.spatialReference
                if desc_layer.dataType == "FeatureLayer":
                    footprint_layer = os.path.join(self.gdb_path, "input_layer")
                    arcpy.CopyFeatures_management(selected_layer, footprint_layer)
                    desc_featurelayer = desc_layer.featureClass
                    data_type = desc_featurelayer.shapeType
                    feature_count = int(arcpy.GetCount_management(footprint_layer).getOutput(0))
                    if feature_count == 0:
                        pythonaddins.MessageBox("Input has no features. Please try again...", "Warning", 0)
                        sys.exit()
                    elif feature_count > 1:
                        if purpose == "General Reference" or purpose == "Thematic Vector":
                            pythonaddins.MessageBox("For {0} purpose, the geographic footprint characteristics can only be determined adequately from a single "\
                                                    "polygon which completely encompasses the desired footprint of interest. Please select an input dataset which "
                                                    "has a single polygon feature.".format(purpose), "Warning", 0)
                            sys.exit()
                elif desc_layer.dataType == "RasterLayer" or desc_layer.dataType == "RasterCatalogLayer":
                    footprint_layer = raster_polygon(selected_layer, self.gdb_path, input_sr, project_name.this_map, project_name.this_df)
                    desc_rasterlayer = desc_layer.dataElement
                    data_type = desc_rasterlayer.datasetType
                    feature_count = 1

    # Run projection selection function to determine geographic footprint characteristics and candidate projection options to enable creation of
    # spatial reference objects and .prj files for coordinate system selection drop-down list.
                candidates, proj_notes, input_type, output_size_str, input_lon, input_lat, input_lat_max, input_lat_min, input_k,\
                output_sp_str, output_extent_str, input_line, input_2pt = projection_selection(self.gdb_path, purpose, footprint_layer, data_type, feature_count,\
                                                                                                input_sr, project_name.this_map, project_name.this_df)
                pcs_list = []
                pcs_id = 1
                for candidate in candidates:
                    pcs_sr, pcs_string = candidate.create_sr()
                    pcs_prj = "{0}_PCS_{1}.prj".format(proj_name, str(pcs_id))
                    pcs_prj_filename = os.path.join(project_name.output_path, pcs_prj)
                    pcs_prj_file = open(pcs_prj_filename, "w")
                    pcs_prj_file.write(pcs_string)
                    pcs_prj_file.close()
                    pcs_details = (pcs_sr, pcs_id, pcs_string, pcs_prj)
                    pcs_list.append(pcs_details)
                    pcs_id = pcs_id + 1
                del candidate

    # Create list of coordinate systems for distortion assessment using the selection tool output. Run distortion features function to generate
    # dataset of fibonacci lattice points then run distortion assessment measures function for each item in the list. Create list of output details
    # needed for results file and coordinate system selection drop-down list.
                self.output_list = []
                distortion_fc, point_count = distortion_features(input_type, footprint_layer, data_type, self.gdb_path, input_sr, project_name.this_map, project_name.this_df)
                for prj in pcs_list:
                    dist = distortion_measures(distortion_fc, self.gdb_path, prj[0], prj[1], point_count, input_purpose.dist_weights, input_sr,\
                                                project_name.this_map, project_name.this_df)
                    if prj[0].name == "Custom_Lambert_Conformal_Conic" or prj[0].name == "Custom_Equidistant_Conic" or prj[0].name == "Custom_Albers":
                        output_details = (prj[0], prj[1], prj[2], prj[3], dist[0], dist[1], dist[2], dist[3], dist[4], output_sp_str)
                    else:
                        output_details = (prj[0], prj[1], prj[2], prj[3], dist[0], dist[1], dist[2], dist[3], dist[4], "")
                    self.output_list.append(output_details)
                del prj

    # Sort the list of candidates into a ranked order based on the Combined Index distortion measure.
                def getKey(item):
                    return item[8]
                self.output_list.sort(key = getKey)

    # For certain projection selections which are ranked first (to avoid overwhelming the user with candidates), optimise parameter values which minimise
    # Combined Index by carrying out additional distortion assessment for a small number of variations of the key projection parameter. Create additional
    # list of pcs and distortion measure details and add to output_list.
                rank1 = self.output_list[0]
                sr_name = rank1[0].name
                if sr_name == "Custom_Lambert_Conformal_Conic" or sr_name == "Custom_Equidistant_Conic" or sr_name == "Custom_Albers":
                    optimal_list = optimal_conic(input_sr, proj_name, pcs_id, sr_name, input_k, input_lon, input_lat, input_lat_max, input_lat_min,\
                                                    distortion_fc, project_name.output_path, self.gdb_path, point_count, input_purpose.dist_weights,\
                                                    project_name.this_map, project_name.this_df)
                    for details in optimal_list:
                        self.output_list.append(details)
                if sr_name == "Custom_Transverse_Mercator" or sr_name == "Custom_Stereographic":
                    optimal_list = optimal_sf(input_sr, proj_name, pcs_id, sr_name, input_lon, input_lat, distortion_fc, project_name.output_path,\
                                                self.gdb_path, point_count, input_purpose.dist_weights, project_name.this_map, project_name.this_df)
                    for details in optimal_list:
                        self.output_list.append(details)

    # Sort the output list again to get a new ranked order including the optimisation candidates.
                def getKey(item):
                    return item[8]
                self.output_list.sort(key = getKey)

    # Check button to show tool has been run, and enable user selections of spatial reference object.
                self.checked = True
                spatial_reference.enabled = True

    # Remove input layer copy from data frame.
                layer_list = arcpy.mapping.ListLayers(project_name.this_map, "", project_name.this_df)
                for layer in layer_list:
                    if layer.name == "input_layer" or layer.name == "raster_extent":
                        arcpy.mapping.RemoveLayer(project_name.this_df, layer)

    # Set up results text file and write standard input information.
                result_filename = os.path.join(project_name.output_path, self.result_name)
                result_file = open(result_filename, "w")
                result_file.write("Selected Purpose = {0}\nProject name = {1}\nInput feature class = {2}\nInput spatial reference = {3}\n"\
                                  "Output file path = {4}\nGDB name = {5}\n\n"\
                                  .format(purpose, proj_name, layer_name, input_sr.name, project_name.output_path, gdb_name))

    # Write footprint characteristics to results text file.
                result_file.write("Input footprint type = {0}\n".format(input_type))
                if input_type == "1 Feature":
                    result_file.write("Input footprint feature lon/lat = {0:.6f}, {1:.6f}\n".format(input_lon, input_lat))
                elif input_type == "2 Features":
                    result_file.write("Input footprint feature 1 lon/lat = {0:.6f}, {1:.6f}\nInput footprint feature 2 lon/lat = {2:.6f}, {3:.6f}\n"\
                                      .format(input_2pt[0], input_2pt[1], input_2pt[2], input_2pt[3]))
                elif input_type == "1 Line":
                    result_file.write("Input line start lon/lat = {0:.6f}, {1:.6f}\nInput line centre lon/lat = {2:.6f}, {3:.6f}\n"\
                                      "Input line end lon/lat = {4:.6f}, {5:.6f}\nInput line geodesic length (metres) = {6:.6f}\n"\
                                      .format(input_line[0], input_line[1], input_line[2], input_line[3], input_line[4], input_line[5], input_line[6]))
                elif input_type == "Single Start Point":
                    result_file.write("Input lines start lon/lat = {0:.6f}, {1:.6f}\n".format(input_lon, input_lat))
                else:
                    result_file.write("Input footprint centre lon/lat = {0:.6f}, {1:.6f}\nInput footprint latitude max/min = {2:.6f}, {3:.6f}\n"\
                                      .format(input_lon, input_lat, input_lat_max, input_lat_min))
                result_file.write("Input footprint size (km squared) = {0}\nInput footprint extent group = {1}\nInput footprint standard parallels for conic projections = {2}\n"\
                                  .format(output_size_str, output_extent_str, output_sp_str))

    # For distortion assessment, write list of candidates and distortion measure calculations to results text file. For purposes using area distortion only
    # use 10 decimal places to distinguish results, else use 5 decimal places.
                result_file.write("\nNumber of points selected for distortion assessment = {0}\n".format(str(point_count)))
                result_file.write("Distortion assessment weights: Distance = {0}; Area = {1}; Shape = {2}\n\n"\
                                    .format(str(input_purpose.dist_weights[0]), str(input_purpose.dist_weights[1]), str(input_purpose.dist_weights[2])))
                result_file.write("Distortion assessment results (listed in order of least Combined Index):\n\n")
                rank = 1
                for output in self.output_list:
                    result_file.write("Candidate {0} = {1} ({2}){3}\n   Distortion Measures:\n".format(rank, output[3], output[0].name, output[9]))
                    if input_purpose.dist_weights == (0, 1, 0):
                        result_file.write("     Distance:   Not measured\n")
                        result_file.write("     Area:       EA = {0:.10f}\n".format(output[5]))
                        result_file.write("     Shape:      Not measured\n")
                        result_file.write("     Combined Index = {0:.10f}\n\n".format(output[8]))
                    else:
                        if input_purpose.dist_weights[0] == 0:
                            result_file.write("     Distance:   Not measured\n")
                        else:
                            result_file.write("     Distance:   EP = {0:.5f}\n".format(output[4]))
                        if input_purpose.dist_weights[1] == 0:
                            result_file.write("     Area:       Not measured\n")
                        else:
                            result_file.write("     Area:       EA = {0:.5f}\n".format(output[5]))
                        if input_purpose.dist_weights[2] == 0:
                            result_file.write("     Shape:      Not measured\n")
                        else:
                            result_file.write("     Shape:   Index = {0:.5f}; No of points selected for shape calculation = {1}\n".format(output[6], str(output[7])))
                        result_file.write("     Combined Index = {0:.5f}\n\n".format(output[8]))
                    rank = rank + 1
                del output

    # Write projection notes to results text file.
                if proj_notes != []:
                    result_file.write("Projection notes:\n")
                    for note in proj_notes:
                        result_file.write(note)
                    del note

    # Close the results text file and provide a message to confirm output .prj files created and advise user on their next steps.
                result_file.close()
                pythonaddins.MessageBox("Completed. Please review results text file for full details.\n\nSelect a .prj file name from the Projected Coord System "\
                                        "drop-down list and click Apply to modify the data frame coordinate system. Please be aware that this function uses the "\
                                        "ArcGIS project-on-the-fly capability, so layers have retained their coordinate system settings but can be re-projected "\
                                        "as necessary using the Project tool.\n\nEach candidate can be reviewed before making a final selection, ideally using "\
                                        "Layout View with an appropriate graticule to best view the impact of the selection on the geographical footprint of interest.",
                                        "Message", 0)

    # Use Try/Except statement to catch any exceptions with running the selection and distortion assessment tools.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem running the tool. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)

# Class to handle user selection of .prj file name from ranked list of candidates generated by running the selection tool.
class Spatial_Reference(object):
    """Implementation for AutomatedProjectionSelection_addin.spatial_reference (ComboBox)"""
    def __init__(self):
        self.editable = True
        self.enabled = False
        self.dropdownWidth = 'WWWWWWWWWWWW'
        self.width = 'WWWWWWWWWWWW'
        self.value = ""

    def onSelChange(self, selection):

    # Set spatial reference object from the same item in the output_list to that chosen by the user using the .prj file name.
        try:
            prj_file = selection
            for details in run_tool.output_list:
                if details[3] == prj_file:
                    self.spatial_ref = details[0]
                    self.wkt = details[2]
            apply_projection.enabled = True

    # Use Try/Except statement to catch any exceptions with selecting spatial reference.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem selecting spatial reference. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)

    def onFocus(self, focused):

    # Set drop-down list to be the ranked list of .prj file names previously generated by selection and distortion assessment tools.
        try:
            if focused:
                self.items = []
                for output in run_tool.output_list:
                    self.items.append(output[3])

    # Use Try/Except statement to catch any exceptions with selecting spatial reference.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem selecting spatial reference. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)

    def refresh(self):
        pass

# Class to handle modification of data frame coordinate system setting.
class Apply_Projection(object):
    """Implementation for AutomatedProjectionSelection_addin.apply_projection (Button)"""
    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):

    # Set spatial reference object based on user selection of coordinate system from drop-down list.
        try:
            if spatial_reference.value != "":
                apply_sr = spatial_reference.spatial_ref
                apply_wkt = spatial_reference.wkt
            else:
                pythonaddins.MessageBox("Something has gone wrong with the selection of spatial reference. Please try again.", "Warning", 0)
                apply_sr = None

    # Modify data frame coordinate system setting and zoom to extent of selected dataset.
            if apply_sr != None:
                select_layer = input_layer.select_layer
                project_name.this_df.spatialReference = apply_sr
                project_name.this_df.extent = select_layer.getExtent()

    # Use Try/Except statement to catch any exceptions with applying the selected spatial reference.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem applying the selected spatial reference. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)

# Class to handle button click which resets toolbar to enable user to start again.
class Reset_Inputs(object):
    """Implementation for AutomatedProjectionSelection_addin.reset_inputs (Button)"""
    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        try:
            run_tool.checked = False
            run_tool.enabled = False
            project_name.name = ""
            project_name.value = ""
            project_name.refresh()
            project_name.enabled = True
            if input_purpose.value != "":
                input_purpose.purpose = ""
                input_purpose.value = ""
                input_purpose.refresh()
            input_purpose.enabled = False
            if input_layer.value != "":
                input_layer.value = ""
                input_layer.refresh()
            input_layer.enabled = False
            if spatial_reference.value != "":
                spatial_reference.value = ""
                spatial_reference.refresh()
            spatial_reference.enabled = False
            apply_projection.enabled = False
            self.enabled = False

    # Use Try/Except statement to catch any exceptions with resetting inputs.
        except Exception as e:
            pythonaddins.MessageBox("There has been a problem resetting inputs. Python error message:\n\n{0}".format(e),\
                                    "Error", 0)
