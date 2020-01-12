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
#-----------------------------------------------------------------------------------------------------------------------------
# Title: Automated Projection Selection Tools(ArcGIS Desktop Python Toolbox)
#
# Purpose: Toolbox version of functions incorporated into the Automated Projection Selection for ArcGIS (ArcMap Python Add-in).
#          Provides the user with the opportunity to run the following tools: Projection Selection only for an area of interest
#          and purpose, i.e. without assessing the distortion of the candidate projections; Distortion Assessment for a chosen
#          .prj file and the user's area of interest and purpose; Geographic Centre determination using the iterative process
#          based on the azimuthal equidistant projection; Fibonacci Lattice generation of equally spaced assessment points.
#
# Author: Paul Gosling
# Version: 1.0
# Published: 12 Jan 2020
#------------------------------------------------------------------------------------------------------------------------------

import arcpy
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
def footprint_type(in_lyr, out_path, in_sr):

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
    equidistant_fc = os.path.join(out_path, "area_centroid")
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

    # If footprint geometry encompasses single Pole or centroid in Polar reion then set appropriate type, reset latitude to Pole and longitude to Greenwich Meridian
    # if footprint encompasses the entire Polar region.
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

    # If area geometry crosses Equator or centroid is within certain latitude then set appropriate type and reset latitude centre to Equator.
    elif (lat_max > 0 and lat_min < 0) or abs(lat_centre) < 15:
        in_type = "Equatorial"
        lat_centre = 0.0

    # By process of elimination all remaining areas are at Middle latitudes so set appropriate type.
    else:
        in_type = "Middle"

    return in_type, in_size, lon_centre, lat_centre, lat_max, lat_min

# Function to define area extent characteristics grouping.
def footprint_extent(in_lyr, out_path, in_sr, in_lon, in_lat):

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

    return extent_group, extent_ratio

# Function to determine input type grouping where function uses set of individual point/line/polygon features to define footprint of interest.
# Returns footprint characteristics to be used in PCS definitions.
def location_type(in_purpose, in_lyr, out_path, in_sr):

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
    # if footprint encompasses the entire Polar region.
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
def raster_polygon(in_lyr, out_path, in_sr):

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
def projection_selection(out_gdb, in_purpose, in_layer, in_type, in_count, in_sr):

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
        geo_type, geo_size, geo_lon, geo_lat, geo_lat_max, geo_lat_min = footprint_type(in_layer, out_gdb, in_sr)
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
            geo_extent_group, geo_extent_ratio = footprint_extent(in_layer, out_gdb, in_sr, geo_lon, geo_lat)
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
        geo_type, geo_size, geo_lon, geo_lat, geo_lat_max, geo_lat_min = footprint_type(in_layer, out_gdb, in_sr)
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
            geo_extent_group, geo_extent_ratio = footprint_extent(in_layer, out_gdb, in_sr, geo_lon, geo_lat)
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
        geo_type, geo_size, geo_lon, geo_lat, geo_lat_max, geo_lat_min = footprint_type(in_layer, out_gdb, in_sr)
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
        geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr)
        if geo_type == "World":
            pythonaddins.MessageBox("ArcGIS analysis techniques are not applicable for an area of this size. Please select a different input dataset.",\
                                    "Warning", 0)
            sys.exit()
        elif geo_type == "North Pole" or geo_type == "South Pole":
            projections.append(ProjectedCS("B", "Azimuthal_Equidistant", in_sr, geo_lon, geo_lat, "", "", "", ""))
            projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
            notes.append(sf_note)
        else:
            geo_extent_group, geo_extent_ratio = footprint_extent(in_layer, out_gdb, in_sr, geo_lon, geo_lat)
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
            geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr)
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
            geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr)
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
            geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr)
            projections.append(ProjectedCS("F", "Stereographic", in_sr, geo_lon, 1.0, geo_lat, "", "", ""))
            notes.append(sf_note)

# Flow Patterns: Maps which generally display a symbolised arrow to highlight the movement of objects between locations, e.g. migration,
# commercial distribution, airline routes.
    elif in_purpose == "Flow Patterns":
        if in_count == 1:
            geo_type = "1 Feature"
            cursor = arcpy.da.SearchCursor(in_layer, ["SHAPE@TRUECENTROID"])
            for row in cursor:
                geo_lon, geo_lat = row[0]
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
            geo_type, geo_lon, geo_lat, geo_lat_max, geo_lat_min = location_type(in_purpose, in_layer, out_gdb, in_sr)
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

    return projections, notes, geo_type, geo_size_str, geo_lon, geo_lat, geo_lat_max, geo_lat_min, geo_k,\
           geo_sp_str, geo_extent_str, line_details, point_details

# Function to create feature class of fibonacci lattice points required by distortion assessment measures.
def distortion_features(geo_type, in_lyr, in_type, out_gdb, in_sr):

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

    return dist_fc, pt_count

# Function to calculate distance distortion measures for input points.
def distance_distortion(in_lyr, proj_in_lyr, out_path, in_count, in_sr):

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
def shape_distortion(geo_fc, circ_fc, in_count):

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

    return shp_index, shp_count

# Function to run distortion assessment measures for each SpatialReference object in candidate list.
def distortion_measures(in_lyr, out_path, sr, pt_count, in_weight, in_sr):

    # Create 4 new feature classes for each spatial reference object: project input points to spatial ref; create geodesic buffers from
    # original points; project geodesic buffers to spatial ref; create Euclidean buffers (circles) using same buffer distance for projected points.
    distortion_proj_fc = os.path.join(out_path, "projected_points")
    arcpy.Project_management(in_lyr, distortion_proj_fc, sr)
    geo_buffer_fc = os.path.join(out_path, "geo_buffers")
    arcpy.Buffer_analysis(in_lyr, geo_buffer_fc, "Buffer_Dist", method = "GEODESIC")
    arcpy.AddGeometryAttributes_management(geo_buffer_fc, "AREA_GEODESIC", Area_Unit = "SQUARE_METERS")
    proj_buffer_fc = os.path.join(out_path, "projected_buffers")
    arcpy.Project_management(geo_buffer_fc, proj_buffer_fc, sr, preserve_shape = "PRESERVE_SHAPE")
    circle_fc = os.path.join(out_path, "buffer_projected_points")
    arcpy.Buffer_analysis(distortion_proj_fc, circle_fc, "Buffer_Dist", method = "PLANAR")

    # Run distortion measure functions and write results to list of tuples.
    if in_weight[0] == 1:
        out_EP = distance_distortion(in_lyr, distortion_proj_fc, out_path, pt_count, in_sr)
    else:
        out_EP = 0.0
    if in_weight[1] == 1:
        out_EA = area_distortion(proj_buffer_fc, pt_count)
    else:
        out_EA = 0.0
    if in_weight[2] == 1:
        shape_index, shape_count = shape_distortion(proj_buffer_fc, circle_fc, pt_count)
    else:
        shape_index, shape_count = 0.0, 0
    combined_index = out_EP + out_EA + shape_index
    dist_measures = (out_EP, out_EA, shape_index, shape_count, combined_index)

    return dist_measures

# Function to determine optimal conic projection based on the value of k which minimises the Combined Index of distortion.
def optimal_conic(in_sr, in_proj, id, pcs_name, in_k, in_lon, in_lat, in_lat_max, in_lat_min, dist_fc, out_path, gdb_path, pt_count, in_weights):

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
        dist = distortion_measures(dist_fc, gdb_path, pcs_sr, id, pt_count, in_weights, in_sr)
        output_details = (pcs_sr, id, pcs_string, pcs_prj, dist[0], dist[1], dist[2], dist[3], dist[4], sp_str)
        conic_list.append(output_details)
        id = id + 1
    del k

    return conic_list

# Function to determine optimal scale factor which minimises Combined Index of distortion for Transverse Mercator and Stereographic projections.
def optimal_sf(in_sr, in_proj, id, pcs_name, in_lon, in_lat, dist_fc, out_path, gdb_path, pt_count, in_weights):

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
        dist = distortion_measures(dist_fc, gdb_path, pcs_sr, id, pt_count, in_weights, in_sr)
        output_details = (pcs_sr, id, pcs_string, pcs_prj, dist[0], dist[1], dist[2], dist[3], dist[4], "")
        sf_list.append(output_details)
        sf = sf + 0.0001
        id = id + 1
    del k

    return sf_list

# -----------------------------
# Definition of Python toolbox.
class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Projection Selection Toolbox"
        self.alias = "ProjectionSelection"

        # List of tool classes associated with this toolbox
        self.tools = [Selection, Distortion, GeoCentre, Lattice]

# Definition of tool within toolbox to carry out Projection Selection only without distortion assessment.
class Selection(object):
    def __init__(self):
        self.label = "Projection Selection"
        self.description = "Outputs a list of candidate projections based on user input of Purpose and Footprint of Interest, "+\
                           "generating a projection definition (.prj) file for each candidate. The same selection "+\
                           "process undertaken by the Automated Projection Selection for ArcGIS (ArcMap Python Add-in) tool "+\
                           "is used, but in this case there is no quantitative distortion assessment and the user is unable "+\
                           "to interactively amend the data frame coordinate system to each candidate projection. The"+\
                           "projection definition (.prj) files can be imported as normal using the data frame properties dialog.\n\n"+\
                           "The user also selects an existing File GeoDatabase which is used for feature class objects "+\
                           "generated by this tool."
        self.canRunInBackground = False

    def getParameterInfo(self):

    # Purpose input parameter from drop-down list
        input_purpose = arcpy.Parameter(
            displayName="Purpose",
            name="input_purpose",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        input_purpose.filter.type = "ValueList"
        input_purpose.filter.list = ["General Reference", "Thematic Vector", "Thematic Raster", "Geospatial Analysis (Distance)",\
                                    "Geospatial Analysis (Area)","Navigation Routes (Geodesic)", "Navigation Routes (Loxodrome)",\
                                    "Ranges of Activity", "Flow Patterns", "World Index"]

    # Feature class or raster data which defines area of interest input parameter
        in_dataset = arcpy.Parameter(
            displayName="Input Footprint (Feature Class or Raster Dataset/Catalog which defines footprint of interest)",
            name="in_dataset",
            datatype="DEDatasetType",
            parameterType="Required",
            direction="Input")

    # Location for PRJ files generated for each candidate projection input parameter
        in_folder = arcpy.Parameter(
            displayName="Folder location for PRJ files generated for each candidate projection",
            name="in_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

    # Results text file output parameter
        out_results = arcpy.Parameter(
            displayName="Results Text File",
            name="out_results",
            datatype="DETextfile",
            parameterType="Required",
            direction="Output")

    # File geodatabase input parameter
        in_gdb = arcpy.Parameter(
            displayName="File Geodatabase for datasets created by tool",
            name="in_gdb",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        in_gdb.filter.list = ["Local Database"]

        params = [input_purpose, in_dataset, in_folder, out_results, in_gdb]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

    # Set variables from input parameters.
        inPurpose = parameters[0].value
        inDataset = parameters[1].value
        inFolder = parameters[2].value
        outResults = parameters[3].value
        inGdb = parameters[4].value

    # Create list of appropriate data types for each purpose.
        if inPurpose == "General Reference":
            data_types = ["Polygon"]
        elif inPurpose == "Thematic Vector":
            data_types = ["Polygon"]
        elif inPurpose == "Thematic Raster":
            data_types = ["RasterDataset", "RasterCatalog"]
        elif inPurpose == "Geospatial Analysis (Distance)":
            data_types = ["Point", "Polygon", "Polyline", "RasterDataset"]
        elif inPurpose == "Geospatial Analysis (Area)":
            data_types = ["Polygon", "RasterDataset"]
        elif inPurpose == "Navigation Routes (Geodesic)":
            data_types = ["Polyline"]
        elif inPurpose == "Navigation Routes (Loxodrome)":
            data_types = ["Polyline"]
        elif inPurpose == "Ranges of Activity":
            data_types = ["Point"]
        elif inPurpose == "Flow Patterns":
            data_types = ["Point", "Polygon"]
        elif inPurpose == "World Index":
            data_types = ["Polygon"]

    # Set input spatial reference and validate that the selected input dataset is referenced to a geographic coordinate system.
        desc_inDataset = arcpy.Describe(inDataset)
        inDataset_sr = desc_inDataset.spatialReference
        if inDataset_sr.type != "Geographic":
            messages.addErrorMessage("Input dataset must be referenced to a Geographic Coordinate System.")
            raise arcpy.ExecuteError

    # Create file path strings for input geodatabase (for temp feature classes) and input folder (for generated prj files).
        desc_gdb = arcpy.Describe(inGdb)
        gdb_path = desc_gdb.catalogPath
        desc_folder = arcpy.Describe(inFolder)
        folder_path = desc_folder.catalogPath

    # For selected input layer, either make copy of feature class in GDB or use raster_polygon function to create extent layer. Set data type and feature
    # count variables required by selection tools. Validate that an input Feature Class has features, and one feature where necessary for Purpose.
        if desc_inDataset.dataType == "FeatureClass":
            footprint_fc = os.path.join(gdb_path, "footprint_fc")
            arcpy.CopyFeatures_management(inDataset, footprint_fc)
            data_type = desc_inDataset.shapeType
            feature_count = int(arcpy.GetCount_management(inDataset).getOutput(0))
            if feature_count == 0:
                messages.addErrorMessage("Input dataset has no features.")
                raise arcpy.ExecuteError
            elif feature_count > 1:
                if purpose == "General Reference" or purpose == "Thematic Vector":
                    messages.addErrorMessage("For this Purpose the geographic footprint characteristics can only be determined adequately from a single "\
                                            "polygon which completely encompasses the desired footprint of interest. Please select an input dataset which "\
                                            "has a single polygon feature.")
                    raise arcpy.ExecuteError
        elif desc_inDataset.dataType == "RasterDataset" or desc_inDataset.dataType == "RasterCatalog":
            footprint_fc = raster_polygon(inDataset, gdb_path, inDataset_sr)
            data_type = desc_inDataset.dataType
            feature_count = 1
        else:
            messages.addErrorMessage("Input dataset is invalid.")
            raise arcpy.ExecuteError

    # Validate that the input dataset is of appropriate data type for the selected purpose.
        if data_type not in data_types:
            messages.addErrorMessage("Input dataset for this Purpose must be of type: {0}".format(data_types))
            raise arcpy.ExecuteError

    # Run projection selection function to determine geographic footprint characteristics and candidate projection options to enable creation of
    # spatial reference objects and .prj files for coordinate system selection drop-down list.
        candidates, proj_notes, input_type, output_size_str, input_lon, input_lat, input_lat_max, input_lat_min, input_k, output_sp_str, \
        output_extent_str, input_line, input_2pt = projection_selection(gdb_path, inPurpose, footprint_fc, data_type, feature_count, inDataset_sr)
        pcs_list = []
        pcs_id = 1
        for candidate in candidates:
            pcs_sr, pcs_string = candidate.create_sr()
            pcs_prj = "Candidate_PCS_{0}.prj".format(str(pcs_id))
            pcs_prj_filename = os.path.join(folder_path, pcs_prj)
            pcs_prj_file = open(pcs_prj_filename, "w")
            pcs_prj_file.write(pcs_string)
            pcs_prj_file.close()
            pcs_details = (pcs_sr, pcs_id, pcs_string, pcs_prj)
            pcs_list.append(pcs_details)
            pcs_id = pcs_id + 1
        del candidate

    # Set up results text file and write standard input information.
        result_filename = arcpy.Describe(outResults).catalogPath
        result_file = open(result_filename, "w")
        result_file.write("Selected Purpose = {0}\nInput Footprint Dataset = {1}\nOutpue GDB = {2}\n\n"\
                          .format(inPurpose, inDataset, gdb_path))

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

    # Write candidate options and projection notes to results text file.
        result_file.write("\nProjection selection options (see User Guide for selection diagrams and references):\n")
        options_str = ""
        for pcs in pcs_list:
            option_str = "   Candidate {0} = {1} ({2})\n".format(str(pcs[1]), pcs[3], pcs[0].name)
            options_str = options_str + option_str
        del pcs
        result_file.write(options_str)
        if proj_notes != []:
            result_file.write("\nProjection notes:\n")
            for note in proj_notes:
                result_file.write(note)
            del note

    # Close the results text file.
        result_file.close()

        return

# Definition of tool within toolbox to carry out Distortion Assessment only from input projection definition.
class Distortion(object):
    def __init__(self):
        self.label = "Distortion Assessment"
        self.description = "Enables distortion assessment of projection options not included as a candidate from the "+\
                           "Automated Projection Selection for ArcGIS (ArcMap Python Add-in) tool. Outputs a quantitative "+\
                           "assessment measure for an input projection definition (.prj file) and user input of Purpose "+\
                           "and Footprint of Interest. The user also selects an existing File GeoDatabase which is used for "+\
                           "feature class objects generated by this tool."
        self.canRunInBackground = False

    def getParameterInfo(self):

    # Projection definition input parameter
        in_prj = arcpy.Parameter(
            displayName="Projection Definition (.prj file)",
            name="in_prj",
            datatype="DEPrjFile",
            parameterType="Required",
            direction="Input")

    # Purpose input parameter from drop-down list
        input_purpose = arcpy.Parameter(
            displayName="Purpose",
            name="input_purpose",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        input_purpose.filter.type = "ValueList"
        input_purpose.filter.list = ["General Reference", "Thematic Vector", "Thematic Raster", "Geospatial Analysis (Distance)",\
                                    "Geospatial Analysis (Area)","Navigation Routes (Geodesic)", "Navigation Routes (Loxodrome)",\
                                    "Ranges of Activity", "Flow Patterns", "World Index"]

    # Feature class or raster data which defines footprint of interest input parameter
        in_dataset = arcpy.Parameter(
            displayName="Input Footprint (Feature Class or Raster Dataset/Catalog which defines footprint of interest)",
            name="in_dataset",
            datatype="DEDatasetType",
            parameterType="Required",
            direction="Input")

    # Results text file output parameter
        out_results = arcpy.Parameter(
            displayName="Results Text File",
            name="out_results",
            datatype="DETextfile",
            parameterType="Required",
            direction="Output")

    # File geodatabase input parameter
        in_gdb = arcpy.Parameter(
            displayName="File Geodatabase for datasets created by tool",
            name="in_gdb",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        in_gdb.filter.list = ["Local Database"]

        params = [in_prj, input_purpose, in_dataset, out_results, in_gdb]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

    # Set variables from input parameters.
        inPrj = parameters[0].value
        inPurpose = parameters[1].value
        inDataset = parameters[2].value
        outResults = parameters[3].value
        inGdb = parameters[4].value

    # Create spatial reference object from input prj file.
        input_sr = arcpy.SpatialReference()
        input_sr.createFromFile(inPrj)

    # Create list of appropriate data types and set distortion weights (distance / area / shape) for each purpose.
        if inPurpose == "General Reference":
            data_types = ["Polygon"]
            dist_weights = (0, 1, 1)
        elif inPurpose == "Thematic Vector":
            data_types = ["Polygon"]
            dist_weights = (0, 1, 0)
        elif inPurpose == "Thematic Raster":
            data_types = ["RasterDataset", "RasterCatalog"]
            dist_weights = (0, 1, 0)
        elif inPurpose == "Geospatial Analysis (Distance)":
            data_types = ["Point", "Polygon", "Polyline", "RasterDataset"]
            dist_weights = (1, 0, 0)
        elif inPurpose == "Geospatial Analysis (Area)":
            data_types = ["Polygon", "RasterDataset"]
            dist_weights = (0, 1, 0)
        elif inPurpose == "Navigation Routes (Geodesic)":
            data_types = ["Polyline"]
            dist_weights = (0, 1, 1)
        elif inPurpose == "Navigation Routes (Loxodrome)":
            data_types = ["Polyline"]
            dist_weights = (0, 1, 1)
        elif inPurpose == "Ranges of Activity":
            data_types = ["Point"]
            dist_weights = (1, 0, 0)
        elif inPurpose == "Flow Patterns":
            data_types = ["Point", "Polygon"]
            dist_weights = (0, 1, 1)
        elif inPurpose == "World Index":
            data_types = ["Polygon"]
            dist_weights = (0, 1, 1)

    # Set input spatial reference and validate that the selected input dataset is referenced to a geographic coordinate system.
        desc_inDataset = arcpy.Describe(inDataset)
        inDataset_sr = desc_inDataset.spatialReference
        if inDataset_sr.type != "Geographic":
            messages.addErrorMessage("Input dataset must be referenced to a Geographic Coordinate System.")
            raise arcpy.ExecuteError

    # Create file path string for input geodatabase.
        desc_gdb = arcpy.Describe(inGdb)
        gdb_path = desc_gdb.catalogPath

    # For selected input layer, either make copy of feature class in GDB or use raster_polygon function to create extent layer. Set data type and feature
    # count variables required by selection tools. Validate that an input Feature Class has features, and one feature where necessary for Purpose.
        if desc_inDataset.dataType == "FeatureClass":
            footprint_fc = os.path.join(gdb_path, "footprint_fc")
            arcpy.CopyFeatures_management(inDataset, footprint_fc)
            data_type = desc_inDataset.shapeType
            feature_count = int(arcpy.GetCount_management(inDataset).getOutput(0))
            if feature_count == 0:
                messages.addErrorMessage("Input dataset has no features.")
                raise arcpy.ExecuteError
            elif feature_count > 1:
                if purpose == "General Reference" or purpose == "Thematic Vector":
                    messages.addErrorMessage("For this Purpose the geographic footprint characteristics can only be determined adequately from a single "\
                                            "polygon which completely encompasses the desired footprint of interest. Please select an input dataset which "\
                                            "has a single polygon feature.")
                    raise arcpy.ExecuteError
        elif desc_inDataset.dataType == "RasterDataset" or desc_inDataset.dataType == "RasterCatalog":
            footprint_fc = raster_polygon(inDataset, gdb_path, inDataset_sr)
            data_type = desc_inDataset.dataType
            feature_count = 1
        else:
            messages.addErrorMessage("Input dataset is invalid.")
            raise arcpy.ExecuteError

    # Validate that the input dataset is of appropriate data type for the selected purpose.
        if data_type not in data_types:
            messages.addErrorMessage("Input dataset for this Purpose must be of type: {0}".format(data_types))
            raise arcpy.ExecuteError

    # For Geo Analysis Area purpose, reset to relevant Thematic purpose depending on data type.
        if inPurpose == "Geospatial Analysis (Area)":
            if data_type == "RasterDataset":
                inPurpose = "Thematic Raster"
            else:
                inPurpose = "Thematic Vector"
                if feature_count != 1:  # If input has multiple polygons create single input using Minimum Bounding Geometry tool.
                    footprint_mbg = os.path.join(gdb_path, "footprint_mbg")
                    arcpy.MinimumBoundingGeometry_management(footprint_fc, footprint_mbg, "CONVEX_HULL", "ALL")
                    footprint_fc = footprint_mbg

   	# Determine input geographical characteristics in terms of size / location / extent and other special types for Purpose based on input dataset.
        if inPurpose == "Geospatial Analysis (Distance)":
            input_type, input_lon, input_lat, input_lat_max, input_lat_min = location_type(inPurpose, footprint_fc, gdb_path, inDataset_sr)
        elif inPurpose == "Navigation Routes (Geodesic)" or inPurpose == "Navigation Routes (Loxodrome)":
            single_start, input_lon, input_lat = nav_condition(footprint_fc, feature_count)
            if feature_count == 1:
                input_type = "1 Line"
                arcpy.AddGeometryAttributes_management(footprint_fc, "LENGTH_GEODESIC; LINE_START_MID_END")
                cursor = arcpy.da.SearchCursor(footprint_fc,["LENGTH_GEO", "START_X", "START_Y", "MID_X", "MID_Y", "END_X", "END_Y"])
                for row in cursor:
                    length_geo = row[0]
                    lon1, lat1, lon2, lat2 = row[1], row[2], row[5], row[6]
                    if lon2 > 180:   # Reset longitude for geodesic line features created using the XY to Line function which cross 180.
                        lon2 = lon2 - 360
                    centre_lon, centre_lat = row[3], row[4]
                del row, cursor
                input_line = (lon1, lat1, centre_lon, centre_lat, lon2, lat2, length_geo)
            elif single_start == "TRUE":
                input_type = "Single Start Point"
            else:
                input_type, input_lon, input_lat, input_lat_max, input_lat_min = location_type(inPurpose, footprint_fc, gdb_path, inDataset_sr)
        elif inPurpose == "Ranges of Activity":
            if feature_count == 1:
                input_type = "1 Feature"
            else:
                input_type, input_lon, input_lat, input_lat_max, input_lat_min = location_type(inPurpose, footprint_fc, gdb_path, inDataset_sr)
        elif inPurpose == "Flow Patterns":
            if feature_count == 1:
                input_type = "1 Feature"
                cursor = arcpy.da.SearchCursor(footprint_fc, ["SHAPE@TRUECENTROID"])
                for row in cursor:
                    input_lon, input_lat = row[0]
                del row, cursor
            elif feature_count == 2:
                input_type = "2 Features"
                array = []
                cursor = arcpy.da.SearchCursor(footprint_fc, ["SHAPE@TRUECENTROID"])
                for row in cursor:
                    point = row[0]
                    array.append(point)
                del row, cursor
                lon1, lat1 = array[0]
                lon2, lat2 = array[1]
                input_2pt = (lon1, lat1, lon2, lat2)
            else:
                input_type, input_lon, input_lat, input_lat_max, input_lat_min = location_type(inPurpose, footprint_fc, gdb_path, inDataset_sr)
        elif inPurpose == "World Index":
            input_type = "World"
        else:
            input_type, input_size, input_lon, input_lat, input_lat_max, input_lat_min = footprint_type(footprint_fc, gdb_path, inDataset_sr)

    # Run distortion features function to generate dataset of fibonacci lattice points then run distortion assessment measures function for input prj file.
    # Create list of output details needed for results file and coordinate system selection drop-down list.
        distortion_fc, point_count = distortion_features(input_type, footprint_fc, data_type, gdb_path, inDataset_sr)
        dist = distortion_measures(distortion_fc, gdb_path, input_sr, point_count, dist_weights, inDataset_sr)

    # Set up results text file and write standard input information.
        result_filename = arcpy.Describe(outResults).catalogPath
        result_file = open(result_filename, "w")
        result_file.write("Input Projection Definition (.prj) File = {0}\nSelected Purpose = {1}\nInput Footprint Dataset = {2}\nOutpue GDB = {3}\n\n"\
                          .format(inPrj, inPurpose, inDataset, gdb_path))

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

    # Write distortion assessment outputs to results text file. For purposes using area distortion only use 10 decimal places to distinguish results, else use 5 decimal places.
        result_file.write("\nNumber of points selected for distortion assessment = {0}\n".format(str(point_count)))
        result_file.write("Distortion assessment weights (as defined by Purpose): Distance = {0}; Area = {1}; Shape = {2}\n\n"\
                          .format(str(dist_weights[0]), str(dist_weights[1]), str(dist_weights[2])))
        result_file.write("Distortion assessment results:\n\n   Distortion Measures:\n")
        if dist_weights == (0, 1, 0):
            result_file.write("     Distance:   Not measured\n")
            result_file.write("     Area:       EA = {0:.10f}\n".format(dist[1]))
            result_file.write("     Shape:      Not measured\n")
            result_file.write("     Combined Index = {0:.10f}\n\n".format(dist[4]))
        else:
            if input_purpose.dist_weights[0] == 0:
                result_file.write("     Distance:   Not measured\n")
            else:
                result_file.write("     Distance:   EP = {0:.5f}\n".format(dist[0]))
            if input_purpose.dist_weights[1] == 0:
                result_file.write("     Area:       Not measured\n")
            else:
                result_file.write("     Area:       EA = {0:.5f}\n".format(dist[1]))
            if input_purpose.dist_weights[2] == 0:
                result_file.write("     Shape:      Not measured\n")
            else:
                result_file.write("     Shape:   Index = {0:.5f}; No of points selected for shape calculation = {1}\n".format(dist[2], str(dist[3])))
            result_file.write("     Combined Index = {0:.5f}\n\n".format(dist[4]))

    # Close the results text file.
        result_file.close()

        return

class GeoCentre(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Geographic Centre"
        self.description = "Creates a point feature class for the approximate geographic centre of an input "+\
                           "dataset which defines a footprint of interest. Method is based on the iterative approach "+\
                           "by Rogerson (2015) using the Azimuthal Equidistant projection."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

    # Feature class or raster data which defines footprint of interest
        in_dataset = arcpy.Parameter(
            displayName="Input Footprint (Feature Class or Raster Dataset/Catalog which defines footprint of interest)",
            name="in_dataset",
            datatype="DEDatasetType",
            parameterType="Required",
            direction="Input")

    # File geodatabase input parameter
        in_gdb = arcpy.Parameter(
            displayName="File Geodatabase for datasets created by tool",
            name="in_gdb",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        in_gdb.filter.list = ["Local Database"]

    # Name of output geo centre feature class parameter
        out_fc = arcpy.Parameter(
            displayName="Name of geographic centre point feature class created by tool",
            name="out_fc",
            datatype="GPString",
            parameterType="Required",
            direction="Output")

        params = [in_dataset, in_gdb, out_fc]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

    # Set variables from input parameters.
        inDataset = parameters[0].value
        inGdb = parameters[1].value
        outFc = parameters[2].valueAsText

    # Set input spatial reference and validate that the selected input dataset is referenced to a geographic coordinate system.
        desc_inDataset = arcpy.Describe(inDataset)
        inDataset_sr = desc_inDataset.spatialReference
        if inDataset_sr.type != "Geographic":
            messages.addErrorMessage("Input dataset must be referenced to a Geographic Coordinate System.")
            raise arcpy.ExecuteError

    # Create file path string for input geodatabase.
        desc_gdb = arcpy.Describe(inGdb)
        gdb_path = desc_gdb.catalogPath
        env.workspace = gdb_path

    # For selected input layer, either make copy of feature class in GDB or use raster_polygon function to create extent layer. Set data type and feature
    # count variables, and validate that an input Feature Class has features.
        if desc_inDataset.dataType == "FeatureClass":
            input_fc = os.path.join(gdb_path, "input_fc")
            arcpy.CopyFeatures_management(inDataset, input_fc)
        elif desc_inDataset.dataType == "RasterDataset" or desc_inDataset.dataType == "RasterCatalog":
            input_fc = raster_polygon(inDataset, gdb_path, inDataset_sr)
        else:
            messages.addErrorMessage("Input dataset is invalid.")
            raise arcpy.ExecuteError
        feature_count = int(arcpy.GetCount_management(input_fc).getOutput(0))
        if feature_count == 0:
            messages.addErrorMessage("Input dataset has no features.")
            raise arcpy.ExecuteError

    # For single feature input, assumed to be polygon, read input feature centroid lat/long as initial approximation.
        elif feature_count == 1:
            cursor = arcpy.da.SearchCursor(input_fc, ["SHAPE@"])
            for row in cursor:
                centre = row[0].trueCentroid
                lon_centre_gcs, lat_centre_gcs  = centre.X, centre.Y
            del row, cursor

    # For multiple feature input, determine approximate mean centre from GCS input features as initial approximation.
        else:
            mean_centre_gcs = os.path.join(gdb_path, "mean_centre_gcs")
            arcpy.MeanCenter_stats(input_fc, mean_centre_gcs)
            cursor = arcpy.da.SearchCursor(mean_centre_gcs, ["SHAPE@XY"])
            for row in cursor:
                lon_centre_gcs, lat_centre_gcs = row[0]
            del row, cursor

    # Use approximate centroid to project input feature using Azimuthal Equidistant projection based on that location.
        equidistant_fc = os.path.join(gdb_path, "equidistant_fc")
        pcs_equidistant = ProjectedCS("B", "Azimuthal_Equidistant", inDataset_sr, lon_centre_gcs, lat_centre_gcs, "", "", "", "")
        centroid_sr, equidistant_string = pcs_equidistant.create_sr()
        arcpy.Project_management(input_fc, equidistant_fc, centroid_sr)

    # For single feature input, determine centroid of projected feature and create new feature class containing that point using Azimuthal Equidistant projection.
        if feature_count == 1:
            cursor = arcpy.da.SearchCursor(equidistant_fc, ["SHAPE@TRUECENTROID"])
            for row in cursor:
                x_centre, y_centre = row[0]
            del row, cursor
            arcpy.CreateFeatureclass_management(gdb_path, "equidistant_centroid", "POINT", spatial_reference = centroid_sr)
            equidistant_centroid = os.path.join(gdb_path, "equidistant_centroid")
            centroid = (x_centre, y_centre)
            cursor = arcpy.da.InsertCursor(equidistant_centroid, ["SHAPE@XY"])
            cursor.insertRow([centroid])
            del cursor

    # For multiple feature input, recalculate mean centre using projected features.
        else:
            equidistant_centroid = os.path.join(gdb_path, "equidistant_centroid")
            arcpy.MeanCenter_stats(equidistant_fc, equidistant_centroid)

    # Project to GCS and read lat/long as better approximation of geographic centre of input dataset.
        geo_centre = os.path.join(gdb_path, outFc)
        arcpy.Project_management(equidistant_centroid, geo_centre, inDataset_sr)
        cursor = arcpy.da.SearchCursor(geo_centre, ["SHAPE@XY"])
        for row in cursor:
            lon_centre, lat_centre = row[0]
        del row, cursor

    # Add field for lat/long and populate with geo centre values
        arcpy.AddField_management(geo_centre, "Longitude", "DOUBLE")
        arcpy.AddField_management(geo_centre, "Latitude", "DOUBLE")
        cursor = arcpy.da.UpdateCursor(geo_centre, ["Longitude", "Latitude"])
        for row in cursor:
            cursor.updateRow([lon_centre, lat_centre])
        del row, cursor

        messages.addMessage("See {0} for {1} point feature class which is calculated for the "\
                            "input dataset {2}.\nLatitude = {3}\nLongitude = {4}"\
                            .format(desc_gdb.name, outFc, desc_inDataset.name, lat_centre, lon_centre))

        return

class Lattice(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fibonacci Lattice"
        self.description = "Creates a feature class of points which are equally spaced using a fibonaaci lattice approach, see Gonzalez (2010) and Baselga (2018)."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

    # Input value for N which defines how many points will be generated
        in_n = arcpy.Parameter(
            displayName="N value (number of points generated = 2N + 1)",
            name="in_n",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

    # Spatial reference input parameter
        in_sr = arcpy.Parameter(
            displayName="Geographic coordinate system of lattice points generated by tool",
            name="in_sr",
            datatype="GPSpatialReference",
            parameterType="Required",
            direction="Input")

    # File geodatabase input parameter
        in_gdb = arcpy.Parameter(
            displayName="File Geodatabase for lattice points feature class created by tool",
            name="in_gdb",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        in_gdb.filter.list = ["Local Database"]

    # Name of output lattice points feature class parameter
        out_fc = arcpy.Parameter(
            displayName="Name of lattice points feature class created by tool",
            name="out_fc",
            datatype="GPString",
            parameterType="Required",
            direction="Output")

        params = [in_n, in_sr, in_gdb, out_fc]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

    # Set variables from input parameters.
        N = parameters[0].value
        inSr = parameters[1].value
        inGdb = parameters[2].value
        outFc = parameters[3].valueAsText

    # Set input spatial reference and validate that the selected input dataset is referenced to a geographic coordinate system.
        if inSr.type != "Geographic":
            messages.addErrorMessage("Input dataset must be referenced to a Geographic Coordinate System.")
            raise arcpy.ExecuteError

    # Create file path string for input geodatabase.
        desc_gdb = arcpy.Describe(inGdb)
        gdb_path = desc_gdb.catalogPath
        env.workspace = gdb_path

    # Create feature class and add points based on algorithm by Gonzalez (2010) and Baselga (2018).
        lattice_pts = os.path.join(gdb_path, outFc)
        arcpy.CreateFeatureclass_management(gdb_path, outFc, "Point", spatial_reference = inSr)
        P = (2 * N) + 1
        phi = (1 + math.sqrt(5)) / 2   # Calculates the Golden Ratio.
        i = 0 - N
        cursor = arcpy.da.InsertCursor(lattice_pts, ["SHAPE@XY"])
        while i <= N:
            lat_rad = math.asin((2 * i) / float(P))
            lat_deg = math.degrees(lat_rad)
            lon_rad = (2 * math.pi * math.fmod(i, phi)) / phi
            lon_deg = math.degrees(lon_rad)
            if lon_deg > 180.0:
                lon_deg = lon_deg - 360.0
            elif lon_deg < -180.0:
                lon_deg = lon_deg + 360.0
            lattice_pt = arcpy.Point(lon_deg, lat_deg)
            cursor.insertRow([lattice_pt])
            i = i + 1
        del cursor

        return
