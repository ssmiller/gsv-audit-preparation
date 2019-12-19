"""
Start to (hopefully) finish process to generate street segments, routes,
and combined segment and route, shortened to maximum distance.
These will be used for Google Audits in REGARDSBSE.
Created by: Stephanie Miller
Last revision: 2019-12-18 temporary remove of extra code; may turn into tests
Biggest TODO -- fix feature class naming conventions and locations
Adjust some manual steps if possible
"""

# Import packages
import os
import random
import csv
import itertools
import cmath
import arcpy
from gsv.participant import Participant

# Check out Network Analyst license if available.
# Fail if the Network Analyst license is not available.
if arcpy.CheckExtension("network") == "Available":
    arcpy.CheckOutExtension("network")
else:
    raise arcpy.ExecuteError("Network Analyst Extension license is not available.")

# Global variables
# spatial reference for input geocodes or business locations
sr = arcpy.SpatialReference(4326) #GCS WGS84
# spatial reference for reprojections
srCEC = arcpy.SpatialReference("USA Contiguous Equidistant Conic")
# network analysis info
networkgdb = r"C:\Users\smillers\Downloads\StreetmapPremium2017\FGDB\StreetMap_Data\NorthAmerica.gdb"
network = os.path.join(networkgdb, "Routing", "Routing_ND")
travel_mode="Walking Distance"





def getRouteToClosestDestination(businesscsv, incidents_cf, datesuffix, fcprefix):
    """Create routes to closest business (e.g. supermarket)"""
    businesslyr = "business_lyr"
    xfield_sma =  "adr_gis_xwgs84_x_2015" #WGS84
    yfield_sma = "adr_gis_ywgs84_x_2015" #WGS84
    print("Creating feature class from businesses")
    arcpy.MakeXYEventLayer_management(businesscsv, xfield_sma, yfield_sma, businesslyr, sr)
    businessfcnameWGS84 = fcprefix + "_" + datesuffix + "_WGS84"
    try:
        arcpy.CopyFeatures_management(businesslyr, businessfcnameWGS84)
        print(arcpy.GetMessages())
    except arcpy.ExecuteError as e:
        print(e)
    print("Creating Closest Facility Analysis Layer")
    result_object_cf = arcpy.na.MakeClosestFacilityAnalysisLayer(network,
        "Closest {}".format(fcprefix), travel_mode, line_shape="ALONG_NETWORK")
    # default travel from incident to facilities, no cutoff. Use default number of facilities to find (1)
    facilities_cf = businessfcnameWGS84
    #Get the layer object from the result object. The closest facility layer can
    #now be referenced using the layer object.
    layer_object_cf = result_object_cf.getOutput(0)
    #Get the names of all the sublayers within the closest facility layer.
    sublayer_names_cf = arcpy.na.GetNAClassNames(layer_object_cf)
    #Stores the layer names that we will use later
    facilities_layer_name_cf = sublayer_names_cf["Facilities"]
    incidents_layer_name_cf = sublayer_names_cf["Incidents"]
    #Load the supermarkets as Facilities using the default field mappings and
    # Map the Name property to the uhcid field
    print("Loading Locations")
    field_mappings_cf1 = arcpy.na.NAClassFieldMappings(
       layer_object_cf, facilities_layer_name_cf)
    field_mappings_cf1["Name"].mappedFieldName = "uhcid"
    arcpy.na.AddLocations(layer_object_cf, facilities_layer_name_cf,
                          facilities_cf, field_mappings_cf1, "", append="CLEAR")
    #Load the people as Incidents. Map the Name property from the kmzid field
    #using field mappings
    field_mappings_cf = arcpy.na.NAClassFieldMappings(
        layer_object_cf, incidents_layer_name_cf)
    field_mappings_cf["Name"].mappedFieldName = "kmzid"
    arcpy.na.AddLocations(layer_object_cf, incidents_layer_name_cf,
                          incidents_cf, field_mappings_cf, "", append="CLEAR")
    #Solve the closest facility layer
    print("Solving Closest Facility Layer")
    arcpy.na.Solve(layer_object_cf)
    print(arcpy.GetMessages())
    routes_layer_name_cf = sublayer_names_cf["CFRoutes"]
    # Export the route lines for further analysis
    print("Copying Routes to feature class")
    out_cfroutes = "Closest{}_googleaudit{}".format(fcprefix, datesuffix)
    arcpy.management.CopyFeatures(routes_layer_name_cf, out_cfroutes)
    # Add field to output routes
    print("Adding kmzid to routes")
    arcpy.AddField_management(out_cfroutes, "kmzid", "TEXT", field_length=10)
    print(arcpy.GetMessages())
    arcpy.CalculateField_management(out_cfroutes, "kmzid", "!Name!.split(' - ')[0]")
    print(arcpy.GetMessages())
    print("Adding uhcid to routes")
    arcpy.AddField_management(out_cfroutes, "uhcid", "TEXT", field_length=20)
    print(arcpy.GetMessages())
    arcpy.CalculateField_management(out_cfroutes, "uhcid", "!Name!.split(' - ')[1]")
    print(arcpy.GetMessages())
    print("Reprojecting routes")
    # Reproject to USCEC
    cfroute_cec = out_cfroutes + "_cec"
    arcpy.Project_management(out_cfroutes, cfroute_cec, srCEC)
    print(arcpy.GetMessages())
    # TODO add other return data if needed
    # This is all routes created for the users.
    return cfroute_cec #pass back the name of the route output


def checklinedir(inputLine_lyr, comparisonEndpoints_dir, fieldtoupdate, disttolerance):
    with arcpy.da.UpdateCursor(
        inputLine_lyr, ["kmzid", "START_X", "START_Y", fieldtoupdate]) as cursor:
        for row in cursor:
            print(row)
            kmz = row[0]
            if kmz not in comparisonEndpoints_dir.keys():
                flip = 99
            else:
                startx_combinedrte = row[1]
                starty_combinedrte = row[2]
                startx_segment = comparisonEndpoints_dir[kmz]["startx"]
                starty_segment = comparisonEndpoints_dir[kmz]["starty"]
                endx_segment = comparisonEndpoints_dir[kmz]["endx"]
                endy_segment = comparisonEndpoints_dir[kmz]["endy"]
                if abs(startx_combinedrte - startx_segment) <= disttolerance and abs(starty_combinedrte - starty_segment) <= disttolerance:
                    flip = 0
                elif abs(startx_combinedrte - endx_segment) <= disttolerance and abs(starty_combinedrte - endy_segment) <= disttolerance:
                    flip = 0
                else:
                    flip = 1
            row[3] = flip
            cursor.updateRow(row)


def get_adjusted_segment_endpoints(gcloc, length_of_segment, halfmax, testcase):
    shiftseglst = [0,25,-25,40,-40,25, -25, 0]
    starting_shift = shiftseglst[testcase]
    new_start = gcloc - (halfmax + starting_shift)
    new_end = gcloc + (halfmax - starting_shift)
    while testcase < len(shiftseglst) and (new_start < 0 or new_end > length_of_segment):
        print("ERROR - bad START or END location ({}, {})".format(new_start, new_end))
        starting_shift = shiftseglst[testcase]
        new_start = gcloc - (halfmax + starting_shift)
        new_end = gcloc + (halfmax - starting_shift)
        testcase += 1
    return (new_start, new_end)


def geocode_in_segment(segtuple, gcdist):
    if gcdist <= segtuple[1] and gcdist >= segtuple[0]:
        return True
    else:
        return False


def addFlipFields(fc):
    arcpy.AddField_management(fc, "flip1", "SHORT")
    arcpy.AddField_management(fc, "flip2", "SHORT")
    return


def cut_second_route(rte1, rte2, k, max_dist=400):
    print("updating second route for kmzid {}".format(k))
    with arcpy.da.UpdateCursor(rte2, ['kmzid', "SHAPE@"]) as ucursor:
        for row in ucursor:
            route2feat = row[1]
            segfeat = [lrow[0]for lrow in arcpy.da.SearchCursor(rte1, ["SHAPE@"])][0]
            print("firstRouteLength: {}".format(segfeat.length))
            print("secondRouteLength: {}".format(route2feat.length))
            route2featdiff = route2feat.difference(segfeat)
            print("Route 2 length not overlapping with first route: {}".format(route2featdiff.length))
            r2lennew = max_dist - segfeat.length
            print(r2lennew)
            row[1] = route2featdiff.segmentAlongLine(0,r2lennew)
            ucursor.updateRow(row)
    print("completed update of {}".format(k))
    return


def createSideOfStreetIndicators(geocodes, roads):
    output_gdb = os.path.basename(arcpy.env.workspace)
    #out_location = arcpy.env.workspace
    out_location = os.path.dirname(arcpy.env.workspace)
    #ate_suf = output_gdb.split("_")[1].split(".")[0] #TODO revert!
    midpoint_csv = os.path.join(out_location, "googleAuditRoadSide_REGARDSBSE_{}.csv".format(date_suf))

    # ensure that the roads have the midpoint location and bearing
    arcpy.AddGeometryAttributes_management(roads, ["LINE_START_MID_END","LINE_BEARING"])
    # Get details for both the lines and points
    #pnts = dict([(p[0],(p[1], p[2])) for p in arcpy.da.SearchCursor(geocodes,["kmzid","SHAPE@","auditgr"])])
    pnts = dict([(p[0],(p[1])) for p in arcpy.da.SearchCursor(geocodes,["kmzid","SHAPE@"])]) #TODO REVERT
    lines = dict([(l[0], (l[1], l[2], l[3], l[4])) for l in arcpy.da.SearchCursor(roads,["kmzid","SHAPE@","BEARING", u"MID_X", u"MID_Y"])])

    # Set a distance that you want to shift the midpoint by, in units of the layer (meters)
    indic_distance = 10

    with open(midpoint_csv, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        #writer.writerow(['kmzid', 'auditgr', 'INDIC_X', 'INDIC_Y', "LineBearing", "PointBearing"])
        writer.writerow(['kmzid', 'INDIC_X', 'INDIC_Y', "LineBearing", "PointBearing"])
        for pl in itertools.product(pnts,lines):
            pntID = pl[0]
            lineID = pl[1]
            if str(pntID) == lineID:
                #auditgroup = pnts[pntID][1]
                linebearing = lines[lineID][1]
                qpdoutput = lines[lineID][0].queryPointAndDistance(pnts[pntID][0])
                rside = qpdoutput[3]
                if rside == True:
                    pointbearing = (linebearing + 90) % 360
                elif rside == False:
                    pointbearing = (linebearing - 90) % 360
                else:
                    print("issue with line bearing for kmzid {}".format(lineID))
                midx = lines[lineID][2]  # 'MIDPOINT_X'
                midy = lines[lineID][3]  # 'MIDPOINT_Y'
                #convert bearing to arithmetic angle in radians
                angle = 90 - pointbearing
                if angle < -180:
                    angle = 360 + angle
                angle = math.radians(angle)
                start = complex(midx, midy)
                movement = cmath.rect(indic_distance, angle)
                end = start + movement
                endx = end.real
                endy = end.imag
                #temprow = [pntID, auditgroup, endx, endy, linebearing, pointbearing]
                temprow = [pntID, endx, endy, linebearing, pointbearing]
                writer.writerow(temprow)
                print("Point {0} is on the right side of line {1}: {2} (bearing: {3})".format(pntID, lineID, rside, linebearing))

    x_coords = 'INDIC_X'
    y_coords = 'INDIC_Y'
    out_layer = "indicators_layer"

    saved_fc = 'streetside_indicator'

    # Make the XY event layer...
    arcpy.MakeXYEventLayer_management(midpoint_csv, x_coords, y_coords, out_layer, srCEC)

    # Print the total rows
    print("There are {} records created from the midpoint indicator CSV".format(arcpy.GetCount_management(out_layer)))

    # print("Fields in the layer:")
    # fields = arcpy.ListFields(out_layer)
    # for field in fields:
    #     print(field.name, field.type, arcpy.ValidateFieldName(field.name))

    # Execute FeatureClassToFeatureClass
    arcpy.FeatureClassToFeatureClass_conversion(out_layer, out_location, saved_fc)
    print("Created Midpoint Feature Class {}".format(saved_fc))
    print("Need to manually update indicator for 'ShortLoops'")
    return

def main(date_suf, maxdistance=400, nets_csv = r"M:\EPL-GEO_BSERE\Data\WorkingData\google.audit\NETS_Supermarkets_yr2014_REGARDSBSE_2019.04.04_v1_SM.csv", gccsv_file="googleAuditLatLong_REGARDSBSE_2019.06.26.csv"):
    """ Overall Process -- input the desired date suffix, distance, and point csv datasources (destinations, origins) """
    # Create working gdb
    output_dir = r"M:\EPL-GEO_BSERE\Data\WorkingData\google.audit\intermediate-output"
    outputgdb = "googleaudit_{}.gdb".format(date_suf)
    salinesfd = "ServiceAreaLines_CEC"
    nearFinalFD = "auditSegmentRoutes_nearFinal_CEC"
    if not arcpy.Exists(os.path.join(output_dir, outputgdb)):
        arcpy.CreateFileGDB_management(output_dir, outputgdb)
    arcpy.env.workspace = os.path.join(output_dir, outputgdb)
    arcpy.env.overwriteOutput = True
    if not arcpy.Exists(os.path.join(arcpy.env.workspace, nearFinalFD)):
        arcpy.CreateFeatureDataset_management(arcpy.env.workspace, nearFinalFD, srCEC)
    else:
        print("FD {} already exists".format(nearFinalFD))
    #xfield = "x"
    #yfield = "y"
    #gccsvfolder = r"M:\EPL-GEO_BSERE\Data\WorkingData\google.audit\original_data"
    #gccsv = os.path.join(gccsvfolder, gccsv_file)
    #gclayername = "googleauditGC_lyr"
    #arcpy.MakeXYEventLayer_management(gccsv, xfield, yfield, gclayername, sr)
    #gcfcnameWGS84 = "GoogleAuditGeocodes_" + date_suf + "_WGS84"
    #try:
    #    arcpy.CopyFeatures_management(gclayername, gcfcnameWGS84)
    #    print(arcpy.GetMessages())
    #except arcpy.ExecuteError as e:
    #    print(e)
    participants = Participant(r"M:\EPL-GEO_BSERE\Data\WorkingData\google.audit\original_data", sr)
    # set up list for IDs that will require some manual intervention
    adjustAuditRouteManually = list()
    ######################
    ### STREET SEGMENT ###
    ######################
    # Get the closest service area line as the individual's street segment
    # returns a result in USCEC projection
    # Use a larger service area than the maximum distance to ensure as much of the
    # individual's street segment is covered as possible.
    # There will be some segments which are too long and need to be cut down.
    segmentfc = getClosestLineFromServiceArea(maxdistance)
    gcCEC = "_".join(gcfcnameWGS84.split("_")[:-1] + ["cec"]) #TODO maybe just pass this name back from the function
    arcpy.MakeFeatureLayer_management(gcCEC, "geocodes_lyr")
    serviceareaCEC = "SAWDlines_googleaudit{}_{}m_dslv_cec".format(date_suf, int(maxdistance*1.5)) #TODO maybe just pass this name back from the function
    arcpy.MakeFeatureLayer_management(serviceareaCEC, "servicearea_lyr")
    # Add part count and start/end (and length?) values to segment
    arcpy.AddGeometryAttributes_management(
        segmentfc, ["PART_COUNT", "LINE_START_MID_END"])
    # Get the endpoints for each segment and store in a dictionary
    with arcpy.da.SearchCursor(segmentfc, ["kmzid", "START_X", "START_Y", "END_X", "END_Y"]) as cursor:
        segmentendpoints = {row[0]: {"startx": row[1], "starty": row[2], "endx": row[3], "endy": row[4]} for row in cursor}
    # Verify all geocodes have a segment
    if int(arcpy.GetCount_management(gcfcnameWGS84).getOutput(0)) != len(segmentendpoints.keys()):
        print("Some segments may be missing!")
    #### Reason 1 for multipart combined segments and routes: The segment is a loop ###
    issue1_segmentloops = []
    for k, v in segmentendpoints.items():
        if v['startx'] == v['endx'] and v['starty'] == v['endy']:
            v['startendsame'] = 1
            issue1_segmentloops.append(k)
        else:
            v['startendsame'] = 1

    # Limit based on original segment length?
    segmentlengths = {}
    addroutes = []
    cutsegment = []
    segnochange = []
    with arcpy.da.SearchCursor(segmentfc, ["kmzid", "Shape_Length"]) as cursor:
        for row in cursor:
            segmentlengths[row[0]] = row[1]
            if row[1] < maxdistance - 1:
                addroutes.append(row[0])
            elif row[1] > maxdistance:
                cutsegment.append(row[0])
            else:
                segnochange.append(row[0])

    #####################
    ### LONG SEGMENTS ###
    #####################
    # Cut down long segments that do not form loops
    # copy the segmentfc to a new feature class before processing
    segmentfc_copy = segmentfc + "_max{}m".format(maxdistance)
    arcpy.CopyFeatures_management(segmentfc, segmentfc_copy)

    halfmaxdist = maxdistance//2
    with arcpy.da.UpdateCursor(segmentfc_copy, ["kmzid", "Shape_Length", "SHAPE@"], spatial_reference=srCEC) as cursor:
        random.seed(30) # set a seed to enable the initial value chosen for the segment offset
        for row in cursor:
            # select the associated line
            k = row[0]
            if k in cutsegment and k not in issue1_segmentloops:
                within_halfmax_of_start = False
                within_halfmax_of_end = False
                within_max_of_start = False
                within_max_of_end = False
                #print(type(row[2]))
                linefeat = row[2]
                linelen = linefeat.length
                # print("kmzid: {}; Line Length: {}; Length from Shape_Length: {}; These match: {}".format(k, linelen, row[1], linelen==row[1]))
                # Get geocode geometry - assumes only one feature returned
                geocode = arcpy.da.SearchCursor("geocodes_lyr", "SHAPE@", "{} = '{}'".format(arcpy.AddFieldDelimiters("geocodes_lyr", "kmzid"),k),spatial_reference=srCEC).next()[0] #TODO REVERT!!!
                # Get the distance this point is along the line
                qpd = linefeat.queryPointAndDistance(geocode)
                # lineoffset = qpd[2]
                distalong = qpd[1]
                if distalong - halfmaxdist < 0:
                    within_halfmax_of_start = True
                if linelen - distalong < halfmaxdist:
                    within_halfmax_of_end = True
                if distalong - maxdistance < 0:
                    within_max_of_start = True
                if linelen - distalong < maxdistance:
                    within_max_of_end = True
                print("""Geocode {} is located {}m along the line, which is {} in total length.""".format(k, round(distalong,2), round(linelen,2)))
                if within_halfmax_of_start == True and within_halfmax_of_end == True:
                    print("ERRROR! {} within halfmax of both ends -- should not be cut".format(k))
                elif within_halfmax_of_start == True or within_max_of_start == True:
                    startpt = 0
                    endpt = maxdistance
                elif within_halfmax_of_end == True or within_max_of_end == True:
                    startpt = linelen - maxdistance
                    endpt = linelen
                else: #(if within_halfmax_of_start == False and within_halfmax_of_end == False:)
                    # really long segment, probably rural area
                    # randomize the potential distances so that it's not always 200m on either side of participant
                    # will need to do additional tests each time to ensure this works...
                    randcase = random.randint(0,4)
                    startpt, endpt = get_adjusted_segment_endpoints(distalong, linelen, halfmaxdist,randcase)
                print("Suggested segment: ({}, {})".format(round(startpt, 2), round(endpt, 2)))
                print("TEST 1: Point within segment: {}".format(distalong <= endpt and distalong >= startpt))
                print("TEST 2: segment ~400m: {}".format(endpt - startpt <= maxdistance + 1 and endpt - startpt >= maxdistance - 1))
                row[2] = linefeat.segmentAlongLine(startpt, endpt)
                #print(type(row[2]))
                cursor.updateRow(row)
    # Cut down long segments that form loops
    with arcpy.da.UpdateCursor(segmentfc_copy, ["kmzid", "Shape_Length", "SHAPE@"], spatial_reference=srCEC) as cursor:
        random.seed(30) # set a seed to enable the initial value chosen for the segment offset
        for row in cursor:
            test1, test2, test3, test4a, test4b = False, False, False, False, False
            test5, test6a, test6b = False, False, False
            shortened_segment = ""
            k = row[0]
            if k in cutsegment and k in issue1_segmentloops:
                geocode = arcpy.da.SearchCursor("geocodes_lyr", "SHAPE@", "{} = '{}'".format(arcpy.AddFieldDelimiters("geocodes_lyr", "kmzid"),k),spatial_reference=srCEC).next()[0] #TODO REVERT!!!
                linefeat = row[2]
                print(linefeat)
                linelen = linefeat.length
                qpd = linefeat.queryPointAndDistance(geocode)
                distalong = qpd[1]
                # Where is the intersection located where this loop rejoins the street network?
                arcpy.SelectLayerByAttribute_management("servicearea_lyr", "NEW_SELECTION", """{} = '{}'""".format(arcpy.AddFieldDelimiters("servicearea_lyr", "kmzid"), k))
                # subset to only the items that aren't the same length (should eliminate the segment from the service area selection)
                arcpy.SelectLayerByAttribute_management("servicearea_lyr", "REMOVE_FROM_SELECTION", """{} = {}""".format(arcpy.AddFieldDelimiters("servicearea_lyr", "SHAPE_LENGTH"), linelen))
                with arcpy.da.SearchCursor("servicearea_lyr", "SHAPE@") as sa_cursor:
                    potential_intersects = []
                    for srow in sa_cursor:
                        roadfeat = srow[0]
                        if linefeat.disjoint(roadfeat) is False:
                            potential_intersects.append(linefeat.intersect(roadfeat, 1))
                first_intersect = potential_intersects[0]
                if first_intersect.type == "multipoint": # get the first point from the multipoint array
                    first_intersect = first_intersect[0]
                intersectpt_dist = linefeat.measureOnLine(first_intersect)
                print("Looped segment {} intersects service area at {}m along the line which is {}m long. The geocode is at {}m along the line".format(k, round(intersectpt_dist,2), round(linelen,2), round(distalong, 2)))
                # Does the intersection occur at 0m along line?
                if intersectpt_dist == 0:
                    # start at intersection / start point and go maximum distance
                    potential_segment1 = (0, maxdistance)
                    # start at intersection / end point and go (negative) maximum distance
                    potential_segment2 = (linelen - maxdistance, linelen)
                    test1 = geocode_in_segment(potential_segment1, distalong)
                    test2 = geocode_in_segment(potential_segment2, distalong)
                    # print("Tests 1 and 2, Intersection at 0m along line: {} {}".format(test1, test2))
                    if test1:
                        print("Suggested segment: ({}, {})".format(round(potential_segment1[0], 2), round(potential_segment1[1], 2)))
                        shortened_segment = linefeat.segmentAlongLine(potential_segment1[0], potential_segment1[1])
                    elif test2:
                        print("Suggested segment: ({}, {})".format(round(potential_segment2[0], 2), round(potential_segment2[1], 2)))
                        shortened_segment = linefeat.segmentAlongLine(potential_segment2[0], potential_segment2[1])
                if shortened_segment == "" and intersectpt_dist + maxdistance <= linelen:
                    potential_segment3 = (intersectpt_dist, intersectpt_dist + maxdistance)
                    test3 = geocode_in_segment(potential_segment3, distalong)
                    # print("Test 3: {}".format(test3))
                    if test3:
                        print("Suggested segment: ({}, {})".format(round(potential_segment3[0], 2), round(potential_segment3[1], 2)))
                        shortened_segment = linefeat.segmentAlongLine(potential_segment3[0], potential_segment3[1])
                    # else:
                    #     print("test3 failed for {}".format(k))
                if shortened_segment == "":
                    potential_segment4a = (intersectpt_dist, linelen)
                    potential_segment4b = (0, maxdistance - (linelen - intersectpt_dist))
                    test4a = geocode_in_segment(potential_segment4a, distalong)
                    test4b = geocode_in_segment(potential_segment4b, distalong)
                    if test4a or test4b:
                        print("Suggested segment: ({}, {}) plus ({}, {})".format(
                            round(potential_segment4a[0], 2), round(potential_segment4a[1], 2),
                            round(potential_segment4b[0], 2), round(potential_segment4b[1], 2)))
                        partialseg_a = linefeat.segmentAlongLine(potential_segment4a[0], potential_segment4a[1])
                        partialseg_b = linefeat.segmentAlongLine(potential_segment4b[0], potential_segment4b[1])
                        fixed_segment_a = arcpy.Polyline(partialseg_a.getPart(0), partialseg_a.spatialReference)
                        fixed_segment_b = arcpy.Polyline(partialseg_b.getPart(0), partialseg_b.spatialReference)
                        shortened_segment = fixed_segment_a.union(fixed_segment_b)
                    elif intersectpt_dist - maxdistance > 0:
                        potential_segment5 = (intersectpt_dist - maxdistance, intersectpt_dist)
                        test5 = geocode_in_segment(potential_segment5, distalong)
                        if test5:
                            print("Suggested segment: ({}, {})".format(round(potential_segment5[0], 2), round(potential_segment5[1], 2)))
                            shortened_segment = linefeat.segmentAlongLine(potential_segment5[0], potential_segment5[1])
                        # else:
                        #     print("test5 failed for {}".format(k))
                    else:
                        potential_segment6a = (linelen - (maxdistance - intersectpt_dist), linelen)
                        potential_segment6b = (0, intersectpt_dist)
                        test6a = geocode_in_segment(potential_segment6a, distalong)
                        test6b = geocode_in_segment(potential_segment6b, distalong)
                        if test6a or test6b:
                            print("Suggested segment: ({}, {}) plus ({}, {})".format(
                                round(potential_segment6a[0], 2), round(potential_segment6a[1], 2),
                                round(potential_segment6b[0], 2), round(potential_segment6b[1], 2)))
                            partialseg_a = linefeat.segmentAlongLine(potential_segment6a[0], potential_segment6a[1])
                            partialseg_b = linefeat.segmentAlongLine(potential_segment6b[0], potential_segment6b[1])
                            fixed_segment_a = arcpy.Polyline(partialseg_a.getPart(0), partialseg_a.spatialReference)
                            fixed_segment_b = arcpy.Polyline(partialseg_b.getPart(0), partialseg_b.spatialReference)
                            shortened_segment = fixed_segment_a.union(fixed_segment_b)
                        else:
                            print("failed all tests; examine {} manually.".format(k))
                #print(shortened_segment, shortened_segment.type, shortened_segment.length)
                if shortened_segment != "":
                    print("Updating line segment.")
                    row[2] = shortened_segment
                    cursor.updateRow(row)

    # select out just the short segments
    # TODO maybe don't include the short loops?
    output_short_segments = "ShortSegments_NeedtoAddRoutes"
    arcpy.Select_analysis(segmentfc, output_short_segments,
                          "{} < {}".format(arcpy.AddFieldDelimiters(segmentfc, "Shape_Length"),maxdistance-1))

    # select out just the short segments which also form loops -- will need special processing
    shortloops = [k for k in addroutes if k in issue1_segmentloops]
    if len(shortloops) > 0:
        output_loopedShortSegments = "ShortLoops"
        arcpy.Select_analysis(segmentfc, output_loopedShortSegments, """{} IN ({})""".format(arcpy.AddFieldDelimiters(segmentfc, "kmzid"), str(shortloops)[1:-1]))

    # Select out the completed (cut down) segments)
    output_trimmed_segments = os.path.join(nearFinalFD, "auditSegmentRoutes_TrimmedSegments")
    arcpy.Select_analysis(segmentfc_copy, output_trimmed_segments,
                          "{} >= {}".format(arcpy.AddFieldDelimiters(segmentfc, "Shape_Length"), maxdistance-1))

    # Start a list of feature classes to merge later.
    segroutestomerge = [output_trimmed_segments]

    ######################
    ### ROUTE ############
    ######################
    arcpy.MakeFeatureLayer_management(gcfcnameWGS84, "geocode_lyr")
    arcpy.SelectLayerByAttribute_management("geocode_lyr", "NEW_SELECTION", """{} IN ({})""".format(arcpy.AddFieldDelimiters("geocode_lyr", "kmzidint"), ", ".join(addroutes))) #TODO REVERT

    # Get route from participant to supermarket
    featureClassPrefix = "".join(os.path.basename(nets_csv).split("_")[:3]) #NETS_Supermarkets_yr2014
    routefc = getRouteToClosestDestination(nets_csv, "geocode_lyr", date_suf,
                                           featureClassPrefix)
    # copy routefc to new fc before modifying route
    routelimit_dist = maxdistance + 40
    routefc_copy = routefc + "_max{}m".format(routelimit_dist)
    arcpy.CopyFeatures_management(routefc, routefc_copy)
    # Cut routes down to maxdistance + 40m (leaves a little room just in case)
    with arcpy.da.UpdateCursor(routefc_copy, ["SHAPE@"]) as cursor:
        for row in cursor:
            feat = row[0].segmentAlongLine(0, routelimit_dist)
            row[0] = feat
            cursor.updateRow(row)
    print("Finished cutting down routes to {}m".format(routelimit_dist))
    #####################################
    ### COMBINED ROUTE 1 - SINGLEPART ###
    #####################################
    # combine routes (limited to maxdistance + 40) and segments
    mergedsegmentroute = "MergedRouteandSegment_cec_{}".format(date_suf)
    arcpy.Merge_management([output_short_segments, routefc_copy], mergedsegmentroute)
    # Dissolve by kmzid, multipart
    mergedsegmentroute_dslv = mergedsegmentroute + "_dslv"
    arcpy.Dissolve_management(mergedsegmentroute, mergedsegmentroute_dslv, "kmzid")
    # Add part count and start/end values to combined route and segment
    arcpy.AddGeometryAttributes_management(mergedsegmentroute_dslv,
                                           ["PART_COUNT", "LINE_START_MID_END"])
    # select out just the single-part route segments
    outputcombinedrte = "MergedRouteandSegment_mergedOK"
    arcpy.Select_analysis(mergedsegmentroute_dslv, outputcombinedrte,
                          "{} < 2".format(arcpy.AddFieldDelimiters(mergedsegmentroute_dslv, "PART_COUNT")))
    # Identify the multi-part route segments for later processing
    with arcpy.da.SearchCursor(
        mergedsegmentroute_dslv,
        ["kmzid", "PART_COUNT"], "{} > 1".format(arcpy.AddFieldDelimiters(mergedsegmentroute_dslv, "PART_COUNT"))) as cursor:
        issue2_multipartsegroutes = [row[0] for row in cursor]
    # select out the multipart segment routes for later processing
    outputcombinedrte_multipart = "MergedRouteandSegment_multipart"
    arcpy.Select_analysis(mergedsegmentroute_dslv, outputcombinedrte_multipart,
                          "{} > 1".format(arcpy.AddFieldDelimiters(mergedsegmentroute_dslv, "PART_COUNT")))

    ### CHECK DIRECTION OF SINGLEPART SUBSET ###
    # Add fields for later use
    addFlipFields(outputcombinedrte)
    # Create layer for selection
    arcpy.MakeFeatureLayer_management(outputcombinedrte, "segmentroutesSinglepart_lyr")
    # all should be singlepart in this feature class
    checklinedir("segmentroutesSinglepart_lyr", segmentendpoints, "flip1", 5)
    # Identify number of features to change and check outputs
    with arcpy.da.SearchCursor("segmentroutesSinglepart_lyr", ["kmzid", "flip1"], "{} = 1".format(arcpy.AddFieldDelimiters("segmentroutesSinglepart_lyr", "flip1"))) as cursor:
        needtoflippass1_sp = [row[0] for row in cursor]
    with arcpy.da.SearchCursor("segmentroutesSinglepart_lyr", ["kmzid", "flip1"], "{} = 0".format(arcpy.AddFieldDelimiters("segmentroutesSinglepart_lyr", "flip1"))) as cursor:
        goodpass1_sp = [row[0] for row in cursor]
    print("There are {} lines that need to be corrected in the singlepart merged route and segments".format(len(needtoflippass1_sp)))
    print("There are {} lines that go the correct direction in the singlepart merged route and segments".format(len(goodpass1_sp)))
    ### UPDATE DIRECTION ##
    arcpy.SelectLayerByAttribute_management(
       "segmentroutesSinglepart_lyr", "NEW_SELECTION", "{} = 1".format(arcpy.AddFieldDelimiters("segmentroutesSinglepart_lyr", "flip1")))
    arcpy.FlipLine_edit("segmentroutesSinglepart_lyr")
    ### CHECK DIRECTION AGAIN ##
    arcpy.SelectLayerByAttribute_management("segmentroutesSinglepart_lyr", "CLEAR_SELECTION")
    arcpy.AddGeometryAttributes_management("segmentroutesSinglepart_lyr", ["LINE_START_MID_END"])
    # Check whether start matches/is within a small distance of either segment endpoint
    checklinedir("segmentroutesSinglepart_lyr", segmentendpoints, "flip2", 5)
    # 3. If still not matching after flip:
    #    Set aside for later evaluation
    arcpy.SelectLayerByAttribute_management(
        "segmentroutesSinglepart_lyr", "NEW_SELECTION", "{} IN (1, 99)".format(arcpy.AddFieldDelimiters("segmentroutesSinglepart_lyr","flip2")))
    if int(arcpy.GetCount_management("segmentroutesSinglepart_lyr").getOutput(0)) > 0:
        print("There are still some lines that don't start near a segment endpoint.")
        arcpy.CopyFeatures_management("segmentroutesSinglepart_lyr", "problemSinglePartSegRoutes")
        with arcpy.da.SearchCursor("problemSinglePartSegRoutes", "kmzid") as cursor:
            problemSinglePartSegRoutes_list = [row[0] for row in cursor]
            print(problemSinglePartSegRoutes_list)
    else:
        print("All merged segment route lines (singlepart) start near a segment endpoint.")
    # Select the lines that are going in the proper direction and cut them down to maxdistance
    arcpy.SelectLayerByAttribute_management(
        "segmentroutesSinglepart_lyr",
        "NEW_SELECTION", "{} = 0".format(arcpy.AddFieldDelimiters("segmentroutesSinglepart_lyr","flip2")))
    print("A total of {} singlepart segment routes have the correct direction.".format(int(arcpy.GetCount_management("segmentroutesSinglepart_lyr").getOutput(0))))
    # copy to new fc just in case
    singlepartsegroutes_maxdist = "SinglePartSegmentRoutes_max{}m".format(maxdistance)
    arcpy.CopyFeatures_management("segmentroutesSinglepart_lyr", singlepartsegroutes_maxdist)
    print("Cutting SinglePart merged segments and routes to maximum distance of {}m.s".format(maxdistance))
    totalsinglepartsegroutes = int(arcpy.GetCount_management(singlepartsegroutes_maxdist).getOutput(0))
    with arcpy.da.UpdateCursor(singlepartsegroutes_maxdist, ["SHAPE@"]) as cursor:
        for row in cursor:
            feat = row[0].segmentAlongLine(0,maxdistance)
            row[0] = feat
            cursor.updateRow(row)
    if int(arcpy.GetCount_management(singlepartsegroutes_maxdist).getOutput(0)) == totalsinglepartsegroutes:
        print("Output from segmentalongline has same number of features after trimming to {}m".format(maxdistance))
    ### IDENTIFY SHORT ROUTES ###
    # All of the single part merges have completed google audit routes unless they are too short.
    issue3_shortsegroute = {}
    shortroutebreak = maxdistance - 1
    with arcpy.da.SearchCursor(singlepartsegroutes_maxdist, ["SHAPE@", "kmzid"]) as cursor:
        for row in cursor:
            featlength = row[0].length
            if featlength < shortroutebreak: # e.g. < 399m
                issue3_shortsegroute[row[1]] = row[0].length
    print("There are {} segment routes shorter than {}m in the singlepart processing step.".format(len(issue3_shortsegroute.keys()), shortroutebreak))
    # Select out just the segment routes which are long enough in the first pass
    singlepartsegroutes_maxdist_notshort = os.path.join(
        nearFinalFD, "auditSegmentRoutes_" + singlepartsegroutes_maxdist + "_lengthOK")
    arcpy.Select_analysis(singlepartsegroutes_maxdist, singlepartsegroutes_maxdist_notshort,
                          "{} >= {}".format(arcpy.AddFieldDelimiters(singlepartsegroutes_maxdist_notshort, "Shape_Length"),shortroutebreak))
    # Select the too short segment routes for later processing
    singlepartsegroutes_maxdist_short = singlepartsegroutes_maxdist + "_shortAddSecondDest"
    arcpy.Select_analysis(singlepartsegroutes_maxdist, singlepartsegroutes_maxdist_short,
                          "{} < {}".format(arcpy.AddFieldDelimiters(singlepartsegroutes_maxdist_short, "Shape_Length"),shortroutebreak))
    # Add the single part segment routes that are long enough to the list of feature classes to merge later
    segroutestomerge.append(singlepartsegroutes_maxdist_notshort)

    #####################################
    ### COMBINED ROUTE 2 - MULTIPART ###
    #####################################
    ### Another reason for multipart segment routes: Segment is not a loop, but slight offset in segment/route
    # Use the list from earlier in processing
    multipart_segroute_nonloop = [x for x in issue2_multipartsegroutes if x not in issue1_segmentloops]
    multipart_segroute_nonloop.sort()
    #TODO add if statement to skip processing
    # if len(multipart_segroute_nonloop) > 0:
    # Copy input feature classes, subset to just the ones we are integrating
    for inputfc in [segmentfc, routefc_copy]:
        arcpy.Select_analysis(inputfc, inputfc + "_integrate", """{} in ({})""".format(arcpy.AddFieldDelimiters(inputfc, "kmzid"),str(multipart_segroute_nonloop)[1:-1]))
    integrated_segments = segmentfc + "_integrate"
    integrated_routes = routefc_copy + "_integrate"
    arcpy.MakeFeatureLayer_management(integrated_segments, "segmentsInt_lyr")
    arcpy.MakeFeatureLayer_management(integrated_routes, "routeInt_lyr")
    # loop through all of the kmzids in the multipart list (except loops)
    # Could probably process all at once but would prefer to handle individually
    # to ensure that some routes do not "snap to" neighboring route.
    progresscount = 0
    # keep track of temporary layers for later merge
    mergelist = []
    for k in multipart_segroute_nonloop:
        whereclause = "{} = '{}'".format(arcpy.AddFieldDelimiters("routeInt_lyr", "kmzid"),k)
        arcpy.SelectLayerByAttribute_management("routeInt_lyr", "NEW_SELECTION", whereclause)
        whereclausesegments = "{} = '{}'".format(arcpy.AddFieldDelimiters("segmentsInt_lyr", "kmzid"),k)
        arcpy.SelectLayerByAttribute_management("segmentsInt_lyr", "NEW_SELECTION", whereclausesegments)
        # Integrate the layers if within 5 meters of each other, using Segment as the primary layer
        arcpy.Integrate_management([["segmentsInt_lyr", 1],["routeInt_lyr",2]], "5 Meters")
        print(arcpy.GetMessages())
        # store in_memory for now
        segmenttmp = r"in_memory\segments_{}".format(k)
        routetmp = r"in_memory\routes_{}".format(k)
        arcpy.CopyFeatures_management("segmentsInt_lyr", segmenttmp)
        arcpy.CopyFeatures_management("routeInt_lyr", routetmp)
        mergelist.extend([segmenttmp, routetmp])
        progresscount += 1
        print("Finished integrating {} of {}".format(progresscount, len(multipart_segroute_nonloop)))
    outmerge = "MergedRouteandSegmentIntegratedSubset_cec_{}".format(date_suf)
    arcpy.Merge_management(mergelist, outmerge)
    outdslv = outmerge + "_dslv"
    arcpy.Dissolve_management(outmerge, outdslv, "kmzid")
    arcpy.AddGeometryAttributes_management(outdslv, ["PART_COUNT", "LINE_START_MID_END"])
    # What's still multipart?
    with arcpy.da.SearchCursor(
        outdslv,
        ["kmzid", "PART_COUNT"], "{} > 1".format(arcpy.AddFieldDelimiters(outdslv, "PART_COUNT"))) as cursor:
        issue4_multipartAfterIntegrate = [row[0] for row in cursor]
    # Extract single part features after integrate/dissolve and prepare to cut down
    multipartsegroutes_maxdist = os.path.join(
        nearFinalFD, "auditSegmentRoutes_MultiPartSegmentRoutesIntegratedOK_max{}m".format(maxdistance))
    arcpy.Select_analysis(
        outdslv,
        multipartsegroutes_maxdist,
        "{} < 2".format(arcpy.AddFieldDelimiters(outdslv, "PART_COUNT")))
    ### CHECK DIRECTION OF INTEGRATED SUBSET ##
    # linestartend already added in previous step.
    addFlipFields(multipartsegroutes_maxdist)
    arcpy.MakeFeatureLayer_management(multipartsegroutes_maxdist, "integrated_segroutes")
    checklinedir("integrated_segroutes", segmentendpoints, "flip1", 5)
    with arcpy.da.SearchCursor("integrated_segroutes", ["kmzid", "flip1"], "{} = 1".format(arcpy.AddFieldDelimiters("segmentroutesSinglepart_lyr", "flip1"))) as cursor:
        needtoflippass1_mp = [row[0] for row in cursor]
    with arcpy.da.SearchCursor("integrated_segroutes", ["kmzid", "flip1"], "{} = 0".format(arcpy.AddFieldDelimiters("segmentroutesSinglepart_lyr", "flip1"))) as cursor:
        goodpass1_mp = [row[0] for row in cursor]
    print("There are {} lines that need to be corrected in the integrated multipart merged route and segments".format(len(needtoflippass1_mp)))
    print("There are {} lines that go the correct direction in the integrated multipart merged route and segments".format(len(goodpass1_mp)))
    ### UPDATE DIRECTION IF NEEDED ##
    arcpy.SelectLayerByAttribute_management(
       "integrated_segroutes", "NEW_SELECTION", "{} = 1".format(arcpy.AddFieldDelimiters("integrated_segroutes", "flip1")))
    arcpy.FlipLine_edit("integrated_segroutes")
    ### CHECK DIRECTION AGAIN ##
    arcpy.SelectLayerByAttribute_management("integrated_segroutes", "CLEAR_SELECTION")
    arcpy.AddGeometryAttributes_management("integrated_segroutes", ["LINE_START_MID_END"])
    checklinedir("integrated_segroutes", segmentendpoints, "flip2", 5)
    # 3. If still not matching after flip:
    #    Set aside for later evaluation
    arcpy.SelectLayerByAttribute_management(
        "integrated_segroutes", "NEW_SELECTION", "{} IN (1, 99)".format(arcpy.AddFieldDelimiters("integrated_segroutes","flip2")))
    with arcpy.da.SearchCursor("integrated_segroutes", "kmzid") as cursor:
        problemMultiPartSegRoutes_list = [row[0] for row in cursor]
    if int(arcpy.GetCount_management("integrated_segroutes").getOutput(0)) > 0:
        print("There are still some (previously multipart) lines that don't start near a segment endpoint.")
        arcpy.CopyFeatures_management("integrated_segroutes", "problemMultiPartSegRoutes")
    else:
        print("All merged segment route lines (integrated) start near a segment endpoint.")
    # Cut down route segments
    arcpy.SelectLayerByAttribute_management(
       "integrated_segroutes", "NEW_SELECTION", "{} = 0".format(arcpy.AddFieldDelimiters("integrated_segroutes", "flip2")))
    with arcpy.da.UpdateCursor("integrated_segroutes", ["SHAPE@", "kmzid"]) as cursor:
        for row in cursor:
            if row[0].length < shortroutebreak: # within e.g. 399m rather than 400m
                issue3_shortsegroute[row[1]] = row[0].length
            feat = row[0].segmentAlongLine(0,maxdistance)
            row[0] = feat
            cursor.updateRow(row)
    print("finished cutting down previously multipart segment routes which merged OK")
    segroutestomerge.append(multipartsegroutes_maxdist)

    # handle the short routes
    ####################################
    ### SHORT - NEXT SMA Destination ###
    ####################################
    shortkmzids = issue3_shortsegroute.keys()
    # TODO skip processing if there aren't any results
    print("Analyzing second destination for short (singlepart) routes ({} total)...".format(len(shortkmzids)))
    firstdest = dict()
    with arcpy.da.SearchCursor(routefc, ["kmzid", "Name", "uhcid"]) as cursor:
        for row in cursor:
            if row[0] in shortkmzids:
                # print(row[0], row[1], row[2])
                firstdest[row[0]] = row[2]

    # Which, if any, items can be left as is? (process_rca == 0)
    no_additional_processing_noDest2 = []
    # shortroutes
    with arcpy.da.SearchCursor(gcCEC, ["kmzid", "process_rca"]) as cursor:
        for row in cursor:
            if str(row[0]) in firstdest.keys():
                if row[1] == 0:
                    print(row[0])
                    no_additional_processing.append(str(row[0]))
    # print the lengths of items that are too short but don't need a second destination
    if len(no_additional_processing_noDest2) > 0:
        for i in no_additional_processing_noDest2:
            print(i, firstdest[i])
    else:
        print("All short routes need a second destination if possible.")
    print("There are {} short segment routes (need a second destination)".format(len(issue3_shortsegroute) - len(no_additional_processing_noDest2)))

    # make closest facility layer and create route to next supermarket (closest to first supermarket)
    # First need to delete the in-memory layers that were generated for the previous Closest Facility Analysis Layer
    arcpy.Delete_management("Incidents")
    arcpy.Delete_management("Facilities")
    arcpy.Delete_management("Routes")

    # Create a layer to use to limit the supermarkets
    businessfcnameWGS84 = featureClassPrefix + "_" + date_suf + "_WGS84"
    arcpy.MakeFeatureLayer_management(businessfcnameWGS84, "businesses_lyr")
    # Create a Closest Facility Analysis layer
    result_object_cf2 = arcpy.na.MakeClosestFacilityAnalysisLayer(network,
        "Closest_{}_seconddest".format(featureClassPrefix), travel_mode, line_shape="ALONG_NETWORK")
    #Get the layer object from the result object. The closest facility layer can
    #now be referenced using the layer object.
    layer_object_cf2 = result_object_cf2.getOutput(0)
    #Get the names of all the sublayers within the closest facility layer.
    sublayer_names_cf2 = arcpy.na.GetNAClassNames(layer_object_cf2)
    #Stores the layer names that we will use later

    facilities_layer_name_cf2 = sublayer_names_cf2["Facilities"]
    incidents_layer_name_cf2 = sublayer_names_cf2["Incidents"]
    routes_layer_name_cf2 = sublayer_names_cf2["CFRoutes"]
    #TODO maybe move these next 3 lines
    facilities_sublayer = layer_object_cf2.listLayers(facilities_layer_name_cf2)[0]
    incidents_sublayer = layer_object_cf2.listLayers(incidents_layer_name_cf2)[0]
    routes_sublayer = layer_object_cf2.listLayers(routes_layer_name_cf2)[0]

    # Map the Name property to the uhcid field
    print("setting up field mappings")
    arcpy.na.AddFieldToAnalysisLayer(layer_object_cf2, facilities_layer_name_cf2,
                                     "uhcid", "DOUBLE")
    field_mappings_cf2_potentialdest = arcpy.na.NAClassFieldMappings(
       layer_object_cf2, facilities_layer_name_cf2)
    field_mappings_cf2_potentialdest["Name"].mappedFieldName = "uhcid"
    field_mappings_cf2_potentialdest["uhcid"].mappedFieldName = "uhcid"

    # Incidents -- starting store location
    arcpy.na.AddFieldToAnalysisLayer(layer_object_cf2, sublayer_names_cf2["Incidents"],
                                     "uhcid", "DOUBLE")
    arcpy.na.AddFieldToAnalysisLayer(layer_object_cf2, sublayer_names_cf2["Incidents"],
                                     "kmzid", "TEXT", field_length=10)

    field_mappings_cf2_start = arcpy.na.NAClassFieldMappings(
        layer_object_cf2, sublayer_names_cf2["Incidents"])
    field_mappings_cf2_start["uhcid"].mappedFieldName = "uhcid"
    fieldMapName = field_mappings_cf2_start["Name"]
    fieldMapK = field_mappings_cf2_start["kmzid"]

    arcpy.na.AddFieldToAnalysisLayer(layer_object_cf2, sublayer_names_cf2["CFRoutes"],
                                     "kmzid", "TEXT", field_length=10)
    arcpy.na.AddFieldToAnalysisLayer(layer_object_cf2, sublayer_names_cf2["CFRoutes"],
                                     "uhcid", "TEXT", field_length=20)
    # Load all of the businesses to start
    print("Loading Locations - All supermarkets")
    arcpy.na.AddLocations(layer_object_cf2, facilities_layer_name_cf2,
                              "businesses_lyr", field_mappings_cf2_potentialdest, "", append="CLEAR")
    # Load the problem supermarket locations one at a time; remove matching supermarket locations before solving

    seconddestfcs = []
    seconddestfcs_cec = []
    firstroutes = list(firstdest.items())
    firstroutes.sort()
    for i, rte in enumerate(firstroutes):
        k = rte[0]
        v = rte[1]
        print("Processing item {}: {}".format(i, k))
        outputroute = "secondDestination_{}".format(k)
        outputroute_cec = "secondDestination_{}_cec".format(k)
        # select just the single supermarket location by uhcid
        arcpy.SelectLayerByAttribute_management(
            "businesses_lyr", "NEW_SELECTION",
            """{} = {}""".format(arcpy.AddFieldDelimiters("businesses_lyr", "uhcid"),v))
        fieldMapName.defaultValue = k
        fieldMapK.defaultValue = k
        # default travel from incident to facilities, no cutoff. Use default number of facilities to find (1)
        # Load the selected uhcid as the starting location
        arcpy.na.AddLocations(layer_object_cf2, sublayer_names_cf2["Incidents"],
                              "businesses_lyr", field_mappings_cf2_start, "", append="CLEAR")
        # remove matching locations from the Facilities sublayer
        arcpy.SelectLayerByLocation_management(facilities_sublayer, "ARE_IDENTICAL_TO", incidents_sublayer)
        if int(arcpy.GetCount_management(facilities_sublayer).getOutput(0)) > 0:
            arcpy.DeleteFeatures_management(facilities_sublayer)
        arcpy.SelectLayerByAttribute_management(facilities_sublayer, "CLEAR_SELECTION")
        # Solve
        arcpy.na.Solve(layer_object_cf2)
        print(arcpy.GetMessages())
        # copy output
        print("Copying Routes to feature class")
        arcpy.CopyFeatures_management(routes_sublayer, outputroute)
        print(arcpy.GetMessages())
        # Add the kmzid field to the route
        arcpy.CalculateField_management(outputroute, "kmzid", "!Name!.split(' - ')[0]", "PYTHON3")
        arcpy.CalculateField_management(outputroute, "uhcid", "!Name!.split(' - ')[1]")
        seconddestfcs.append(outputroute)
        print("Reprojecting routes")
        # Reproject to USCEC
        arcpy.Project_management(outputroute, outputroute_cec, srCEC)
        print(arcpy.GetMessages())
        seconddestfcs_cec.append(outputroute_cec)
        # Add the incident back to potential destinations (facilities sublayer)
        print("adding Incident back to potential destinations (Facilities)")
        arcpy.na.AddLocations(layer_object_cf2, facilities_layer_name_cf2,
                              incidents_sublayer)
    print("{} new routes from first supermarket to second supermarket created -- REVIEW".format(len(firstdest)))

    # Try to link the initial segment routes to the route from first destination to second destination
    # If route2 plus segment & route1 creates a polygon, set aside for further evaluation.
    tworoutes_list = []
    # set environment variable to be used within dissolve
    arcpy.env.XYTolerance = "0.5 Meters"
    for secondroute in seconddestfcs_cec:
        # cut down overall second route to 2.5 times the maximum distance (1000m)
        with arcpy.da.UpdateCursor(secondroute, "SHAPE@") as ucursor:
            for row in ucursor:
                row[0] = row[0].segmentAlongLine(0, maxdistance*2.5)
                ucursor.updateRow(row)
        kmzid = secondroute.split("_")[1] # e.g. secondDestination_<kmzid>_cec
        tworouteoutput = "auditRouteTwoDests_{}".format(kmzid)
        merge_tmp = r"in_memory\mergetmp_{}".format(kmzid)
        firstroute_tmp = r"in_memory\firstrouteseg_{}".format(kmzid)
        mergeAfterRte2Cut_tmp = r"in_memory\mergetmpAfterCut_{}".format(kmzid)
        arcpy.Select_analysis(singlepartsegroutes_maxdist_short, firstroute_tmp, """{} = '{}'""".format(arcpy.AddFieldDelimiters(singlepartsegroutes_maxdist_short, "kmzid"),kmzid))
        arcpy.Merge_management([secondroute, firstroute_tmp], merge_tmp)
        polygontest_tmp = r"in_memory\polygontest_{}".format(kmzid)
        arcpy.FeatureToPolygon_management(merge_tmp, polygontest_tmp)
        numpolys = int(arcpy.GetCount_management(polygontest_tmp).getOutput(0))
        if numpolys == 0:
            print("route2 can be shortened normally by comparison with segmentRoute1: {}".format(kmzid))
            cut_second_route(firstroute_tmp, secondroute, kmzid, maxdistance)
            # merge and dissolve the results
            arcpy.Merge_management([firstroute_tmp, secondroute], mergeAfterRte2Cut_tmp)
            arcpy.Dissolve_management(mergeAfterRte2Cut_tmp, tworouteoutput, "kmzid")
            tworoutes_list.append(tworouteoutput)
        else: # merge created polygon(s) and test
            print("Dissolving {}...".format(kmzid))
            dslv_tmp = r"in_memory\tworoutes_dslv_{}".format(kmzid)
            arcpy.Dissolve_management(merge_tmp, dslv_tmp, "kmzid")
            print(arcpy.GetMessages())
            print("Adding part count")
            arcpy.AddGeometryAttributes_management(dslv_tmp, "PART_COUNT")
            print(arcpy.GetMessages())
            print("Getting report of part count")
            numparts = int([row[0] for row in arcpy.da.SearchCursor(dslv_tmp, "PART_COUNT")][0])
            if numparts == 1:
                print("Route2 needs to be dissolved then evaluated for length: {}".format(kmzid))
                with arcpy.da.UpdateCursor(dslv_tmp, "SHAPE@") as ucursor2:
                    for row in ucursor2:
                        row[0] = row[0].segmentAlongLine(0,maxdistance)
                        cursor.updateRow(row)
                # copy the in memory dataset
                arcpy.CopyFeatures_management(dslv_tmp, tworouteoutput)
                tworoutes_list.append(tworouteoutput)
            else: #still multipart after dissolve with high tolerance
                print("Route2 plus SegmentRoute1 causes a loop. Process separately. kmizd: {}".format(kmzid))
                adjustAuditRouteManually.append(kmzid)
                multipart_dissolve = "manualReview_tmpDslvMultipart_{}".format(kmzid)
                arcpy.MultipartToSinglepart_management(dslv_tmp, multipart_dissolve)
                print('Created multipart output for {} to review manually'.format(multipart_dissolve))
                # TODO experiment to see if we can capture the correct segments or manually adjust this audit route.
        # Note: when using ArcGIS Pro, can delete a list. In 10.5 cannot have in list.
        arcpy.Delete_management(firstroute_tmp)
        arcpy.Delete_management(merge_tmp)
        arcpy.Delete_management(polygontest_tmp)
        arcpy.Delete_management(mergeAfterRte2Cut_tmp)
    output_segmentPlusTwoRoutes = os.path.join(
        nearFinalFD, "auditSegmentRoutes_short_secondDestAdded_{}".format(len(tworoutes_list)))
    arcpy.Merge_management(tworoutes_list, output_segmentPlusTwoRoutes)
    segroutestomerge.append(output_segmentPlusTwoRoutes)

    # multipart audit routes
    print("Checking for multipart audit routes that can be left as-is...")
    no_additional_processing_multipart = []
    no_additional_processing_shortloops = []
    with arcpy.da.SearchCursor(gcCEC, ["kmzid", "process_rca"]) as cursor:
        for row in cursor:
            if str(row[0]) in issue4_multipartAfterIntegrate:
                if row[1] == 0:
                    print(row[0])
                    no_additional_processing_multipart.append(str(row[0]))
            if str(row[0]) in shortloops:
                if row[1] == 0:
                    print(row[0])
                    no_additional_processing_shortloops.append(str(row[0]))
    print("There are {} ids which are OK to leave as-is (segment only; would be multipart if route added): {}".format(len(no_additional_processing_multipart), no_additional_processing_multipart))
    # Extract segments (only) that don't need additional steps
    segments_noRCA = os.path.join(
        nearFinalFD, "auditSegmentRoutes_Segment_noAdditionalProcessingNoRCA")
    arcpy.Select_analysis(segmentfc, segments_noRCA,
                          """{} IN ({})""".format(arcpy.AddFieldDelimiters(segmentfc, "kmzid"),str(no_additional_processing_multipart)[1:-1]))
    segroutestomerge.append(segments_noRCA)

    #############################
    ### LOOPED SHORT SEGMENTS ###
    #############################
    # These are loops that were smaller than the maximum distance; the larger loops were already processed
    # Handle the route segments for segments that are a loop (cul de sac, development, etc.)
    # Take the total length of the segment loop and subtract it from the maximum distance (newroutelength).
    # Extract the part of the route that doesn't overlap the segment and cut it to the new routelength.
    # TODO: Q - Is this still confusing for Google Earth Audit purposes?
    # get the geometry of the route which is unique to the route (not part of the looped segment)

    print("There are {} short loops which need to be extended.".format(len(shortloops)))
    routefc_forloops = "routes_for_shortloops"
    arcpy.Select_analysis(routefc_copy, routefc_forloops, """{} IN ({})""".format(arcpy.AddFieldDelimiters(segmentfc, "kmzid"), str(shortloops)[1:-1]))
    with arcpy.da.UpdateCursor(routefc_forloops, ["kmzid", "SHAPE@"]) as routeloopcursor:
        for rlrow in routeloopcursor:
            k_id = rlrow[0]
            print(k_id, type(k_id))
            routefeat = rlrow[1]
            loopfeat = [lrow[0] for lrow in arcpy.da.SearchCursor(output_loopedShortSegments, ["SHAPE@"], """{} = '{}'""".format(arcpy.AddFieldDelimiters(output_loopedShortSegments, "kmzid"), k_id))][0]
            print("Loop Length: {}".format(loopfeat.length))
            routefeat2 = routefeat.difference(loopfeat)
            print("Starting Route Length: {}".format(routefeat.length))
            print("Route Length not counting overlap with segment loop: {}".format(routefeat2.length))
            rlrow[1] = routefeat2
            routeloopcursor.updateRow(rlrow)
    # Figure out whether the new route starts at or ends at the looped segment -- manual review: these are correct.
    # update direction if needed -- manual review: not needed.
    # cut the route down to maximum maxdist - segmentlength
    with arcpy.da.UpdateCursor(routefc_forloops, ["kmzid", "SHAPE@"]) as routeloopcursor:
        for rlrow in routeloopcursor:
            k_id = rlrow[0]
            routefeat = rlrow[1]
            loopfeat = [lrow[0] for lrow in arcpy.da.SearchCursor(output_loopedShortSegments, ["SHAPE@"], """{} = '{}'""".format(arcpy.AddFieldDelimiters(output_loopedShortSegments, "kmzid"), k_id))][0]
            print("Loop Length: {}".format(loopfeat.length))
            print("Starting Route Length: {}".format(routefeat.length))
            routeonlylengthneeded = maxdistance - loopfeat.length # e.g. 400 - 86 ~= 313
            print("new route length (trimmed): {}".format(routeonlylengthneeded))
            routefeat2 = routefeat.segmentAlongLine(0, routeonlylengthneeded)
            rlrow[1] = routefeat2
            routeloopcursor.updateRow(rlrow)
    # combine the segment loops and routes
    # these will be multipart, but that's OK in this case.
    arcpy.Merge_management([routefc_forloops, output_loopedShortSegments], r"in_memory\loopedsegments_routes_merged")
    output_loopedShortSegmentsPlusRoutes = os.path.join(nearFinalFD,
                                                        "auditSegmentRoutes_LoopedShortSegmentsPlusRoutes")
    arcpy.Dissolve_management(r"in_memory\loopedsegments_routes_merged", output_loopedShortSegmentsPlusRoutes,"kmzid")
    segroutestomerge.append(output_loopedShortSegmentsPlusRoutes)

    ####################################
    ### Multipart to Audit           ###
    ####################################
    # Similar to segment loop section above, what if we select out the route portion that doesn't overlap the segment?
    # The output may or may not be multipart
    print("{} segments and routes form a loop or multipart feature when combined.".format(len(issue4_multipartAfterIntegrate) - len(no_additional_processing_multipart)))
    routefc_formultipartsubset = "routes_for_multipartsubset"
    mpsubset = [x for x in issue4_multipartAfterIntegrate if x not in no_additional_processing_multipart]
    arcpy.Select_analysis(routefc_copy, routefc_formultipartsubset,
                          """{} IN ({})""".format(arcpy.AddFieldDelimiters(segmentfc, "kmzid"), str(mpsubset)[1:-1]))

    # Confirm these would be multipart because they would form polygons; if not, need to address manually
    arcpy.MakeFeatureLayer_management(outputcombinedrte_multipart, "mp_segroute_lyr")

    postMPprocessing_list = list()
    with arcpy.da.UpdateCursor(routefc_formultipartsubset, ["kmzid", "SHAPE@"]) as cursor:
        for row in cursor:
            k_id = row[0]
            print(k_id)
            # Select the original combined route and segment
            arcpy.SelectLayerByAttribute_management(
                "mp_segroute_lyr", "NEW_SELECTION",
                """{} = '{}'""".format(arcpy.AddFieldDelimiters("mp_segroute_lyr", "kmzid"), k_id))
            polygontest_tmp = r"in_memory\polygontest_{}".format(k_id)
            arcpy.FeatureToPolygon_management("mp_segroute_lyr", polygontest_tmp)
            numpolys = int(arcpy.GetCount_management(polygontest_tmp).getOutput(0))
            if numpolys >= 1:
                postMPprocessing_list.append(k_id)
                routefeat = row[1]
                segfeat = [lrow[0] for lrow in arcpy.da.SearchCursor(segmentfc, ["SHAPE@"], """{} = '{}'""".format(arcpy.AddFieldDelimiters(segmentfc, "kmzid"), k_id))][0]
                print("Segment Length: {}".format(segfeat.length))
                routefeat2 = routefeat.difference(segfeat)
                print("Starting Route Length: {}".format(routefeat.length))
                print("Route Length not counting overlap with segment loop: {}".format(routefeat2.length))
                routeonlylengthneeded = maxdistance - segfeat.length # e.g. 400 - 86 ~= 313
                row[1] = routefeat2.segmentAlongLine(0, routeonlylengthneeded)
            else:
                print("The segment route combination for {} is multipart but does not form a polygon.".format(k_id))
                print("Check this manually.")
                adjustAuditRouteManually.append(k_id)
            cursor.updateRow(row)

    # select the routes that were just shortened (leave out the segroutes that were multipart for a different reason
    # than looping)
    arcpy.MakeFeatureLayer_management(routefc_formultipartsubset, "mp_route_lyr")
    arcpy.SelectLayerByAttribute_management(
        "mp_route_lyr", "NEW_SELECTION",
        """{} IN ({})""".format(arcpy.AddFieldDelimiters("mp_route_lyr", "kmzid"), str(postMPprocessing_list)[1:-1]))

    # select the associated segments
    arcpy.MakeFeatureLayer_management(segmentfc, "segment_lyr")
    arcpy.SelectLayerByAttribute_management(
        "segment_lyr", "NEW_SELECTION",
        """{} IN ({})""".format(arcpy.AddFieldDelimiters("segment_lyr", "kmzid"), str(postMPprocessing_list)[1:-1]))

    # merge these and dissolve
    mergeMP_tmp = r"in_memory\mergedmultipartsubset"
    arcpy.Merge_management(["mp_route_lyr", "segment_lyr"], mergeMP_tmp)
    output_mp_segroute_max = os.path.join(nearFinalFD, "auditSegmentRoutes_SegmentRoutesClippedBasedOnRouteOverlap")
    arcpy.Dissolve_management(mergeMP_tmp, output_mp_segroute_max, "kmzid")
    segroutestomerge.append(output_mp_segroute_max)

    # exit main function, returning the name of the feature class which has the segments, maximum distance (400m)
    # This will be used in the creation of the "Start and end" segment indicators as well as the midpoint indicator for segment
    print("All current Google Audit processing complete")
    print("Complete manual processing for {}".format(adjustAuditRouteManually))
    print("Then create segment midpoint indicators, segment start/endpoint indicators, and side of street indicator.")
    return segmentfc_copy, segroutestomerge, adjustAuditRouteManually, gcCEC

dateprocessed = "20190628"
participant_segments, completed_audit_route_fcs, idsForManualFixes, geocodesCEC = main(dateprocessed)
fcsprocessed = [int(arcpy.GetCount_management(fc).getOutput(0)) for fc in completed_audit_route_fcs]
output_semiFinal = "GoogleAuditSegRoutes_NearFinal_{}".format(sum(fcsprocessed))
arcpy.Merge_management(completed_audit_route_fcs, output_semiFinal)

#TODO incorporate into overall process
##############################################
### MIDPOINT INDICATORS FOR SIDE OF STREET ###
##############################################
# https://gis.stackexchange.com/questions/189902/determining-whether-point-is-on-left-or-right-side-of-road-using-arcpy
createSideOfStreetIndicators(geocodesCEC, participant_segments)


if __name__ == "__main__":
    dateprocessed = "20190628"
    main(dateprocessed)


# Check extension back in
arcpy.CheckInExtension("network")








