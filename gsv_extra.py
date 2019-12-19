"""
Created by: Stephanie Miller
Last revision: 2019-12-18 temporary remove of extra code; may turn into tests
"""

##################################
# #TODO START HERE and incorporate again into MAIN if possible
#
# # TODO - Process Manually or figure out how to script this part...
# # One last problem segment / route combination -- will need to be adjusted manually.
# routefc_formultipartsubset2 = "routes_for_multipartsubset2"
# mpsubset2 = ["6601484"]
# arcpy.Select_analysis(routefc_copy, routefc_formultipartsubset2, """{} IN ({})""".format(arcpy.AddFieldDelimiters(segmentfc, "kmzid"), str(mpsubset2)[1:-1]))
# # Manually select the segments from with service area to form a route
# manualaudit = "ManualAuditRoute_6601484"
# manualaudit_dslv = manualaudit + "_dslv"
# arcpy.Dissolve_management(manualaudit, manualaudit_dslv, "kmzid")
# # Confirmed direction after dissolve manually in ArcMap
# # cut down route
# with arcpy.da.UpdateCursor(manualaudit_dslv, ["SHAPE@"]) as cursor:
#     for row in cursor:
#         feat = row[0].segmentAlongLine(0, maxdistance)
#         row[0] = feat
#         cursor.updateRow(row)
# segroutestomerge.append(manualaudit_dslv)
#
#
#
#
# #####################################
# ### START AND ENDPOINT INDICATORS ###
# #####################################
# # TODO create start and points (Feature Vertices to Points?)
# ##########################################
# ### Midpoint side of street INDICATORS ###
# ##########################################
# # TODO look back at old code; general idea:
# # get midpoint of SEGMENT (FVTP?)
# # Add geometry indicator for line bearing of SEGMENT
# # using line bearing (0 = north) and geocodes_cec, determine
# # new position for midpoint (offset from line ~ x m in perpendicular direction)
# ####################
# ### KMZ creation ###
# ####################
# # TODO join back to original geocode layer by kmzid
# # This will add back on the "auditgr" column.
# # Either a) joinfield and add auditgr, b) addjoin and calculate the
# # group or c) (probably easiest) add join and export to a new fc
# # Create Map with routesegments, segments, endpoints, side of street.
# # Check / update google audit documentation for symbols to use for each.
# # Select each unique audtitgroup, save mxd, and export to kmz.
# # max need to experiment with this or create a tool to process.
# # Once kmz is created, using location that will not be audited:
# # --get screenprints of segments on (and endpoints off)
# # --get screenprints of endpoints on (and segments off)
# # Ask for feedback


>>> with arcpy.da.SearchCursor("GoogleAuditSegRoutes_Final_2825", ["SHAPE@", "kmzid"]) as routecursor:
...     for row in routecursor:
...         k = row[1]
...         distancedict[k] = {}
...         routefeat = row[0]
...         with arcpy.da.SearchCursor("Segment_Endpoints", ["SHAPE@", "kmzid", "OBJECTID"], "{} = '{}'".format(arcpy.AddFieldDelimiters("Segment_Endpoints", "kmzid"), k)) as segcursor:
...             for srow in segcursor:
...                 objid = srow[2]
...                 print(k, objid)
...                 segfeat = srow[0]
...                 meas = routefeat.measureOnLine(segfeat)
...                 distancedict[k][objid] = meas


pointdict = {}
>>> for k in distancedict.keys():
...     min_key = min(distancedict[k], key=distancedict[k].get)
...     pointdict[min_key] = 1
...     max_key = max(distancedict[k],key=distancedict[k].get)
...     pointdict[max_key] = 2
...
>>> len(pointdict)


>>> pointdictmax = {}
>>> for k in distancedict.keys():
...     max_key = max(distancedict[k],key=distancedict[k].get)
...     pointdictmax[max_key] = 2
...
>>> len(pointdictmax)
2825
>>> pointdictmin = {}
>>> for k in distancedict.keys():
...     min_key = min(distancedict[k], key=distancedict[k].get)
...     pointdictmin[min_key] = 1

# some keys (objectIds) are the same for min and max (e.g. loops).
>>> dupkeys = [x for x in pointdictmin.keys() if x in pointdictmax.keys()]
>>> dupkeys
[664, 883, 2333, 2729, 3803, 4521, 5616]
>>> len(dupkeys)
7

with arcpy.da.UpdateCursor("segment_endpoints", ["OJBJECTID", "PointNum"]) as ucursor:
    loop_processed = []
    probably_loop = []
    for urow in ucursor:
        objid = urow[0]
        if objid not in pointdict.keys[1]:
            probably_loop.append(objid)
            urow[1] = 2
        elif objid in dupkeys:
            loop_processed.append(objid)
            urow[1] = 1
        else:
            urow[1] = pointdict[objid]
        ucursor.updateRow(urow)



def getfieldnames(fc):
    fieldlist = [field.name for field in arcpy.ListFields(fc)]
    return fieldlist

redcapdict = dict()
arcpy.env.workspace = r"M:\EPL-GEO_BSERE\Data\WorkingData\google.audit\intermediate-output\googleaudit_20190628.gdb\Final_audit_layers"

fcs = arcpy.ListFeatureClasses("*_sorted")

with arcpy.da.SearchCursor("streetside_indicators_2828_sorted2", ["kmzid", "auditgr", "process_rca"]) as cursor:
    for row in cursor:
        k = row[0]
        auditgr = row[1]
        rcaAudit = row[2]
        redcapdict[k] = dict()
        redcapdict[k]["auditgr"] = auditgr
        redcapdict[k]["process_rca"] = rcaAudit

with arcpy.da.SearchCursor("segment_endpoints_2828_sorted2", ["kmzid", "auditgr", "process_rca"]) as cursor:
    # verify auditgroup and kmzid
    for row in cursor:
        k = int(row[0])
        auditgr = row[1]
        rcaAudit = row[2]
        # check auditgroup
        auditgr_test = redcapdict[k]["auditgr"]
        rcaAudit_test = redcapdict[k]["process_rca"]
        if (auditgr_test != auditgr) or (rcaAudit_test != rcaAudit):
            print(k, rcaAudit_test == rcaAudit, auditgr_test == auditgr)


with arcpy.da.SearchCursor("street_segments_2828_sorted2", ["kmzid", "auditgr", "process_rca", "oddRoute", "loopedSeg", "twoDests", "Shape_Length"]) as cursor:
    for row in cursor:
        k = int(row[0])
        auditgr = row[1]
        rcaAudit = row[2]
        # check auditgroup
        auditgr_test = redcapdict[k]["auditgr"]
        rcaAudit_test = redcapdict[k]["process_rca"]
        # verify auditgroup and kmzid
        if (auditgr_test != auditgr) or (rcaAudit_test != rcaAudit):
            print(k, rcaAudit_test == rcaAudit, auditgr_test == auditgr)
            # TODO raise an error
        oddRouteFlag = row[3]
        loopFlag = row[4]
        twoDestsFlag = row[5]
        segLength = round(row[6], 2)
        redcapdict[k]["oddRoute"] = oddRouteFlag
        redcapdict[k]["loopedSeg"] = loopFlag
        redcapdict[k]["twoDests"] = twoDestsFlag
        redcapdict[k]["Segment_Length"] = segLength


with arcpy.da.SearchCursor("audit_routes_2828_sorted2", ["kmzid", "auditgr", "process_rca", "oddRoute", "loopedSeg", "twoDests", "Shape_Length"]) as cursor:
    for row in cursor:
        k = int(row[0])
        auditgr = row[1]
        rcaAudit = row[2]
        # check auditgroup
        auditgr_test = redcapdict[k]["auditgr"]
        rcaAudit_test = redcapdict[k]["process_rca"]
        # verify auditgroup and kmzid
        if (auditgr_test != auditgr) or (rcaAudit_test != rcaAudit):
            print(k, rcaAudit_test == rcaAudit, auditgr_test == auditgr)
            # TODO raise an error
        oddRouteFlag = row[3]
        loopFlag = row[4]
        twoDestsFlag = row[5]
        routeLength = round(row[6], 2)
        if redcapdict[k]["oddRoute"] != oddRouteFlag:
            print("oddroute difference: {}".format(k))
        if redcapdict[k]["loopedSeg"] != loopFlag:
            print("loopedseg difference: {}".format(k))
        if redcapdict[k]["twoDests"] != twoDestsFlag:
            print("twodests difference: {}".format(format(k)))
        redcapdict[k]["Route_Length"] = routeLength


import csv

# https://stackoverflow.com/questions/29400631/python-writing-nested-dictionary-to-csv
def mergedict(a,b):
    a.update(b)
    return a


folder = r"M:\EPL-GEO_BSERE\Data\WorkingData\google.audit\RedCapUpload"
csv_columns = ["kmzid"] + list(outputdict[list(outputdict.keys())[0]].keys())
csv_file = os.path.join(folder, "IndicatorsForRedCapUpload_GoogleAudit_REGARDSBSE_2019.07.11_v1.csv")
with open(csv_file, 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_columns, delimiter=',', lineterminator='\n')
    writer.writeheader()
    for k,d in sorted(outputdict.items()):
        writer.writerow(mergedict({'kmzid': k},d))






with arcpy.da.SearchCursor("auditGroupUpdate", ["kmzid", "caresaudit", "process_rca", "auditgr", "regardsaudit"]) as cursor:
...     for row in cursor:
...         k = row[0]
...         outputdict[k] = dict()
...         outputdict[k]["caresaudit"] = row[1]
...         outputdict[k]["process_rca"] = row[2]
...         outputdict[k]["auditgr"] = row[3]
...         outputdict[k]["regardsaudit"] = row[4]



with arcpy.da.SearchCursor("Audit Routes", ["kmzidINT", "Shape_Length", "oddRoute", "auditgr", "loopedSeg", "process_rca"]) as cursor:
    for row in cursor:
        k = row[0]
        outputdict[k]["loopedSeg"] = row[4]
        outputdict[k]["oddRoute"] = row[2]
        outputdict[k]["RouteLength"] = row[1]
        testrca = row[5] == outputdict[k]["process_rca"]
        testgroup = row[3] == outputdict[k]["auditgr"]
        if testrca is False:
            print("processRCA problem on {}".format(k))
        if testgroup is False:
            print("auditgroup problem on {}".format(k))


with arcpy.da.SearchCursor("Street Segments", ["kmzidINT", "Shape_Length","auditgr"]) as cursor:
    for row in cursor:
        k = row[0]
        outputdict[k]["SegmentLength"] = row[1]
        testgroup = row[2] == outputdict[k]["auditgr"]
        if testgroup is False:
            print("auditgroup problem on {}".format(k))
