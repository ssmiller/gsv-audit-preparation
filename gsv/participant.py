"""
"""
import arcpy


class Participant:
	def __init__(self, gccsv, sr, xfield="x", yfield="y"):
		# gccsvfolder = r"M:\EPL-GEO_BSERE\Data\WorkingData\google.audit\original_data"
		# gccsv = os.path.join(gccsvfolder, gccsv_file)
		srshortname = "" #TODO get the shortname from the spatial reference
		self.gclayername = "googleauditGC_lyr"
		arcpy.MakeXYEventLayer_management(gccsv, xfield, yfield, self.gclayername, sr)
		self.gcfcname = "GoogleAuditGeocodes_{}".format(srshortname)
		arcpy.CopyFeatures_management(self.gclayername, self.gcfcname)
		self.kmzid = "kmzid"

	def getClosestLineFromServiceArea(self, maxdist):
	    """ Extracts the closest (Euclidian) road segment from a service area around participant"""
	    ## Create Service Area for analysis
	    sa_distm = int(maxdist*1.5)  # use 1.5 times overall route distance (400m auditseg = 600m service area)
	    layer_name = "Google Audit Walking Distance Service Area {}m".format(sa_distm)
	    travel_mode_short = travel_mode.split()[0][0] + travel_mode.split()[1][0]
	    cutoffval = [sa_distm/1000.0] # convert to meters from km
	    outputtype = "LINES"
	    print("Creating Service Area...")
	    saresult_object = arcpy.MakeServiceAreaAnalysisLayer_na(
	        network, layer_name, travel_mode, "FROM_FACILITIES", cutoffval,
	        output_type=outputtype)
	    #Get the layer object from the result object. The service layer can now be
	    #referenced using the layer object.
	    salayer_object = saresult_object.getOutput(0)
	    #Get the names of all the sublayers within the service area layer.
	    sublayer_names = arcpy.na.GetNAClassNames(salayer_object)
	    #Stores the layer names that we will use later
	    facilities_layer_name = sublayer_names["Facilities"]
	    lines_layer_name = sublayer_names["SALines"]
	    #The input data has a field for "kmzid" that we want to transfer to
	    #our analysis layer. Add the field, and then use field mapping to transfer
	    #the values.
	    arcpy.na.AddFieldToAnalysisLayer(salayer_object, facilities_layer_name,
	                                     "kmzid", "TEXT", field_length=10)
	    field_mappings = arcpy.na.NAClassFieldMappings(
	        salayer_object, facilities_layer_name)
	    field_mappings[self.kmzid].mappedFieldName = self.kmzid
	    field_mappings["Name"].mappedFieldName = self.kmzid
	    # Add fields to the output Lines sublayer for later use
	    arcpy.na.AddFieldToAnalysisLayer(salayer_object, lines_layer_name,
	                                     self.kmzid, "TEXT", field_length=10)
	    #Get sublayers to work with later
	    facilities_sublayer = salayer_object.listLayers(facilities_layer_name)[0]
	    lines_sublayer = salayer_object.listLayers(lines_layer_name)[0]
	    # Add locations to service area
	    print("Adding locations...")
	    arcpy.na.AddLocations(salayer_object, facilities_layer_name, facilities,
	                          field_mappings, append="CLEAR")
	    # Solve service area
	    print("Solving service area...")
	    arcpy.na.Solve(salayer_object)
	    print(arcpy.GetMessages())
	    #Transfer the kmzid field from the input Facilities to the output Lines
	    print("Adding fields to lines...")
	    arcpy.management.AddJoin(lines_sublayer, "FacilityID",
	                             facilities_sublayer, "ObjectID")
	    #The joined fields are qualified by the feature class name of the joined
	    #table, so determine the feature class names
	    field_qualifier_lin = os.path.basename(lines_sublayer.dataSource)
	    target_field_name = "%s.kmzid" % field_qualifier_lin #TODO
	    field_qualifier_fac = os.path.basename(facilities_sublayer.dataSource)
	    expression = "!%s.kmzid!" % field_qualifier_fac #TODO
	    print(expression)
	    arcpy.management.CalculateField(lines_sublayer, target_field_name,
	                                            expression, "PYTHON")
	    arcpy.management.RemoveJoin(lines_sublayer)
	    # Export the lines for further analysis
	    print("Copying Service Area Lines to new feature class...")
	    out_lines = "SAWDlines_googleaudit_{}m".format(sa_distm)
	    arcpy.management.CopyFeatures(lines_sublayer, out_lines)
	    # Dissolve and Unsplit the lines so that there are intersection-to-intersection road segments,
	    # Especially near the participant (would have been split at their location)
	    print("Dissolving Lines")
	    out_lines_dslv = out_lines + "_dslv"
	    arcpy.Dissolve_management(out_lines, out_lines_dslv, dissolve_field=self.kmzid,
	                              unsplit_lines="UNSPLIT_LINES")
	    # Reproject dissolved roads maxdist * 1.5 (e.g. 600m) to USCEC
	    print("Reprojecting dissolved lines")
	    srCEC = arcpy.SpatialReference("USA Contiguous Equidistant Conic")
	    SAlines_cec = out_lines_dslv + "_cec"
	    arcpy.Project_management(out_lines_dslv, SAlines_cec, srCEC)
	    print(arcpy.GetMessages())
	    # Create a near table of closest 15 roads within the dissolved dataset
	    print("Identifying nearest road")
	    gcCEC = "_".join(facilities.split("_")[:-1] + ["cec"])
	    if not arcpy.Exists(gcCEC):
	        arcpy.Project_management(facilities, gcCEC, srCEC)
	    neartable = "Nearest15Roads_{}m{}".format(sa_distm, travel_mode_short)
	    arcpy.GenerateNearTable_analysis(gcCEC, SAlines_cec, neartable,
	                                     closest="ALL", closest_count=15)
	    ### Identify the matching closest road for the kmzid
	    # Add a field to the nearest roads table to store selectedroad
	    arcpy.AddField_management(neartable, field_name="selectedrd",
	                              field_type="SHORT")
	    # Add a field to the nearest road table to store kmzid associated with IN_FID
	    arcpy.AddField_management(neartable, field_name="kmzidgc",
	                              field_type="TEXT", field_length=10)
	    # Add a field to the nearest road table to store kmzid associated with NEAR_FID
	    arcpy.AddField_management(neartable, field_name="kmzidrd",
	                              field_type="TEXT", field_length=10)
	    # Create a table view for join
	    tblview = arcpy.MakeTableView_management(neartable, "neartable_view")
	    # Add joins
	    desc_nt = arcpy.Describe(neartable)
	    field_qualifier_nt = desc_nt.baseName
	    desc_gc = arcpy.Describe(gcCEC)
	    field_qualifier_gc = desc_gc.baseName
	    desc_sa = arcpy.Describe(SAlines_cec)
	    field_qualifier_sa = desc_sa.baseName
	    # Add join between near table and geocodes and calculate kmzidgc field based on join
	    arcpy.AddJoin_management(tblview, "IN_FID", gcCEC, "OBJECTID")
	    fieldtocalculate = "{}.kmzidgc".format(field_qualifier_nt)
	    expr_gc = "!{}.kmzid!".format(field_qualifier_gc)
	    arcpy.CalculateField_management(tblview, fieldtocalculate, expr_gc, "PYTHON")
	    arcpy.RemoveJoin_management(tblview)
	    # Add join between near table and roadsa and calculate kmzidrd field based on join
	    arcpy.AddJoin_management(tblview, "NEAR_FID", SAlines_cec, "OBJECTID")
	    fieldtocalculate = "{}.kmzidrd".format(field_qualifier_nt)
	    expr_sa = "!{}.kmzid!".format(field_qualifier_sa)
	    arcpy.CalculateField_management(tblview, fieldtocalculate, expr_sa, "PYTHON")
	    arcpy.RemoveJoin_management(tblview)
	    # Identify first match (near table is ordered by near_rank)
	    fieldlist = ["kmzidgc", "kmzidrd", "selectedrd"]
	    roadfound = []
	    with arcpy.da.UpdateCursor(tblview, fieldlist) as cursor:
	        for row in cursor:
	            if row[0] not in roadfound and row[0] == row[1]:
	                roadfound.append(row[0])
	                row[2] = 1 # road is selected as closest to participant
	            else:
	                row[2] = 0 # road is not selected as closest road
	            cursor.updateRow(row)
	    # Now have the table with a field of selectedrd with 1 if the road should be selected out.
	    whereclause = "{} = 1".format(arcpy.AddFieldDelimiters(neartable, "selectedrd"))
	    with arcpy.da.SearchCursor(neartable, ["IN_FID", "selectedrd", "NEAR_FID"], whereclause) as cursor:
	        roadfids = [row[2] for row in cursor]
	        # TODO maybe remove the additional fields as they are not needed.
	    # print(len(roadfids))
	    print("Copying the nearest roads to new feature class")
	    roadfidsstr = ",".join(map(str, roadfids))
	    outputfd = "ServiceAreaLines_CEC"
	    # TODO create feature dataset
	    if not arcpy.Exists(os.path.join(arcpy.env.workspace, outputfd)):
	        arcpy.CreateFeatureDataset_management(arcpy.env.workspace, outputfd, srCEC)
	    else:
	        print("FD {} already exists".format(outputfd))
	    # Create a feature class for the Google Audit locations
	    whereclauserd = "{} IN ({})".format(arcpy.AddFieldDelimiters(SAlines_cec, "OBJECTID"),roadfidsstr)
	    outfc = os.path.join(
	        outputfd,
	        "SelectedClosestSegmentsfrom{}m{}SA_all".format(sa_distm, travel_mode_short))
	    arcpy.Select_analysis(SAlines_cec, outfc, whereclauserd)
	    # TODO possibly add other variables to pass back out.
	    # This output is all of the closest road segments for the participants.
	    return outfc #pass back the name of the segment output
