# gsv-audit-preparation
Google Street View Audit Preparation

Given inputs of geocode and network dataset, this produces several intermediate outputs and four final outputs:
1) street segment, up to a maximum desired length, closest to each geocode (line features)
2) audit route (includes street segment) up to a maximum desired length. The audit route incorporates one or more
street segments or partial street segments for each geocode (line features)
3) side of street: The side of the street of the original geocode, placed at the midpoint of each street segment 
(point features)
4) endpoints of each street segment (point features)

The code currently assumes specific field names for unique identifiers for the geocodes.

Requirements: ArcGIS Pro 2.x +, Python 3.x with arcpy, license for Network Analyst and
data license for StreetMap Premium. This network dataset does not allow direct access to the
streets included so a workaround is utilized within the code. The network dataset also has 
specific Travel Modes built in (e.g. Walking Distance, Driving Distance...) that are assumed
in the code.



