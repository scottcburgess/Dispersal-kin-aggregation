# Dispersal-kin-aggregation
R code to produce the figures in the manuscript: Burgess SC, Powell J, Bueno M. Dispersal, kin aggregation, and the fitness consequences of not spreading sibling larvae.


This README.txt file was generated on 2022-Jan-09 by Scott Burgess

GENERAL INFORMATION

1. Title of Dataset: Dispersal, kin aggregation, and the fitness consequences of not spreading sibling larvae 

2. Author Information
	A. Principal Investigator Contact Information
		Name: Scott Burgess
		Institution: Florida State University
		Address: 319 Stadium Drive, Tallahassee, FL, USA 32306
		Email: sburgess@bio.fsu.edu


3. Date of data collection (single date, range, approximate date): 2016-2017 

4. Geographic location of data collection: Turkey Point, Florida, USA

5. Information about funding sources that supported the collection of the data: National Science Foundation (NSF; OCE-1948788) and the Florida State University Council on Research and Creativity



DATA & FILE OVERVIEW

1. File List: 
Figure 2.R
Figure 3.R
Figure 4.R
Figure 5.R
Figure 2 data.csv
Figure 3 data.csv
Figure 4 data.csv
Figure 5 data.csv
Clark Evans Corrected function.R


2. Relationship between files: 
Figure 2.R uses Figure 2 data.csv
Figure 3.R uses Figure 3 data.csv and Clark Evans Corrected function.R
Figure 4.R uses Figure 4 data.csv and Clark Evans Corrected function.R
Figure 5.R uses Figure 5 data.csv


3. Metadata

Figure 2 data.csv
deployment: Sequential number for each deployment date
deployment.date: The date on which settlement plates were attached to poles in the field
treatment: Control = No colony in the center of the array; Treatment = Seven colonies placed in the center of the array.
sheet: Unique identifier for each settlement plate (=sheet)
distance.m: Distance, in meters, of the settlement plate to the center of the array
settlers: The number of settlers recorded in each settlement plate after retrieval
days: The number of days between the deployment and retrieval of settlement plates
settlers.day: settlers / days = the number of settlers per day

 


Figure 3 data.csv
D.Date: Date on which settlement plates were Deployed
R.Date: Date on which settlement plates were Retrieved
D.Group: Sequential number for each deployment group 
Plate.ID: Unique identifier for each settlement plate
Deployment: 1 = settlement plates deployed for 3 days; 2 = settlement plates placed back into the water after 3 days, and collected again after another 3 days (capturing larvae that settled on between day 4 and 6).
Point: Unique identifier for each settler
Raw.X: The distance, in millimeters, from the left side of the image
Raw.Y: The distance, in millimeters, from the bottom side of the image
True.X: The distance, in millimeters, from the left side of the focal settlement area 
True.Y: The distance, in millimeters, from the bottom side of the focal settlement area



Figure 4 data.csv
Sheet ID: Unique identifier for each settlement plate (=sheet)
Relatedness: Sib = each settler came from the same mother (half-sib or full-sib); NonSibs = each settler came from a different (unrelated) mother.
Point: Unique identifier for each settler
Raw.X: The distance, in millimeters, from the left side of the image
Raw.Y: The distance, in millimeters, from the bottom side of the image
True.X: The distance, in millimeters, from the left side of the focal settlement area 
True.Y: The distance, in millimeters, from the bottom side of the focal settlement area




Figure 5 data.csv
mother.ID: Unique identifier for each maternal colony
sib.total: A code representing the concatenation of sib.group.size and total.group.size
sib.group.size: The total number of focal siblings in the group
total.group.size: The total number of all settlers in the group
age.d: Days since settlement 
Ind: F = focal individual (an offspring from the corresponding mother.id); NF = Non-focal individual (comes from another maternal colony, NOT from the corresponding mother.id)
sheet: Unique identifier for each settlement plate (=sheet)
survival: 1 = alive; 0 = dead.
bifurcations: The number of bifurcations on the colony.
