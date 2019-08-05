# Program:	LipNetwork.tcl
# Author:	Bassam Haddad
# 
# Portland State University
# Reichow Lab
#
#
#	This program takes the densities of lipids calculated from CryoEM, and looks at which lipids from MD are occupying the CryoEM lipid densities at any one point in time.
#	Ideally, the data retrieved by this code will help create a transition matrix that will objectively detail the interconnectivity of the lipid densities, allowing us to find
#	the most probably orientation and organization of lipids around the protein.
#
#	This script requires the matrix structure in TCL, created by Andreas Kupries:
#
#		Script:		https://raw.githubusercontent.com/tcltk/tcllib/master/modules/struct/matrix.tcl
#		Docs:		https://tools.ietf.org/doc/tcllib/html/matrix.html
#
#	Inputs:
#
#		LipCenters.txt	- This input will contain the info for the location and "ID's" of the lipids around the protein...Unfortunately there will be 12 sections,
#				each with multiple (equivalent) lipid centers.
#		IsoValueRange	- This will just be two values passed either as an argument to running a proc (likely) are something hard-coded into the script...essentially
#				it provides the range of CryoEM density that we are using to isolate the lipids that are occupying the space.
#		PercentWithin	- This value (0 - 100) will determine the minimum criteria for a lipid-tail to be counted as within a density...i.e. if PercentWithin is set
#				to 57, then at least 57% of the atoms of a lipid tail must be within the density for the tail to be marked as within the density. The higher
#				the value, the more strict the results.
#
#	Processes:
#
#	get-centers	- This proc() will take a text file with the the geometric (x,y) coordinates of the lipid densities around the protein.
#			The centers will need to be collected by hand (for now). These centers will be the basis for locating where the lipid densities are.
#
#	get-lipids	- Creates a list of the lipids that will be the subject of analysis. Needs Segid & Resid for each lipid. This will also contain the information
#			regarding whether or not it is POPC or DMPC...Since we mostly care about the tails.
#			- This can first make an atom selection based off of any lipids that are within a certain iso-value at anypoin in the simulation, either that or a
#			distance based procedure...which would probably be faster, however less exact?
#
#	lip-analysis	- Goes through each lipid in the list, checking where their tails are, and checking if a tail meets the criteria of being within a density...
#
#	pop-matrix	- This proc() will be called by Lip-Analysis to populate a transition matrix of lipid tail interactions (network config)
#
#	rep-lipids	- This proc() will be called by Lip-Analysis to keep track of which lipids & frames show a lipid with its tails in the densities.
#
#

source	~/Scripts/TCL/matrix.tcl

# Initialize array that the matrix will be linked to. This will allow us to easily access the matrix values.

array	set	LipArr	{}

# Initialize matrix & link to LipArr.

::struct::matrix	LipMat

LipMat link LipArr

proc get-centers	{infile} {

	set	center_file	[open $infile r]

	set	centers		[split $center_file "\n"]

	close	$center_file

	set	i	0

	foreach line $centers {

		dict set LipDic		$i	x	[lindex $line 0]
		dict set LipDic		$i	y	[lindex $line 1]
		dict set LipDic		$i	id	[lindex $line 2] 
	}

	return $LipDic

}

proc get-lipids		{} {



}

proc lip-analysis	{} {



}

proc pop-matrix		{} {



}

proc LipNetwork		{} {

set LipDic [get-centers]

}
