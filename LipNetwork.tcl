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
#		Output_name	- This is the name of the file that the final connectivity matrix will be written to...It is written in a simplistic space delimited
#				text-file, which will have an accompanying python script to make it into a numpy array.
#		IsoLow		- This is the minimum isovalue (density of volume map) that each lipid tail must occupy (on average) in order for the tail to be 
#				considered occupying a lipid density. This is a tuneable parameter...a lower isovalue accepts more lipids. Two alternatives to entering a 
#				value are: 'none', and 'sep'. 'none' takes away the density dependence and instead assignes every lipid to a center, whereas 'sep' tells
#				LipNetwork to use center-specific Isolows which are in the CentersFile.
#
################################
source	~/Scripts/TCL/Tcl-Scripts/matrix.tcl
source	~/Scripts/TCL/Tcl-Scripts/Lip-Analysis-Tools.tcl
source	~/Scripts/TCL/Tcl-Scripts/Animate_Lipid.tcl
################################

Title

proc LipNetwork		{infile outfile CarbonThreshold {IsoVal "none"} {difsel false}} {

#	This is the main program which is what will be called from the TK-Console.
#	It initalizes the matrix by inserting rows/colums, as well as starts the analysis.
#	variables:
#	
#		LipDict	:	Dictionary containing (x,y,lipID) of each lipid density, pulls data from LipCenters.txt
#		LipMat	:	The pseudo-transition matrix (connectivity matrix) linked to LipArr
#		LipArr	:	The array that the matrix is linked to, allowing me to easily call values from the matrix
#		IsoLow	:	Isovalue that is threshold for judging occupancy
#		NumFrames	Number of frames in the simulation
#		DenNum	:	Number of lipid densities (centers)	
#		MinCarbon	Minimumm number of carbons required for classification

	global	LipDict LipMat LipArr LipList IsoLow NumFrames DenNum MinCarbon OUTFILE

	::struct::matrix	LipMat

	array	set		LipArr	{}

	LipMat	link		LipArr

	set	LipList		""

	set	molid		[molinfo top]

	set	NumFrames	[molinfo top get numframes]

	set	IsoLow		$IsoVal

	set	MinCarbon	$CarbonThreshold

	set	OUTFILE		$outfile

	set	ReturnList	[get_centers $infile]

	set	LipDict		[lindex $ReturnList 0]

	set	DenNum		[lindex $ReturnList 1]
	
	for	{set i 0} {$i < $DenNum} {incr i} {

		LipMat	insert	column	$i
		LipMat	insert	row	$i
	}

	get_lipid_list
	
	lip_analysis $difsel

	pop_matrix 0 0 true $outfile

	if {$IsoLow != "none"} {

		set LipLogOutname	"$OUTFILE[set ender "_LipList.log"]"

		set LipLogOut		[open	$LipLogOutname w]

		puts	$LipLogOut	"ResID\tSegID\tLipCenter-1\tLipCenter-2\tFrame\tCarbonTail-1\tCarbonTail-2"

		foreach lipid $LipList {

			puts	$LipLogOut	"[lindex $lipid 0]\t[lindex $lipid 1]\t[lindex $lipid 2]\t\t[lindex $lipid 3]\t\t[lindex $lipid 4]\t[lindex $lipid 5]\t\t[lindex $lipid 6]"

		}

		close	$LipLogOut

	}

	puts	"final matrix written to $outfile"

	LipMat	destroy

	Title
}

proc get_centers	{infile} {

#	This proc() takes a text file that contains the (x,y) coordinates of the lipid densities that are resolved by CryoEM (or any lipid density).
#	These coordinates are used to fill a dictionary with the (x,y) pair, and the lipid ID (arbitrarily defined); there is 12-fold symmetry, 
#	therefore there will be 12 "1st-lipid"'s. This dictionary will be used to evaluate where the lipids are at any frame of a trajectory.

	set	center_file	[open $infile r]

	set	centers_read	[read -nonewline $center_file]

	close	$center_file

	set	centers		[split $centers_read "\n"]

	set	DenNum		[expr [llength $centers]]

	puts	"DenNum: $DenNum"

	set	i		1

	foreach line $centers {

		dict set LipDict		$i	x	[lindex $line 0]
		dict set LipDict		$i	y	[lindex $line 1]
		dict set LipDict		$i	id	[lindex $line 2] 
		dict set LipDict		$i	Zmin	[lindex $line 3]
		dict set LipDict		$i	Zmax	[lindex $line 4]
		dict set LipDict		$i	Isolow	[lindex $line 5]

		incr i
	}

	return [list $LipDict $DenNum]

}

proc get_lipid_list	{} {

#	This proc() generates a list of lipids that ever come within a certain distance of the protein. This is done in order to 
#	limit the number of unnecessary calculations made. It checks every frame of the .dcd trajectory, selects lipids within the 
#	distance and saves a list, containing the indices of the phosphates.

	global NumFrames Phosp_Ind

	set Phosp_List	""

	for {set n 0} {$n < $NumFrames} {incr n} {

		animate goto $n

#		The range that lipids are selected should be toggled (i.e. between 5 - 20) to make sure that you are not processing
#		an unnecessary number of calculations for a given job.

		set lipid	[atomselect top "resname DMPC and same residue as within 10 of (protein and resid 84 215)"]

		$lipid			set beta 1

		set	hold		[[atomselect top "name P and beta = 1"] get index]

		set	Phosp_List	[concat [lindex $Phosp_List] [lindex $hold]]

		unset	hold
	}

	set Phosp_Ind	[lsort -unique $Phosp_List]
}

proc lip_analysis	{difsel} {

	global Phosp_Ind NumFrames IsoLow LipList

	if {!$difsel} { 
	
		set tail_1_text "lipid and (name C22 to C29 C210 to C214)"
		set tail_2_text "lipid and (name C32 to C39 C310 to C316)"
		
	} elseif {$difsel} {

		set tail_1_text "lipid and (name C21 C23)"
		set tail_2_text "lipid and (name C31 C33)"
	}

	set ind_percent [expr {round([llength $Phosp_Ind] / ([expr [llength $Phosp_Ind] / 20]))}]
		
	set i		0

	foreach index $Phosp_Ind {

		set ResID	[[atomselect top "index $index"] get resid]
		set SegID	[[atomselect top "index $index"] get segid]
		
		set tail_1	[atomselect top "resid $ResID and segid $SegID and $tail_1_text"]
		set tail_2	[atomselect top "resid $ResID and segid $SegID and $tail_2_text"]
		
		for {set n 0} {$n < $NumFrames} {incr n} {

			animate goto $n

			set	LipCenter_1	[which_center	$tail_1] 				
			set	LipCenter_2	[which_center	$tail_2]

			if		{$IsoLow != "none"} {

				set	LipOccupy_1	[eval_density	$tail_1 $LipCenter_1]
				set	LipOccupy_2	[eval_density	$tail_2 $LipCenter_2]
				
				if	{[lindex $LipOccupy_1 0] && [lindex $LipOccupy_2 0]} {
					
					pop_matrix $LipCenter_1 $LipCenter_2

					set     attr    [list $ResID $SegID $LipCenter_1 $LipCenter_2 $n [lindex $LipOccupy_1 1] [lindex $LipOccupy_2 1]]

					lappend	LipList	$attr

					unset	attr

				} else	{variable donothing 0}

			} elseif	{$IsoLow == "none"} then	{pop_matrix $LipCenter_1 $LipCenter_2}
	
		}

		if {($i % $ind_percent) == 0} {puts -nonewline "*"}

		incr	i
	}

	puts "*"
}

proc which_center	{lipid_tail} {

	global LipDict
	
	set hold	1000000 

	set tail_COM	[measure center $lipid_tail]

	set tail_x	[lindex $tail_COM	0]

	set tail_y	[lindex $tail_COM	1]

	foreach DEN	[dict keys $LipDict] {

		set den_x	[dict get $LipDict $DEN x]
		set den_y	[dict get $LipDict $DEN y]

		set dist	[expr {sqrt(pow(($tail_x - $den_x),2) + pow(($tail_y - $den_y),2))}] 

		if {$dist < $hold} {

			set hold	$dist

			set LipDen	[dict get $LipDict $DEN id]

		} else {	set hold	$hold

		}
	}

	return	$LipDen	
}

proc eval_density	{lipid_tail lip_center} {
#
#	Currently this proc() takes the list of carbons in the atom selection, creates a list of the carbons indices and finds the
#	the density value that it occupies. It just calculates what the average density of the carbon tail is...instead I need it 
#	to calculate how many of the carbons are inside the density. This will return true if the number of carbons is greater than 
#	the required minimum.
#	
#	In addition to IsoLow, this classifier will need to know the minimum number of carbons (MinNum) to qualify for designation.
#
#	For every "true" hit, the lipid ID (RESID, not the carbon index val) and the frame of the .dcd file will be saved and stored 
#	in an output file.

	global IsoLow LipDict MinCarbon

	set lipid_index		[$lipid_tail get index]

	set tot_den		0

	set NumCarbon		0

	foreach ind $lipid_index {

		set Z_min	[dict get $LipDict $lip_center Zmin]

		set Z_max	[dict get $LipDict $lip_center Zmax]

		set lip_atom	[atomselect top "index $ind"]

		set atom_den	[$lip_atom get interpvol0]

		set atom_center	[which_center $lip_atom]

		set atom_Z	[lindex [measure center $lip_atom] 2]

		if {$IsoLow == "sep"} {set IsoThresh [dict get $LipDict $lip_center Isolow]} else {set IsoThresh $IsoLow}

		if {($atom_den >= $IsoThresh) && ($atom_center == $lip_center) && (($atom_Z >= $Z_min) &&  ($atom_Z <= $Z_max))} {incr NumCarbon}
	}
	
	if {$NumCarbon >= $MinCarbon} {	

		return [list true $NumCarbon]

	} else	{

		return false
	}

}

proc pop_matrix		{den_id_1 den_id_2 {writematrix false} {outfile false}} {

	global LipMat LipArr DenNum

	if {$writematrix == false} {

		set	den_mat_1	[expr $den_id_1 - 1]
		set	den_mat_2	[expr $den_id_2 - 1]

		LipMat set cell $den_mat_1 $den_mat_2 [expr $LipArr($den_mat_1,$den_mat_2) + 1]

		LipMat set cell $den_mat_2 $den_mat_1 [expr $LipArr($den_mat_2,$den_mat_1) + 1]

	} elseif {$writematrix == true} {

		set outf	[open $outfile w]

	# row: j, column: i ... tcl matrix object has the form (column,row) which is distinct from matrices in python

		for {set j 0} {$j < $DenNum} {incr j} {
	
			for {set i 0} {$i < $DenNum} {incr i} {

				set	val	$LipArr($i,$j)

				puts	-nonewline	$outf	"$val\t"
			}

			puts	$outf	"\n"
		} 

		close	$outf
	}
	
}
