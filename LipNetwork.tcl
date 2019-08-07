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
#	which_center	- This evaluates which density the lipids  are closest to, and returns the density ID
#
#	eval_density	- Evaluates the average density that the lipid-tail occupys...if the avg density is >= Iso_Low, then it returns true.
#
################################
source	~/Scripts/TCL/Tcl-Scripts/matrix.tcl
################################
#
#
proc LipNetwork		{infile outfile IsoVal} {

#	This is the main program which is what will be called from the TK-Console.
#	It initalizes the matrix by inserting rows/colums, as well as starts the analysis.

	global	LipDict LipMat LipArr IsoLow NumFrames DenNum

	::struct::matrix	LipMat

	array	set		LipArr	{}

	LipMat	link		LipArr
	
	set	NumFrames	[molinfo top get numframes]

	set	IsoLow		$IsoVal

	set	ReturnList	[get_centers $infile]

	set	LipDict		[lindex $ReturnList 0]

	set	DenNum		[lindex $ReturnList 1]
	
	for	{set i 0} {$i < $DenNum} {incr i} {

		LipMat	insert	column	$i
		LipMat	insert	row	$i
	}

	get_lipid_list
	
	lip_analysis

	pop_matrix 0 0 true $outfile

	puts	"final matrix written to $outfile"
}

proc get_centers	{infile} {

#	This proc() takes a text file that contains the (x,y) coordinates of the lipid densities that are resolved by CryoEM (or any lipid density).
#	These coordinates are used to fill a dictionary with the (x,y) pair, and the lipid ID (arbitrarily defined); there is 12-fold symmetry, 
#	therefore there will be 12 "1st-lipid"'s. This dictionary will be used to evaluate where the lipids are at any frame of a trajectory.

	set	center_file	[open $infile r]

	set	centers_read	[read -nonewline $center_file]

	close	$center_file

	set	centers		[split $centers_read "\n"]

#	set	DenNum		[expr [expr [llength $centers] - 1] / 12] 
	set	DenNum		13

	set	i		0

	foreach line $centers {

		dict set LipDict		$i	x	[lindex $line 0]
		dict set LipDict		$i	y	[lindex $line 1]
		dict set LipDict		$i	id	[lindex $line 2] 

		incr i
	}

	return [list $LipDict $DenNum]

}

proc get_lipid_list	{} {

#	This proc() generates a list of lipids that ever come within a certain distance of the protein. This is done in order to 
#	limit the number of unnecessary calculations made. It checks every frame of the .dcd trajectory, selects lipids within the 
#	distance and saves a list, containing the indices of the phosphates.

	global NumFrames Phos_Ind

	set Phos_List	""

	for {set n 0} {$n < $NumFrames} {incr n} {

		animate goto $n

		set lipid	[atomselect top "resname DMPC and same residue as within 10 of (protein and resid 84 215)"]

		$lipid			set beta 1

		set	hold		[[atomselect top "name P and beta = 1"] get index]

		set	Phos_List	[concat [lindex $Phos_List] [lindex $hold]]

		unset	hold
	}

	set Phos_Ind	[lsort -unique $Phos_List]
}

proc lip_analysis	{} {

	global Phos_Ind NumFrames
	
	set tail_1_text "lipid and (name C22 to C29 C210 to C214)"
	set tail_2_text "lipid and (name C32 to C39 C310 to C316)"

	set ind_percent [expr {round([llength $Phos_Ind] / ([expr [llength $Phos_Ind] / 20]))}]
	
	set i		0

	foreach index $Phos_Ind {

		set ResID	[[atomselect top "index $index"] get resid]
		set SegID	[[atomselect top "index $index"] get segid]
		
		set tail_1	[atomselect top "resid $ResID and segid $SegID and $tail_1_text"]
		set tail_2	[atomselect top "resid $ResID and segid $SegID and $tail_2_text"]
		
		for {set n 0} {$n < $NumFrames} {incr n} {

			animate goto $n

			set	LipCenter_1	[which_center	$tail_1] 				
			set	LipCenter_2	[which_center	$tail_2]

			set	LipOccupy_1	[eval_density	$tail_1]
			set	LipOccupy_2	[eval_density	$tail_2]

			if		{$LipOccupy_1 && $LipOccupy_2}	then	{pop_matrix $LipCenter_1 $LipCenter_2

			} elseif	{$LipOccupy_1 &&! $LipOccupy_2}	then	{pop_matrix $LipCenter_1 $LipCenter_1

			} elseif	{$LipOccupy_2 &&! $LipOccupy_1}	then	{pop_matrix $LipCenter_2 $LipCenter_2

			} else	{variable donothing 0}	
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

			set LipDen	[expr [dict get $LipDict $DEN id] - 1]

		} else {	set hold	$hold

		}
	}

	return	$LipDen	
}

proc eval_density	{lipid_tail} {

	global IsoLow

	set lipid_index		[$lipid_tail get index]

	set tot_den		0

	foreach ind $lipid_index {

		set lip_atom	[atomselect top "index $ind"]

		set tot_den	[expr $tot_den + [$lip_atom get interpvol0]]
	}

	set avg_den		[expr $tot_den / [llength $lipid_index]]

	if {$avg_den >= $IsoLow} {	

		return true 

	} else	{

		return false
	}

}

proc pop_matrix		{den_id_1 den_id_2 {writematrix false} {outfile false}} {

	global LipMat LipArr DenNum

	if {$writematrix == false} {

		LipMat set cell $den_id_1 $den_id_2 [expr $LipArr($den_id_1,$den_id_2) + 1]

		LipMat set cell $den_id_2 $den_id_1 [expr $LipArr($den_id_2,$den_id_1) + 1]

	} elseif {$writematrix == true} {

		set outf	[open $outfile w]

	# row: j, column: i ... tcl matrix object has the form (column,row) which is distinct from matrices in python

		for {set i 0} {$j < $DenNum} {incr i} {
	
			for {set j 0} {$i < $DenNum} {incr j} {

				set	val	$LipArr($i,$j)

				puts	-nonewline	$outf	"$Val"
			}

			puts	$outf	"\n"
		} 

	}
	
}
