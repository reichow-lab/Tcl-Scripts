##########################################################################
##  An all-encompassing script with different processes of MD analysis	##
##									##
##	Portland State University					##
##	P.I.	: Steve Reichow						##
##	Author	: Matthew Veter						##
##########################################################################

proc title {{opt 0}} {
	puts "This script contains the following analyses:"
	puts " **Command**\t\t**Description**"
	if {$opt} {puts "<align>           Aligns the selection to itself throughout frames"}
	puts "<rmsd>            Measure the RMSD of a selection. Don't forget to <align>"
	puts "<rmsf>            Measure the RMSF of a selection. Don't forget to <align>"
	puts "<radius_gyro>     Measure the radius of gyration of a selection"
	puts "<sasa {srad}>     Measure SASA of a selection. Default radius is 1.4 Angstroms"
	puts "<dihedral>        Measure the dihedral angles of a selection"
	puts "<title {bool}>    To show these commands again."
}
title


## Recursive functions used throughout
proc get_molid {} {

	puts	-nonewline	"Are you working with the 'top' molecule? <\[y] or \[n]> "
	flush	stdout
	set	top_check	[gets stdin]

	if {$top_check == "n"} {
		puts -nonewline "What is the molid? "; flush stdout; set molecule_id [gets stdin]
	}	elseif {$top_check == "y"} {set molecule_id "top"
	}	else {error "The only options are <\[y] or \[n]>"}

	return	$molecule_id
}

proc get_selection {} {

	puts	-nonewline	"Type your atom selection i.e. \[atomselect molid <your input>]: "
	flush	stdout
	set	atom_selection	[gets stdin]

	return	$atom_selection
}

proc get_outfile {} {

	puts	-nonewline	"Name your output file: "
	flush	stdout
	set	outfile_name	[gets stdin]

	return	$outfile_name
}

proc get_rmsinfo {} {

	puts	-nonewline	"Provide molID for reference. "
	flush	stdout
	set	ref_molid	[gets stdin]

	puts	-nonewline	"Provide molID for selection. "
	flush	stdout
	set	sel_molid	[gets stdin]
	set	selection	[get_selection]
	return	[list $ref_molid $sel_molid $selection]
}

proc align {} {

	# Setting initial variables
	lassign	[get_rmsinfo]	ref_molid sel_molid selection
	set	ref_frame	[atomselect $ref_molid $selection frame 0]
	set	sys		[atomselect $sel_molid all]
	set	tot_dec		[expr ($finaframe - $initframe + 1) / 10]
	set	n		1
	set	q		10

	# Aligns each frame of the selection to initial frame
	for {set f $initframe} {$f < $finaframe} {incr f} {

		animate goto $f
		set	align_frame	[atomselect $sel_molid $selection]
		set	trans_matrix	[measure fit $align_frame $ref_frame]
		$sys	move		$trans_matrix

		if {($n % $tot_dec) == 0} {
			puts	"Alignment $q% complete"
			set	q	[expr $q + 10]
		}
		incr	n
	}
	puts	"Alignment complete, ready for calculations"
}


## Measures and writes to file the RMSD of select atoms per timestep
proc rmsd {} {

	# Setting initial variables
	set	molid		[get_molid]
	set	selection	[get_selection]
	set	sel		[atomselect $molid $selection]
	set	outfile		[get_outfile]
	set	outf		[open $outfile w]
	set	frames		[molinfo $molid get numframes]
	set	rmsd_sel	[atomselect $molid $selection]
	set	rmsd_ref	[atomselect $molid $selection frame 0]
	
	# Measuring RMSD through each frame
	for {set f 0} {$f < $frames} {incr f} {

		animate goto $f
		set	rmsd_struc	[atomselect $molid $selection]
		set	rmsd_calc	[measure rmsd $rmsd_struc $rmsd_ref]
		puts	$outf		"$f\t$rmsd_calc"
		unset	rmsd_struc rmsd_calc
	}
	unset molid selection sel outfile frames rmsd_sel rmsd_ref
	close $outf
}


## Measures and writes to file the RMSF of select atoms per timestep
proc rmsf {} {

	# Setting initial variables
	set	molid		[get_molid]
	set	selection	[get_selection]
	set	sel		[atomselect $molid $selection]
	set	outfile		[get_outfile]
	set	outf		[open $outfile w]
	set	indices		[$sel get index]

	foreach index $indices {

		set	atom		[atomselect $molid "index $index"]
		set	rmsf_calc	[measure rmsf $atom]
		puts	$outf		"$atom\t$rmsf_calc"
		unset	atom rmsf_calc
	}
	unset	molid selection sel outfile indices
	close	$outf
}


## Measures and writes to file the solvent accessible surface area (SASA) of select atoms per timestep
proc sasa {{srad 1.4}} {

	# Setting initial variables
	set	molid		[get_molid]
	set	selection	[get_selection]
	set	sel		[atomselect $molid $selection]
	set	outfile		[get_outfile]
	set	outf		[open $outfile w]
	set	frames		[molinfo $molid get numframes]

	# Looping through each frame
	for {set f 0} {$f < $frames} {incr f} {

		animate goto $f
		set	area	[measure sasa $srad $sel]
		puts	$outf	[format "$f\t$area"]
		unset	area
	}
	unset	molid selection sel srad outfile frames
	close	$outf
}


## Measures and writes to file the radius of gyration of select atoms per timestep
proc radius_gyro {} {

	set	molid		[get_molid]
	set	selection	[get_selection]
	set	sel		[atomselect $molid $selection]
	set	outfile		[get_outfile]
	set	outf		[open $outfile w]
	set	frames		[molinfo $molid get numframes]

	if {[$sel num] <= 0} {
        	error "radius_of_gyration: requires a selection of at least one atom"
	}

	for {set i 0} {$i < $frames} {incr i} {

        	animate goto $i
	        set	rad	[measure rgyr $sel]
        	puts	$outf	[format "$i\t$rad"]
	}
	unset	molid selection sel outfile frames
	close	$outf
}


## Measures and writes to file the dihedral angles of select atoms per timestep
proc dihedral {} {

	label	delete	Atoms
	label	delete	Dihedrals
	
	puts "Using the GUI and your mouse:"
	puts "Press '4' on the keyboard and select the four atoms of interest..."
	puts "Once your dihedral angle has been labelled, start the script by typing <run>"
	puts "Continue using the command selecting new angles as desired."

	proc run {} {

		set	outfile		[get_outfile]
		set	outf		[open $outfile w]
		set	dihed		[label graph Dihedrals 0]
		set	n		0

		foreach angle $dihed {

			puts	$outf	[format "$n\t$angle"]
			incr	n
		}

		label	delete	Atoms
		label	delete	Dihedrals
		unset	outfile dihed n
		close	$outf
	}
}


## Measures and/or writes to file the center of mass of select atoms per timestep
proc center_of_mass {{opt 0}} {

	# Setting initial variables
	set	molid		[get_molid]
	set	selection	[get_selection]
	set	sel		[atomselect $molid $selection]
	if	{[$sel num] <= 0} {error "center_of_mass: needs a \selection with atoms"}

	if {$opt} {

		# Writes the trajectory of the selection by means of center of mass
		set	outfile		[get_outfile]
		set	outf		[open $outfile w]
		set	frames		[molinfo $molid get numframes]

		for {set f 0} {$f < $frames} {incr f} {

			animate goto $f
			set	com	[veczero]
			set	mass	0

			foreach coord [$sel get {x y z}] m [$sel get mass] {

				set	mass	[expr $mass + $m]
				set	com	[vecadd $com [vecscale $m $coord]]
			}
			set	com_traj	[vecscale [expr 1.0/$mass] $com]
			puts	$outf		"$com_traj"
		}
	} else {

		set	com	[veczero]
		set	mass	0

		# Computes weighted center by mass of atoms
		foreach coord [$sel get {x y z}] m [$sel get mass] {

			set	mass	[expr $mass + $m]
			set	com	[vecadd $com [vecscale $m $coord]]
		}

		return [vecscale [expr 1.0/$mass] $com]
	}
	unset molid selection sel com mass
}
