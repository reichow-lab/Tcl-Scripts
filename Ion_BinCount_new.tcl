##########################################################################
##		Counts the number of select ions within			##
##			the pore of the channel				##
##									##
##	Portland State University					##
##	P.I.	: Steve Reichow						##
##	Author	: Bassam Haddad						##
##########################################################################


puts "Type the following to run program: 'run <output name> <zmin> <zmax> <holeinput>'"

proc run {ofile zmin zmax holeinput} {

	global ion_name

	set prot [atomselect top protein]

	set numframes [molinfo top get numframes]

	set output [open $ofile w]

	# hole file has position along pore (index: 0) and radius (index: 1)
        # create dictionary of (pore-axis, radius) pairs

	set holefile	[open $holeinput r]
	set holeread	[read -nonewline $holefile]
	set holinputs	[split $holeread "\n"]
	close $holefile

	foreach line $holinputs {

	set pore	[format %.0f [lindex $line 0]]
	set radius	[lindex $line 1]
	
	set holedict	[dict set holedict $pore $radius]
	}
	
	## Loops through each of the frames to get a count at each timestep

	for {set i 0} {$i < [expr $numframes + 1]} {incr i} {

		animate goto 		$i

		set 	num_pot		0
		set	num_cla		0

		for {set j $zmax} {$j > $zmin} {set j [expr $j - 1]} {

			set pot [atomselect top "name POT and ((z > [expr $j - 1] and z < $j) and (x^2 + y^2 < [dict get $holedict $j]^2))"]
	                set cla [atomselect top "name CLA and ((z > [expr $j - 1] and z < $j) and (x^2 + y^2 < [dict get $holedict $j]^2))"]
			
			set num_pot	[expr $num_pot + [$pot num]]	
			set num_cla	[expr $num_cla + [$cla num]]
		}
		puts $output "$num_pot $num_cla"
	}
	close $output
}
