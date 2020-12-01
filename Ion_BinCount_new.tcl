##########################################################################
##		Counts the number of select ions within			##
##			the pore of the channel				##
##									##
##	Portland State University					##
##	P.I.	: Steve Reichow						##
##	Author	: Bassam Haddad						##
##########################################################################

puts "Type the following to run program: 'run <output name> <zmin> <zmax>'"
proc run {ofile zmin zmax} {
	global ion_name
	# align pore to the center of the system
	set prot [atomselect top protein]
	set prot_cen [measure center $prot]
	set x	[expr [lindex $prot_cen 0] * -1]
	set y 	[expr [lindex $prot_cen 1] * -1]
	set z	[expr [lindex $prot_cen 2] * -1]
	set vec	[list $x $y $z]
	$prot moveby $vec
	set numframes [molinfo top get numframes]
	set output [open $ofile w]
	## Loops through each of the frames to get a count at each timestep
	for {set i 0} {$i < [expr $numframes]} {incr i} {
		animate goto 		$i
		set 	num_pot		0
		set	num_cla		0
		for {set j $zmax} {$j > $zmin} {set j [expr $j - 1]} {
			set pot [atomselect top "name POT and ((z > [expr $j - 1] and z < $j) and (x^2 + y^2 < (25)^2))"]
	                set cla [atomselect top "name CLA and ((z > [expr $j - 1] and z < $j) and (x^2 + y^2 < (25)^2))"]
			set b [expr $j -1]
			set num_pot	[expr $num_pot + [$pot num]]
			set num_cla	[expr $num_cla + [$cla num]]
		}
		puts $output "$num_pot $num_cla"
	}
	close $output
}
