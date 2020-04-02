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

	set holfile	[open $holeinput r]

## Loops through each of the frames to get a count at each timestep

	for {set i 0} {$i < [expr $numframes + 1]} {incr i} {

		animate goto $i

		set pot [atomselect top "name POT and ((z > $zmin and z < $zmax) and (x^2 + y^2 < 270))"]
                set cla [atomselect top "name CLA and ((z > $zmin and z < $zmax) and (x^2 + y^2 < 270))"]

		set num_pot [$POT num]
		set num_cla [$CLA num]

		puts $output "$num_pot $num_cla"

	}

	close $output

}
