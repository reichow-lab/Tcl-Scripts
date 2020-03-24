##########################################################################
##		Counts the number of select ions within			##
##			the pore of the channel				##
##									##
##	Portland State University					##
##	P.I.	: Steve Reichow						##
##	Author	: Matthew Veter						##
##########################################################################


#puts -nonewline "Select ion by entering its name (potassium = POT, sodium = SOD etc...).  "
#flush stdout

#set ion_name [gets stdin]

puts "Type the following to run program: 'run <output name> <zmin> <zmax>'"

proc run {ofile zmin zmax} {

	global ion_name

	set prot [atomselect top protein]

	set numframes [molinfo top get numframes]

	set output [open $ofile w]

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
