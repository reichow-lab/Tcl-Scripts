#########################################################################
#		Counts the number of select ions within 		#
#			the pore of the channel				#
#									#
#			Author: Matthew Veter				#
#########################################################################


puts -nonewline "Select ion by entering its name (potassium = POT, sodium = SOD etc...).  "
flush stdout

set ion_name [gets stdin]

puts "Type the following to run program: 'run <output name>'"

proc run {ofile} {

	global ion_name

	set prot [atomselect top protein]

	set numframes [molinfo top get numframes]

	set output [open $ofile w]

## Loops through each of the frames to get a count at each timestep

	for {set i 0} {$i < [expr $numframes + 1]} {incr i} {

		animate goto $i

		set ION [atomselect top "name $ion_name and ((abs(z) < 45) and (x^2 + y^2 < 270))"]

		set num_ion [$ION num]

		puts $output "$num_ion"

	}

	close $output

}
