set frame_list []

puts	"This script defaults to the top molecule. May not work as expected with more than one dcd loaded"
puts	"Use the command 'dihed_animator' to collect the frame data"

proc dihed_animator {} {

	global frame_list

	puts	"Enter four name string of atoms to check i.e. C2 C1 C3 C4: "
	flush	stdout
	gets	stdin	name_list

	puts	"Enter the number which the dihedral must be greater than: "
	flush	stdout
	gets	stdin	lower_limit
	
	puts	"Enter the number which the dihedral must be less than: "
	flush	stdout
	gets	stdin	upper_limit

	if {$lower_limit > $upper_limit} {set OR 1} else {set OR 0}

	foreach name $name_list {lappend indices [[atomselect top "name $name"] get index]}

	set numframes [molinfo top get numframes]

	for {set i [expr $numframes - 1]} {$i >= 0} {incr i -1} {

		animate goto $i

		set angle [measure dihed [split $indices]]

		if {$OR} {

			if {$angle > $lower_limit || $angle < $upper_limit} {continue
			} else {animate delete beg $i end $i top}

		} else {

			if {$angle > $lower_limit && $angle < $upper_limit} {continue
			} else {animate delete beg $i end $i top}
		}

		unset angle
	}

	unset indices upper_limit lower_limit name_list
}
