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

	for {set i 0} {$i < $numframes} {incr i} {

		animate goto $i

		set angle [measure dihed [split $indices]]

		if {$OR} {

			if {$angle > $lower_limit || $angle < $upper_limit} {lappend frame_list $i}

		} else {

			if {$angle > $lower_limit && $angle < $upper_limit} {lappend frame_list $i}
		}

		unset angle
	}

	set dcd_len [llength $frame_list]

	puts "The expected number of frames is $dcd_len]"

	unset indices upper_limit lower_limit

	for {set i [expr $numframes -1]} {$i >= 0} {incr i -1} {

		if {[lsearch -exact $frame_list $i] >= 0} {continue
		} else {animate delete beg $i end $i top
		}
	}

	set sys [atomselect top all]

	set ref_frame [atomselect top "name $name_list" frame 0]

	for {set i 0} {$i < $dcd_len} {incr i} {

		animate goto $i

		set align_frame [atomselect top "name $name_list"]

		set trans_matrix [measure fit $align_frame $ref_frame]

		$sys move $trans_matrix
	}

	unset sys ref_frame name_list
}
