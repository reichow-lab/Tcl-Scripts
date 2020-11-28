# The purpose of this script it to calculate the height dependent diffusion coefficient of an Ion along the pore-axis (z-axis). By generating the distribution
# of single steps we can calculate the MSD (i.e. variance) and from that we can approximate diffusion.
puts "Load molecule prior to running."
puts "To run IonMobile.tcl, type: imob <ION> <window size> <outname>"

proc imob {ION WS outname} {
	# set the system size and subdivide into windows
	set sys [atomselect top "all"]
	set zmin [lindex [measure minmax $sys] 0 2]
	set z [expr [lindex [measure minmax $sys] 1 2] - [lindex [measure minmax $sys] 0 2]]
	set WN [expr int($z/$WS)]
	set out [open $outname.txt w]
	# loop through windows
	for {set i 0} {$i <= $WN} {incr i} {
		puts $out "BinCenter:  [expr ((($zmin + ($WS * $i)) + ($zmin + ($WS * ($i + 1)))) / 2) - [lindex [measure center [atomselect top protein]] 2]]"
		# loop through frames
		for {set j 0} {$j < [molinfo top get numframes]} {incr j} {
			set protz [lindex [measure center [atomselect top protein]] 2]
			animate goto $j
			set ionlist [atomselect top "segname ION and name $ION and (z > [expr $zmin + ($WS * $i) - $protz] and z < [expr $zmin + ($WS * ($i + 1)) - $protz])"]
			set indlist [$ionlist get index]
			foreach ind $indlist {
				animate goto $j
				set ion [atomselect top "index $ind"]
				set pos1 [$ion get {z}]
				animate goto [expr $j + 1]
				set pos2 [$ion get {z}]
				set dist [expr $pos2 - $pos1]
				puts $out "$dist"
				unset ion pos1 pos2 dist
			}
			unset ionlist indlist
		}
	}
	close $out
}
