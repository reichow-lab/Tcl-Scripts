puts "Load molecule prior to running."
puts "To run IonMobile.tcl, type: imob <ION> <window size> <outname>"

proc imob {ION WS outname {FF 0}} {
	# set the system size and subdivide into windows
	set sys [atomselect top "all"]
	set zmin [lindex [measure minmax $sys] 0 2]
	set z [expr [lindex [measure minmax $sys] 1 2] - [lindex [measure minmax $sys] 0 2]]
	set WN [expr int($z/$WS)]
	set out [open $outname.txt w]
	# loop through window ranges
	for {set i 0} {$i <= $WN} {incr i} {
		set ionlist [atomselect top "segname ION and name $ION and (z > [expr $zmin + ($WS * $i)] and z < [expr $zmin + ($WS * ($i + 1))])"]
		set indlist [$ionlist get index]
		puts $out "Bin $i ([expr $zmin + ($WS * $i)] to [expr $zmin + ($WS * ($i + 1))])"
		foreach ind $indlist {
			animate goto $FF
			set ion [atomselect top "index $ind"]
			set pos1 [$ion get {x y z}]
			animate goto [expr $FF + 1]
			set pos2 [$ion get {x y z}]
			set dist [veclength [vecsub [lindex $pos2 0] [lindex $pos1 0]]]
			puts $out "$dist"	
			unset ion pos1 pos2 dist
		}
	}
	close $out
}
