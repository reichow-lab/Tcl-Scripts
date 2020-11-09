puts "Load molecule prior to running."
puts "To run IonMobile.tcl, type: imob <ION> <windo size>"

proc imob {ION WS outname [FF 0]} {
	# set the system size and subdivide into windows
	set sys [atomselect top "water"]
	set zmin [lindex [measure minmax $sys] 0 2]
	set z [expr [lindex [measure minmax $sys] 1 2] - [lindex [measure minmax $sys] 0 2]]
	set WN [expr int($z/$WS)]
	set out [open $outname.txt w]
	# loop through window ranges
	set POS1 [list]
	set POS2 [list]
	for {set i 0} {$i < $WN} {incr i} {
		set ionlist [atomselect top "segname ION and name $ION and (z > [expr $zmin * $i] and z < [expr $zmin * ($i + 1)])"]
		set indlist [$ionlist get index]
		puts $out "Bin $i ([expr $zmin * $i] to [expr $zmin * ($i + 1)])"
		foreach ind $indlist {
			set ion [atomselect top "index $ind"]
			set pos1 [measure center $ion]
			animate goto [expr $FF + 1]
			set pos2 [measure center $ion]
			set dist [veclength [vecsub $pos2 $pos1]]
			puts $out "$dist"	
			unset ion pos1 pos2 dist
		}
	close $out
	}
}
