# The purpose of this script it to calculate the height dependent diffusion coefficient of an Ion along the pore-axis (z-axis). By generating the distribution
# of single steps we can calculate the MSD (i.e. variance) and from that we can approximate diffusion.
source ~/Scripts/TCL/Tcl-Scripts/auto-ionz-IN.tcl
puts "Load molecule prior to running."
puts "To run IonMobile.tcl, type: imob <ION> <window size> <outname>"

proc imob {ION WS outname} {
	set ID [molinfo top get id]
	align 0 $ID
	set protz [lindex [measure center [atomselect top protein]] 2]
	# set the system size and subdivide into windows
	set sys [atomselect 0 "all"]
	set zmin [lindex [measure minmax $sys] 0 2]
	set z [expr [lindex [measure minmax $sys] 1 2] - [lindex [measure minmax $sys] 0 2]]
	set WN [expr int($z/$WS)]
	set out [open $outname.txt w]
	# loop through windows
	for {set i 0} {$i <= $WN} {incr i} {
		puts $out "BinCenter:\t[expr ((($zmin + ($WS * $i)) + ($zmin + ($WS * ($i + 1)))) / 2) - $protz]"
		# loop through frames
		for {set j 0} {$j < [molinfo top get numframes]} {incr j} {
			animate goto $j
			# Select ions in the current bin
			set ionlist [atomselect top "segname ION and name $ION and (z > [expr $zmin + ($WS * $i) - $protz] and z < [expr $zmin + ($WS * ($i + 1) - $protz)])"]
			set indlist [$ionlist get index]
			foreach ind $indlist {
				animate goto $j
				# in order to correct for the PBC without having to unwrap my simulation
				set zcorr [expr [lindex [measure minmax $sys] 1 2] - [lindex [measure minmax $sys] 0 2]]
				set zhalf [expr $zcorr / 2]
				set ion [atomselect top "index $ind"]
				set posZ1 [$ion get {z}]
				set posX1 [$ion get {x}]
				set posY1 [$ion get {y}]
				animate goto [expr $j + 1]
				set posZ2 [$ion get {z}]
				set posX2 [$ion get {x}]
				set posY2 [$ion get {y}]
				set distZ [expr $posZ2 - $posZ1]
				if {abs($distZ) >= $zhalf} {
					if {$distZ >= $zhalf} {set distZ [expr $distZ - $zcorr]
					} elseif {$distZ <= -$zhalf} {set distZ [expr $distZ + $zcorr]}
				}
				set distX [expr $posX2 - $posX1]
				if {abs($distX) >= $zhalf} {
					if {$distX >= $zhalf} {set distX [expr $distX - $zcorr]
					} elseif {$distX <= -$zhalf} {set distX [expr $distX + $zcorr]}
				}
				set distY [expr $posY2 - $posY1]
				if {abs($distY) >= $zhalf} {
					if {$distY >= $zhalf} {set distZ [expr $distY - $zcorr]
					} elseif {$distY <= -$zhalf} {set distY [expr $distY + $zcorr]}
				}
				puts $out "$distZ\t$distX\t$distY"
				unset ion posZ1 posX1 posY1 posZ2 posX2 posY2 distZ distX distY zcorr
			}
			unset ionlist indlist
		}
	}
	close $out
}
