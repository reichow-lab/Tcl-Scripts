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
	set ref [atomselect 0 "all"]
	set sys [atomselect top "all"]
	set zmin [lindex [measure minmax $ref] 0 2]
	set z [expr [lindex [measure minmax $ref] 1 2] - $zmin]
	set WN [expr int($z/$WS)]
	set out [open $outname.txt w]
	# loop through windows
	for {set i 0} {$i <= $WN} {incr i} {
		puts $out "BinCenter:\t[expr ((($zmin + ($WS * $i)) + ($zmin + ($WS * ($i + 1)))) / 2) - $protz]"
		# loop through frames
		set numF [molinfo top get numframes]
		for {set j 0} {$j < $numF} {incr j} {
			animate goto $j
			# in order to correct for the PBC without having to unwrap my simulation
			set sysdim [measure minmax $sys]
			set zcorr [expr [lindex $sysdim 1 2] - [lindex $sysdim 0 2] - 8]
			set xcorr [expr [lindex $sysdim 1 0] - [lindex $sysdim 0 0]]
			set ycorr [expr [lindex $sysdim 1 1] - [lindex $sysdim 0 1]]
			set zlim [expr $zcorr / 1.4]
			set xlim [expr $xcorr / 1.4]
			set ylim [expr $ycorr / 1.4]
			# Select ions in the current bin
			set ionlist [atomselect top "segname ION and name $ION and (z > [expr $zmin + ($WS * $i) - $protz] and z < [expr $zmin + ($WS * ($i + 1) - $protz)])"]
			set indlist [$ionlist get index]
			foreach ind $indlist {
				animate goto $j
				set ion [atomselect top "index $ind"]
				set posZ1 [$ion get {z}]
				set posX1 [$ion get {x}]
				set posY1 [$ion get {y}]
				animate goto [expr $j + 1]
				set posZ2 [$ion get {z}]
				set posX2 [$ion get {x}]
				set posY2 [$ion get {y}]
				set distZ [expr $posZ2 - $posZ1]
				if {abs($distZ) >= $zlim} {
					if {$distZ >= $zlim} {set distZ [expr $distZ - $zcorr]
					} elseif {$distZ <= -$zlim} {set distZ [expr $distZ + $zcorr]}
				}
				set distX [expr $posX2 - $posX1]
				if {abs($distX) >= $zlim} {
					if {$distX >= $xlim} {set distX [expr $distX - $xcorr]
					} elseif {$distX <= -$xlim} {set distX [expr $distX + $xcorr]}
				}
				set distY [expr $posY2 - $posY1]
				if {abs($distY) >= $ylim} {
					if {$distY >= $ylim} {set distY [expr $distY - $ycorr]
					} elseif {$distY <= -$zlim} {set distY [expr $distY + $ycorr]}
				}
				puts $out "$distZ\t$distX\t$distY"
				unset posZ1 posX1 posY1 posZ2 posX2 posY2 distZ distX distY
				$ion delete
			}
			unset indlist sysdim zcorr xcorr ycorr zlim xlim ylim
			$ionlist delete
		}
	}
	close $out
}
