source /home/bassam/Scripts/TCL/Tcl-Scripts/align.tcl
puts "to use this program, type: hcomd <outfile>"

proc hcomd {outfile} {
	set numf [molinfo top get numframes]
	set outf [open $outfile-hcomd.txt w]
	for {set i 0} {$i < $numf} {incr i} {
		animate goto $i
		set sel1 [atomselect top "protein and chain A B C D E F"]
		set sel2 [atomselect top "protein and chain G H I J K L"]
		set COM1 [measure center $sel1]
		set COM2 [measure center $sel2]
		set dist [expr {sqrt(pow(([lindex $COM1 0] - [lindex $COM2 0]),2) + pow(([lindex $COM1 1] - [lindex $COM2 1]),2) + pow(([lindex $COM1 2] - [lindex $COM2 2]),2))}]
		puts $outf "$i\t$dist"
	}
	close $outf
}
