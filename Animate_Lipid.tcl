proc lipid_animator {LipList molid} {
	set numframes [molinfo $molid get numframes]
	set	SegList	""
	set	ResList ""
	foreach lipid $LipList {
		set	hold 	[list [lindex $lipid 0] [lindex $lipid 1]]
		set	ResList	[concat [lindex $ResList] [lindex $hold 0]]
		unset	hold
	}
	set	ResList	[lsort -unique $ResList]
	foreach resid $ResList {
		foreach lipid $LipList {
			if {[lsearch -exact [lindex $lipid 0] $resid] >= 0} {
				lappend SegList [lindex $lipid 1]
				break
			}
		}
	}
	foreach resid $ResList segid $SegList {
		set	frame_list	""
		set	outpsf		"$resid$segid.psf"
		set	outdcd		"$resid$segid.dcd"
		foreach	lipid $LipList {
			set lipres [lindex $lipid 0]
			set lipseg [lindex $lipid 1]
			if {($lipres == $resid) && ($lipseg == $segid)} {
				lappend frame_list [lindex $lipid 4]
			}
		}
		puts [llength $frame_list]
		mol	new	[lindex [molinfo $molid get filename] 0 0]
		mol	addfile	[lindex [molinfo $molid get filename] 0 1] waitfor all
		set	sel	[atomselect top "resid $resid and segid $segid"]
		for {set i [expr $numframes -1]} {$i >= 0} {incr i -1} {
			animate goto $i
			if {[lsearch -exact $frame_list $i] >= 0} {continue
			} else {animate delete beg $i end $i top}
		}
		$sel	writepsf	$outpsf
		animate write dcd $outdcd beg 0 end -1 waitfor all sel $sel top
		puts [molinfo top get numframes]
		mol delete top
	}
}
