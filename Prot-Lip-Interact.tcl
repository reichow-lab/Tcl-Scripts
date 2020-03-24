set     NumFrames       [molinfo top get numframes]


proc get_lipid_list     {} {

#       This proc() generates a list of lipids that ever come within a certain distance of the protein. This is done in order to 
#       limit the number of unnecessary calculations made. It checks every frame of the .dcd trajectory, selects lipids within the 
#       distance and saves a list, containing the indices of the phosphates.

        global NumFrames 
	puts "there are $NumFrames frames"
	set IndList [[atomselect top "protein and resid 189 192 and name CZ"] get index]

#	a list of all interacting lipids is created for each protein residue

	foreach INDEX $IndList {
		puts $INDEX
		set Phosp_List	""

		for {set n 0} {$n < $NumFrames} {incr n} {

			animate goto $n

			set lipid       [atomselect top "name P and within 5 of index $INDEX"]

			$lipid                  set beta 1

			set     hold            [[atomselect top "name P and beta = 1"] get index]

			set	Phosp_List	[concat [lindex $Phosp_List] [lindex $hold]]

			unset   hold

			[atomselect top all] set beta 0
		}

        	set (lip_list_$INDEX)  [lsort -unique $Phosp_List]
	}

	foreach INDEX $IndList {

		set prot_resid [[atomselect top "index $INDEX"] get resid]
		set prot_segid [[atomselect top "index $INDEX"] get segid]

		foreach LIPIND $(lip_list_$INDEX) {
			puts "$INDEX $LIPIND"
			set lip_resid [[atomselect top "index $LIPIND"] get resid]
			set lip_segid [[atomselect top "index $LIPIND"] get segid]

			set outname $prot_resid-$prot_segid-to-$lip_resid-$lip_segid
			set out	[open $outname w]
			puts $outname
			for {set n 0} {$n < $NumFrames} {incr n} {

				animate goto $n
			
	                	set prot_x      [[atomselect top "index $INDEX"] get {x}]
			 	set prot_y      [[atomselect top "index $INDEX"] get {y}]
			 	set prot_z      [[atomselect top "index $INDEX"] get {z}]

				set lip_x       [[atomselect top "index $LIPIND"] get {x}]
				set lip_y       [[atomselect top "index $LIPIND"] get {y}]
				set lip_z       [[atomselect top "index $LIPIND"] get {z}]

				set dist	[expr {sqrt(pow(($prot_x - $lip_x),2) + pow(($prot_y - $lip_y),2) + pow(($prot_z - $lip_z),2))}]

				puts	$out	"$n\t$dist"
			}

			close $out
		}
	}
}
