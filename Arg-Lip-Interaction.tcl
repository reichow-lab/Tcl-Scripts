set     NumFrames       [molinfo top get numframes]
set	chains_upper	"A B C D E F"
set	chains_lower	"G H I J K L"

proc run     {} {

#       This proc() generates a list of lipids that ever come within a certain distance of the protein. This is done in order to 
#       limit the number of unnecessary calculations made. It checks every frame of the .dcd trajectory, selects lipids within the 
#       distance and saves a list, containing the indices of the phosphates.

        global NumFrames chains_upper chains_lower 
	puts "there are $NumFrames frames"
	set IndList 	[[atomselect top "protein and resid 183 192 and name CZ"] get index]
	set all		[atomselect top all]


#	a list of all interacting lipids is created for each protein residue

	for {set i 0} {$i < 12} {incr i} {

		if {$i < 6}	{set from [lindex $chains_upper [expr $i % 6]]; set to [lindex $chains_upper [expr ($i + 1) % 6]]
		} else 		{set from [lindex $chains_lower [expr $i % 6]]; set to [lindex $chains_lower [expr ($i + 1) % 6]]}
		puts "$from to $to"
		set Phosp_List ""

		for {set n 0} {$n < $NumFrames} {incr n} {
		
			animate goto $n

			set lipid	[atomselect top "name P and within 7 of ((chain $from and resid 192) or (chain $to and resid 183))"]

			$lipid		set beta 1

			set	sel		[atomselect top "name P and beta = 1"]

			set	hold		[$sel get index]

			set	Phosp_List	[concat [lindex $Phosp_List] [lindex $hold]]

			$lipid	delete
			$sel	delete

			$all set beta 0
		}

		set (lip_list_$from-$to)	[lsort -unique $Phosp_List]
	}

	for {set i 0} {$i < 12} {incr i} {

		if {$i < 6}     {set from [lindex $chains_upper [expr $i % 6]]; set to [lindex $chains_upper [expr ($i + 1) % 6]]
		} else            {set from [lindex $chains_lower [expr $i % 6]]; set to [lindex $chains_lower [expr ($i + 1) % 6]]}

		foreach LIPIND $(lip_list_$from-$to) {

			puts "$LIPIND"
			set lipid	[atomselect top "index $LIPIND"]
			set lip_resid	[$lipid get resid]
			set lip_segid 	[$lipid get segid]

			set outname Arg192$from-Arg183$to-to-$lip_resid-$lip_segid
			set out	[open $outname w]
			puts $outname
			for {set n 0} {$n < $NumFrames} {incr n} {

				animate goto $n

				set	R192_CZ		[atomselect top "protein and resid 192 and name CZ and chain $from"]
				set	R192_P		[measure bond "$LIPIND [$R192_CZ get index]"]
				set	R183_CZ		[atomselect top "protein and resid 183 and name CZ and chain $from"]
				set	R183_P		[measure bond "$LIPIND [$R183_CZ get index]"]

				puts	$out	"$n\t$R192_P\t$R183_P"
				
				unset		R192_P
				unset		R183_P
				$R192_CZ	delete
				$R183_CZ	delete
			}

			close $out
		}
	}
}
