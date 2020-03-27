set     NumFrames       [molinfo top get numframes]
set	chains_upper	"A B C D E F"
set	chains_lower	"G H I J K L"

proc run     {} {

#       This proc() generates a list of lipids that ever come within a certain distance of the protein. This is done in order to 
#       limit the number of unnecessary calculations made. It checks every frame of the .dcd trajectory, selects lipids within the 
#       distance and saves a list, containing the indices of the phosphates.

        global NumFrames chains_upper chains_lower 
	puts "there are $NumFrames frames"
	set IndList [[atomselect top "protein and resid 183 192 and name CZ"] get index]

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

			set	hold		[[atomselect top "name P and beta = 1"] get index]

			set	Phosp_List	[concat [lindex $Phosp_List] [lindex $hold]]

			unset	hold

			[atomselect top all] set beta 0
		}

		set (lip_list_$from-$to)	[lsort -unique $Phosp_List]
	}

	for {set i 0} {$i < 12} {incr i} {

		if {$i < 6}     {set from [lindex $chains_upper [expr $i % 6]]; set to [lindex $chains_upper [expr ($i + 1) % 6]]
		} else            {set from [lindex $chains_lower [expr $i % 6]]; set to [lindex $chains_lower [expr ($i + 1) % 6]]}

		foreach LIPIND $(lip_list_$from-$to) {

			puts "$LIPIND"
			set lip_resid [[atomselect top "index $LIPIND"] get resid]
			set lip_segid [[atomselect top "index $LIPIND"] get segid]

			set outname Arg192$from-Arg183$to-to-$lip_resid-$lip_segid
			set out	[open $outname w]
			puts $outname
			for {set n 0} {$n < $NumFrames} {incr n} {

				animate goto $n
			
#				set	R192_NH1	[measure bond "$LIPIND [[atomselect top "protein and resid 192 and name NH1 and chain $from"] get index]"]
#				set	R192_NH2	[measure bond "$LIPIND [[atomselect top "protein and resid 192 and name NH2 and chain $from"] get index]"]
#				set	R192_NE		[measure bond "$LIPIND [[atomselect top "protein and resid 192 and name NE and chain $from"] get index]"]
				set	R192_CZ		[measure bond "$LIPIND [[atomselect top "protein and resid 192 and name CZ and chain $from"] get index]"]
#				set	R183_NH1	[measure bond "$LIPIND [[atomselect top "protein and resid 183 and name NH1 and chain $from"] get index]"]
#				set	R183_NH2	[measure bond "$LIPIND [[atomselect top "protein and resid 183 and name NH2 and chain $from"] get index]"]
#				set	R183_NE		[measure bond "$LIPIND [[atomselect top "protein and resid 183 and name NE and chain $from"] get index]"]
				set	R183_CZ		[measure bond "$LIPIND [[atomselect top "protein and resid 183 and name CZ and chain $from"] get index]"]

#				puts	$out	"$n\t$R192_NH1\t$R192_NH2\t$R192_NE\t$R192_CZ\t$R183_NH1\t$R183_NH2\t$R183_NE\t$R183_CZ"
				puts	$out	"$n\t$R192_CZ\t$R183_CZ"
				
				unset	R192_CZ
				unset	R183_CZ
			}

			close $out
		}
	}
}
