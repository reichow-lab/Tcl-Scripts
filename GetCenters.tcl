puts			"Press '1' and select the waters that you used for selecting lipid centers..."
puts			"Once all of your waters have been selected, start the script by typing, 'run <output name>'"

proc	run		{ofile} {

	set	AtomList	[label list Atoms]

	set	out		[open $ofile w]

	set	lipn		1
	
	set	N		1

#	The format of the output will be	[xval	yval	lipn	zmin	zmax]

	foreach	line	$AtomList {
		
		if {$N <= 72} {
	
			set	zmin	115
			set	zmax	135 

		} else {

			set	zmin	177
			set	zmax	197 }


		set	Serial	[expr [lindex $line {0 1}] + 1]
		
		set	Sel	[atomselect top "serial $Serial"]

		set	xval	[$Sel get x]

		set	yval	[$Sel get y]

		puts	$out	"$xval\t$yval\t$lipn\t$zmin\t$zmax"
		
		incr	N

		incr	lipn

		if {$lipn >= 13} {set lipn 1}

	}

	

	close	$out

}
