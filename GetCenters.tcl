puts			"Press '1' and select the waters that you used for selecting lipid centers..."
puts			"Once all of your waters have been selected, start the script by typing, 'run <output name>'"

proc	run		{ofile LipNum IsoList zmin1 zmax1 zmin2 zmax2} {

	set	AtomList	[label list Atoms]

	set	out		[open $ofile w]

	set	lipn		1
	
	set	N		1

#	The format of the output will be	[xval	yval	lipn	zmin	zmax	IsoLow]
	
	foreach iso	$IsoList {

		dict	set	IsoVals	$lipn	Isolow	[lindex $iso]

		incr	lipn
	}

	set	lipn		1

	foreach	line	$AtomList {
		
		if {$N <= [expr $LipNum * 6]} {
	
			set	zmin	$zmin1
			set	zmax	$zmax1

		} else {

			set	zmin	$zmin2
			set	zmax	$zmax2 }


		set	Serial	[expr [lindex $line {0 1}] + 1]
		
		set	Sel	[atomselect top "serial $Serial"]

		set	xval	[$Sel get x]

		set	yval	[$Sel get y]

		set	isolow	[dict get $IsoVals $lipn Isolow]

		puts	$out	"$xval\t$yval\t$lipn\t$zmin\t$zmax\t$isolow"
		
		incr	N

		incr	lipn

		if {$lipn >= [expr $LipNum + 1]} {set lipn 1}

	}

	

	close	$out

}
