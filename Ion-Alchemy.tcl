##########################################################################
##	Transformation of selected atoms to one another			##
##									##
##	Portland State University					##
##	P.I.	: Steve Reichow						##
##	Author	: Matthew Veter						##
##########################################################################

puts			"Press '1' and select the ions that you want mutated in the VMD GUI..."
puts			"Once all of your ions have been selected, start the script by typing, 'run <output name>'"

proc	run		{ofile} {

	puts	-nonewline	"What is the name of your PSF file? "
	flush	stdout
	set	PSF_File	[gets stdin]

	puts    -nonewline	"What is the name of your PDB file? "
	flush   stdout
	set     PDB_File	[gets stdin]

	set	AtomList	[label list Atoms]

	set	n		0

	set	out		[open $ofile w]

	puts	$out		"$PSF_File\n$PDB_File"

	puts	$out		"Name\tResid\tChain\tSegID\tSerial"

	foreach	line	$AtomList {

		set	Serial	[expr [lindex $line {0 1}] + 1]
		
		set	Sel	[atomselect top "serial $Serial"]

		set	ResID	[$Sel get resid]

		set	Chain	[$Sel get chain]

		set	Name	[$Sel get resname]

		set	SegID	[$Sel get segid]

		puts	$out	"$Name\t$ResID\t$Chain\t$SegID\t$Serial"
	}

	

	close	$out

}
