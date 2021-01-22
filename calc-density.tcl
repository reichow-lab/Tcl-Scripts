source	~/Scripts/TCL/Tcl-Scripts/auto-ionz-IN.tcl
set	PO_tail_text "lipid and (name C22 to C29 C210 to C218 C32 to C39 C310 to C316)"
set	DM_tail_text "lipid and (name C22 to C29 C210 to C214 C32 to C39 C310 to C314)"
set PC_head_text "lipid and (name O21 O22 O31 O32 O11 to O14 C1 C2 C21 C3 C31 C11 to C15 P N)"
set	prot_text    "protein and noh"

puts	"To run this program type ~ iondensity <ofile> ~ or ~ popcdensity <ofile> ~ or ~ dmpcdensity <ofile>."

proc iondensity {ofile} {

	set	wat [atomselect top water]
	set	watend "_watden.dx"

#	puts	-nonewline "Select anion by entering its name (chloride = CLA)"
#	flush	stdout
#	set	anion_name [gets stdin]
	set	cl [atomselect top "name CLA and segname ION"]
	set	clend "_anden.dx"

#	puts	-nonewline "Select cation by entering its name (potassium = POT, sodium = SOD)"
#	flush	stdout
#	set	cation_name [gets stdin]
	set	s [atomselect top "name SOD and segname ION"]
	set	k [atomselect top "name POT and segname ION"]
	set	c [atomselect top "name CES and segname ION"]
	set ca [atomselect top "name CAL"]
	set	kend "_potden.dx"
	set	send "_sodden.dx"
	set	cend "_cesden.dx"
	set caend "_calden.dx"
	#puts "Beginning water density calculation."

	#volmap density $wat -allframes -combine avg -res 0.649 -o $ofile$watend

	puts "Starting chloride density calculation."

	volmap density $cl -allframes -combine avg -res 0.649 -o $ofile$clend

	puts "Finished chloride, starting cation density calculation."

	if {[$s num] != 0} {volmap density $s -allframes -combine avg -res 0.649 -minmax [list {-100 -100 -100} {100 100 100}] -o $ofile$send}
	if {[$k num] != 0} {volmap density $k -allframes -combine avg -res 0.649 -minmax [list {-100 -100 -100} {100 100 100}] -o $ofile$kend}
	if {[$c num] != 0} {volmap density $c -allframes -combine avg -res 0.649 -minmax [list {-100 -100 -100} {100 100 100}] -o $ofile$cend}
  if {[$ca num] !=0} {volmap density $ca -allframes -combine avg -res 0.649 -minmax [list {-100 -100 -100} {100 100 100}] -o $ofile$caend}
	puts "Finished all calculations."
}

proc popcdensity {ofile} {

	global PC_head_text PO_tail_text

	set	lip_head [atomselect top $PC_head_text]
	set	lip_head_end "_lheadden.dx"

	set	lip_tail [atomselect top $PO_tail_text]
	set	lip_tail_end "_ltailden.dx"

	volmap density $lip_head -allframes -combine avg -res 0.649 -o $ofile$lip_head_end

	volmap density $lip_tail -allframes -combine avg -res 0.649 -o $ofile$lip_tail_end
}

proc dmpcdensity {ofile} {

        global PC_head_text DM_tail_text

        set     lip_head	[atomselect top $PC_head_text]
        set     lip_head_end	"_lheadden.dx"

        set     lip_tail	[atomselect top $DM_tail_text]
        set     lip_tail_end	"_ltailden.dx"

	set	lip_tot		[atomselect top lipids]
	set	lip_tot_end	"_ltotalden.dx"

#	volmap density $lip_head -allframes -combine avg -res 0.649 -o $ofile$lip_head_end

	volmap density $lip_tail -allframes -combine avg -res 0.649 -o $ofile$lip_tail_end

#	volmap density $lip_tot  -allframes -combine avg -res 0.649 -o $ofile$lip_tot_end
}

proc proteindensity {ofile} {

	global	prot_text

	set	protein [atomselect top $prot_text]
	set	prot_end "_protden.dx"

	volmap	density $protein -allframes -combine avg -res 0.649 -o $ofile$prot_end
}

proc autodensity {in} {

	set     infile  [open $in r]

	set     inread  [read -nonewline $infile]

	set     inputs  [split $inread "\n"]

	close	$infile

	## The input file will contain the following: .psf/.pdb, .psf/.dcd, ofile
	##						   0	      1	      2

	set	m	0

	foreach line	$inputs {

		mol new		[lindex $line 0].psf

		mol addfile     [lindex $line 0].pdb

		mol new		[lindex $line 1].psf

		mol addfile     [lindex $line 1].dcd waitfor all

		align   $m [expr $m + 1]

		iondensity	[lindex $line 2]

		mol 	delete	$m
		mol	delete	[expr $m + 1]

		set	m	[expr $m + 2]
	}
}
