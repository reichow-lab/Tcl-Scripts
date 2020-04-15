#	Program:	HOLE-APBS_gen.tcl
#	Author:		Bassam Haddad
#
#	This program prepares the appropriate input files for both BOLE and PDB2PQR/APBS. The inputs to this program are trajectories for any gap junction
#	with or without system components (i.e. water, lipids, ions). A reference pdb should be provided to ensure the geometric center of the system is
#	the protein pore. This also ensures the pore axis is along the z-axis in case the protein drifts throughout the simulation.
#
#	To best use this program, utilize the auto_holegen proc{} by feeding an input file with paths to the inputs.

puts "To use this program:\n\t\tholegen <outfile> <PDBOutName> <RefMolID>\n\t\tauto_holegen <InputFile>"

proc holegen {ofile PDB RefMolID} {

	# Center reference on channel

	set ref_all	[atomselect $RefMolID all]
	set ref_prot	[atomselect $RefMolID protein]

	set prot_cen [measure center $ref_prot]

	set x   [expr [lindex $prot_cen 0] * -1]
	set y   [expr [lindex $prot_cen 1] * -1]
	set z   [expr [lindex $prot_cen 2] * -1]
	set vec [list $x $y $z]
	$ref_all moveby $vec

	set numframes [molinfo top get numframes]

	# Align Protein (align system to the centered reference

	set ref_frame [atomselect $RefMolID "protein and name CA"]

	for {set i 0} {$i < $numframes} {incr i} {

                animate goto $i

                set align_frame [atomselect top "protein and name CA"]

                set trans_matrix [measure fit $align_frame $ref_frame]

                $all move $trans_matrix

        }

	set n 0

	for {set i 0} {$i < $numframes} {incr i} {

		set HOLEinp "$ofile-HOLE-$n.in"
		set APBSinp "$ofile-APBS-$n.in"

		set pdbfilename "$PDB$n.pdb"

		incr n

		animate goto $i

		$prot writepdb $pdbfilename

		set HOLEout [open $HOLEinp w]
		set APBSout [open $APBSinp w]

		puts	$HOLEout	"coord\t$pdbfilename"
		puts	$HOLEout	"radius\t/home/bassam/hole2/rad/simple.rad"
		puts	$HOLEout	"cvect\t0 0 1"
		puts	$HOLEout	"cpoint\t0 0 0"
		puts	$HOLEout	"sample\t1.0"
		puts	$HOLEout	"ENDRAD\t14"
		close	$HOLEout

		puts	$APBSout	"read\n\tmol pqr $PDB$n.pqr\nend\nelec\n\tmg-auto\n\tdime 225 193 289\n\tcglen 149.6544 141.2020 217.5235\n\tfglen 108.0320 103.0600 147.9550\t\ncgcent mol 1\n\tfgcent mol 1\n\tmol 1\n\tlpbe\n\tbcfl sdh\n\tpdie 2.0000\n\tsdie 80.000\n\tsrfm smol\n\tchgm spl2\n\tsdens 10.00\n\tsrad 1.40\n\tswin 0.30\n\ttemp 310\n\tion 1 0.150 2.0\n\tion -1 0.150 2.0\n\tcalcenergy no\n\tcalcforce no\n\twrite pot dx $PDB$n.dx\nend\nquit"
		close	$APBSout
	}
}

proc auto_holegen {in} {
	
	set     infile  [open $in r]
	set     inread  [read -nonewline $infile]
	set     inputs  [split $inread "\n"]
	close   $infile 

	## The input file will contain:		  Ref-PDB	  PSF/DCD	OutFileName	PDBOutName	
	##					     0		     1		     2		    3
	
	set     m       0

	foreach line    $inputs {

		mol new     	[lindex $line 0].pdb
		mol new         [lindex $line 1].psf
		mol addfile     [lindex $line 1].dcd waitfor all
		
		holegen		[lindex $line 2] [lindex $line 3] $m

		set	m	[expr $m + 2]
}
