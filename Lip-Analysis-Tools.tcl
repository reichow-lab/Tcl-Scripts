proc	Title	{{v ""}} {

	if {$v == "-v"} {
		puts		""
		puts		"Proc Lip-align:	Takes the users selection of any lipids, and aligns them together to a reference lipid.\n"
		puts		"			Press '1' and select the lipids that you want aligned in the VMD GUI...\n"
		puts		"			To run '$ Lip-align <MOLID>'\n"
		puts		"Proc align:		aligns trajectory relative to the protein...i.e. aligns the protein to itself or a separate protein.\n"
		puts		"			Select the molid for the reference and the selection. Uses the first frame from the ref. molid as the reference for the whole system.\n"
		puts		"			To run '$ align'\n"
		puts		"Proc Prot-align:	Takes all of the annular lipids around the protein, and aligns them to a single protein chain.\n"
		puts		"			To run '$ Prot-align <MOLID> <LOGFILE-NAME> <CHAIN>'\n"
		puts            "Proc LipNetwork:	Using the centers file, containing the (x,y) coordinates of the annular\n"
                puts            "			lipid densities (from MD or CryoEM) and calculates the interconnectivity\n" 
                puts            "			of the lipid tails within the densities. This script creates a 'transition\n"
                puts            "			matrix' that shows when any two lipid densities are occupied by a single\n"
                puts            "			lipid. The default minimum IsoLow is 0, meaning the density does not matter\n" 
                puts            "			and thus assumed that a tail is always assigned to a lipid density 'region'.\n"
                puts            "			In order for this program to work, the volume file (.mrc or .dx) must be\n"
                puts            "			loaded into the top molecule (the one containing your .dcd. Furthermore, \n"
                puts            "			you need to run 'align' and 'Prot-align' to ensure the protein/lipids are \n"
                puts            "			at there best-fit to each other.\n"
		puts            "			To run '$ LipNetwork <CENTER FILE> <OUTFILE> <CARBON-THRESHOLD> <ISO-THRESHOLD (default: 0)>'\n"

	} else {

		puts		""
		puts            "			To run '$ Lip-align <MOLID>'\n"
		puts            "			To run '$ align'\n"
		puts            " 			To run '$ Prot-align <MOLID> <LOGFILE-NAME> <CHAIN>'\n"
		puts            " 			To run '$ LipNetwork <CENTER FILE> <OUTFILE> <CARBON-THRESHOLD> <ISO-THRESHOLD (def: none)>'\n"
	}
}


# Aligns lipids together to show how similar/different they are. Demonstrates the range of motion,
# but does not preserve information on their contact with the protein.

proc	Lip-align	{MOLID} {

	set	AtomList	[label list Atoms]

	set	NumFrames	[molinfo $MOLID get numframes]

	set	ref_lip		[atomselect $MOLID "same residue as index [lindex $AtomList {0 0 1}]"]

	set	n		0

	foreach	line	$AtomList {

		set	Index		[lindex $line {0 1}]
		
		set	align_lip	[atomselect top "same residue as (index $Index)"]

		for {set i 0} {$i < $NumFrames} {incr i} {

			animate goto $i

			set trans_matrix	[measure fit $align_lip $ref_lip]

			$align_lip move $trans_matrix

		}

	}
	Title
}


proc align {} {

	puts -nonewline " Provide molID for reference."
	flush stdout

	set ref_molid [gets stdin]

	puts -nonewline " Provide molID for selection."
	flush stdout

	set sel_molid [gets stdin]

	set numframes [molinfo $sel_molid get numframes]

	set ref_frame [atomselect $ref_molid "protein and name CA" frame 0]

	set n 1

	set sys [atomselect $sel_molid all]

	set frame_percent [expr {round($numframes / 30)}]

	for {set i 0} {$i < $numframes} {incr i} {

                animate goto $i

                set align_frame [atomselect $sel_molid "protein and name CA"]

                set trans_matrix [measure fit $align_frame $ref_frame]

                $sys move $trans_matrix

		if {($i % $frame_percent) == 0 } {

			puts -nonewline "*"
			flush stdout
		}

        }

        puts  "\nAlignments complete, ready for RMSD calculations"
	Title
}



proc	Prot-align	{MOLID logfile {RefChain A}} {

	set	log	[open $logfile w]
	puts	$log	"RESID\t\tSEGID\t\tLip-Head RMSF\t\tLip-Tail RMSF\t\tLip-Tot RMSF\n"

	# Select the protein chain (connexin) that will be the reference for transformation

	set ProtChainList	[lsort -unique [[atomselect $MOLID protein] get chain]]	

	set NumFrames		[molinfo $MOLID get numframes]

	set ref_prot		[atomselect $MOLID "protein and chain $RefChain frame 0"]

	foreach chain $ProtChainList {

		set	Phos_List	""

		set	align_prot	[atomselect $MOLID "protein and chain $chain"]

		for {set n 0} {$n < $NumFrames} {incr n} {
			
			animate goto $n

			set	Prot_Lip	[atomselect $MOLID "resname DMPC and same residue as within 10 of (protein and (resid 84 215) and chain $chain)"]

			$Prot_Lip		set beta 1

			set	hold		[[atomselect $MOLID "name P and beta = 1"] get index]

			set	Phos_List	[concat	[lindex $Phos_List] [lindex $hold]]
		
			$Prot_Lip		set beta 0

			unset	hold
		}

		set Phos_Ind	[lsort -unique $Phos_List]

		set Lip_Sel	[atomselect $MOLID "same residue as index [lindex $Phos_Ind]"]

		unset Phos_List

		for {set i 0} {$i < $NumFrames} {incr i} {

			animate goto $i

			set trans_matrix	[measure fit $align_prot $ref_prot]

			$Lip_Sel move $trans_matrix
		}

		foreach phos $Phos_Ind {

			set	IND	[atomselect $MOLID "index $phos"]

			set	ResID	[$IND get resid]
			set	SegID	[$IND get segid]

			set	LipTot	[atomselect $MOLID "resid $ResID and segid $SegID"]
			set	LipTail	[atomselect $MOLID "resid $ResID and segid $SegID and (name C22 to C29 C210 to C214 C32 to C39 C310 to C316)"]
			set	LipHead	[atomselect $MOLID "resid $ResID and segid $SegID and (name O21 O22 O31 O32 O11 to O14 C1 C2 C21 C3 C31 C11 to C15 P N)"]

			set	Headr	[measure rmsf $LipHead]
			set	Tailr	[measure rmsf $LipTail]
			set	Totar	[measure rmsf $LipTot]
			
			set	headhold	0
			set	tailhold	0
			set	totalhold	0

			foreach	head $Headr {set headhold	[expr $headhold + $head]}
			foreach tail $Tailr {set tailhold	[expr $tailhold + $tail]}
			foreach total $Totar {set totalhold	[expr $totalhold + $total]}
	
			set	HeadrA		[expr $headhold / [llength $Headr]]
			set	TailrA		[expr $tailhold / [llength $Tailr]]
			set	TotarA		[expr $totalhold / [llength $Totar]]

		puts	$log	"$ResID\t\t$SegID\t\t$HeadrA\t\t$TailrA\t\t$TotarA\n"
		}

		puts -nonewline "*"

	}
	close	$log

	puts "*"

	Title
}
