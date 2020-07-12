proc align {RefMolID SelMolID} {

	# Center reference on channel

	set ref_all	[atomselect $RefMolID all]
	set ref_prot	[atomselect $RefMolID protein]

	set prot_cen [measure center $ref_prot]

	set x   [expr [lindex $prot_cen 0] * -1]
	set y   [expr [lindex $prot_cen 1] * -1]
	set z   [expr [lindex $prot_cen 2] * -1]
	set vec [list $x $y $z]
	$ref_all moveby $vec

	set numframes [molinfo $SelMolID get numframes]

	# Align Protein (align system to the centered reference)

	set ref_frame [atomselect $RefMolID "protein and name CA"]

	for {set i 0} {$i < $numframes} {incr i} {

                animate goto $i

                set align_frame [atomselect $SelMolID "protein and name CA"]
		set all	[atomselect $SelMolID all]
                set trans_matrix [measure fit $align_frame $ref_frame]

                $all move $trans_matrix

        }
