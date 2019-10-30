#puts "for initial alignment (if needed) run: 'align'"
#puts "type the following to run program: 'run name_of_ofile'"


proc align {rmolid smolid} {

#        puts -nonewline " Provide molID for reference."
#        flush stdout

        set ref_molid $rmolid

#        puts -nonewline " Provide molID for selection."
#        flush stdout

        set sel_molid $smolid

        set numframes [molinfo $sel_molid get numframes]

        set ref_frame [atomselect $ref_molid "protein and name CA" frame 0]

        set n 1

        set sys [atomselect $sel_molid all]

        for {set i 0} {$i < $numframes} {incr i} {

                animate goto $i

                set align_frame [atomselect $sel_molid "protein and name CA"]

                set trans_matrix [measure fit $align_frame $ref_frame]

                $sys move $trans_matrix

                if {($n % 100) == 0 } {

                        puts "alignment $n of $numframes"
                }

                incr n

        }

        puts "Alignments complete, ready for RMSD calculations"
}

proc run {ofile IonName id} {

	## Selection of Atom for trajectory output
  
#	puts -nonewline "Select ion by entering its name (potassium = POT, sodium = SOD etc...).  "
	
#	flush stdout

	set ion_name $IonName

#	puts -nonewline "What is the MolID? "

#	flush stdout

	set molid $id

	set ION [atomselect $molid "name $ion_name and not ((abs(z) < 40) and (x^2 + y^2 > 350))"]

	set ind [$ION get index]

	set num_ion [$ION num]

#	puts "There are $num_ion ions"

	set prot [atomselect $molid protein]

	set n 0

	foreach IND $ind {

		set filename [concat $ofile$n]

		incr n

		set ion [atomselect $molid "index $IND"]

## Get number of frames loaded into top molecule

		set numframes [molinfo $molid get numframes]

## Open output-file

		set output [open $filename w]

## Loop over frames

		for {set f 0} {$f < $numframes} {incr f} {

			animate goto $f

			set protcom	[measure center $prot]

			set protz	[lindex $protcom 2]

			set ioncom	[measure center $ion]

			set ionz	[lindex $ioncom 2]

			set pos		[expr $ionz - $protz]

			puts $output "$f\t$pos"

		}

		close	$output

#		puts "$n of $num_ion ion trajectories calculated"

	}
}
proc	ionz	{in} {

	set	infile	[open $in r]
	
	set	inread	[read -nonewline $infile]
	
	set	inputs	[split $inread "\n"]
	
	close	$infile
	
	## The input file will contain the CryoEM .psf/.pdb, .psf/.dcd, CatOUT, AnOUT, CatIon, AnIon
	##					      0          1         2      3       4      5
	set	m	0
	
	foreach	line	$inputs {

		mol new		[lindex $line 0].psf
		
		mol addfile	[lindex $line 0].pdb

		mol new		[lindex $line 1].psf

		mol addfile	[lindex $line 1].dcd waitfor all

		align	$m [expr $m + 1]

		run 	[lindex $line 2] [lindex $line 4] [expr $m + 1]

		run	[lindex $line 3] [lindex $line 5] [expr $m + 1]

		mol	delete	$m
		mol	delete	[expr $m + 1]

		set	m	[expr $m + 2]
	}
}
