#puts "for initial alignment (if needed) run: 'align'"
#puts "type the following to run program: 'run name_of_ofile'"


proc align {rmolid smolid} {
        set ref_molid $rmolid
        set sel_molid $smolid
        set numframes [molinfo $sel_molid get numframes]
	animate goto 0
	set prot [atomselect $ref_molid "protein and name CA"]
	set sys [atomselect $sel_molid all]
	$sys moveby [vecinvert [measure center $prot]]
        set ref_frame [atomselect $ref_molid "protein and name CA" frame 0]
        for {set i 0} {$i < $numframes} {incr i} {
                animate goto $i
                set align_frame [atomselect $sel_molid "protein and name CA"]
                set trans_matrix [measure fit $align_frame $ref_frame]
		$sys move $trans_matrix
        }
}
proc run {ofile IonName id} {
	set ion_name $IonName
	set molid $id
  if {$IonName == "water"} {
    set ION [atomselect $molid "water and name OH2"]
  } else {set ION [atomselect $molid "name $ion_name and not ((abs(z) < 40) and (x^2 + y^2 > 500))"]}
	set ind [$ION get index]
	set num_ion [$ION num]
	set prot [atomselect $molid protein]
	set numframes [molinfo $molid get numframes]
## Open output-file
	set output [open $ofile.txt w]
	foreach IND $ind {
    		puts $output "IonID: $IND"
		set ion [atomselect $molid "index $IND"]
## Loop over frames
		for {set f 0} {$f < $numframes} {incr f} {
			animate goto $f
			set protcom	[measure center $prot]
      set protx [lindex $protcom 0]
      set proty [lindex $protcom 1]
			set protz	[lindex $protcom 2]
			set ioncom	[measure center $ion]
      set ionx  [lindex $ioncom 0]
      set iony  [lindex $ioncom 1]
			set ionz	[lindex $ioncom 2]
      set posx  [expr $ionx - $protx]
      set posy  [expr $iony - $proty]
			set posz	[expr $ionz - $protz]
			puts $output "$f\t$posx\t$posy\t$posz"
		}
	}
  	close	$output
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
		mol	delete	all
		# mol	delete	$m
		# mol	delete	[expr $m + 1]
		set	m	[expr $m + 2]
	}
}
