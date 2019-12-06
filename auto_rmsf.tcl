puts "type the following to run program: 'run outfile_name'"

proc align {rmolid smolid} {

        set ref_molid $rmolid

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

}
proc run {ofile resi resf} {

	set initframe 0

	set finaframe [expr [molinfo top get numframes] - 1]

	set first $resi

	set final $resf

	set output [open $ofile w]

	set rmsf_res [atomselect top "protein and resid $first to $final and name CA and chain A"]

	set res_id [$rmsf_res get resid]

	foreach RES $res_id {

		set CA [atomselect top "resid $RES and name CA"]

		set rmsf_calc [measure rmsf $CA]

		puts $rmsf_calc

		puts $output "$RES\t$rmsf_calc"

		puts "rmsf of residue $RES out of $final calculated"

	}

	close	$output

}
proc	 autormsf    {in} {

         set     infile  [open $in r]

         set     inread  [read -nonewline $infile]

         set     inputs  [split $inread "\n"]

         close   $infile

         ## The input file will contain the CryoEM .psf/.pdb, .psf/.dcd, OUT, ResI, ResF
         ##                                            0          1       2    3     4      
         set     m       0

         foreach line    $inputs {

                 mol new         [lindex $line 0].psf

                 mol addfile     [lindex $line 0].pdb

                 mol new         [lindex $line 1].psf

                 mol addfile     [lindex $line 1].dcd waitfor all

                 align   $m [expr $m + 1]

                 run     [lindex $line 2] [lindex $line 3] [expr $m + 4]

                 mol     delete  $m
                 mol     delete  [expr $m + 1]

                 set     m       [expr $m + 2]
         }
 }

