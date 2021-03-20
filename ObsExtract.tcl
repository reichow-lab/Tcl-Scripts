source ~/Scripts/TCL/Tcl-Scripts/auto-ionz-IN.tcl
set Uchains [list A B C D E F] 
set Lchains [list G H I J K L]
set DiHe [list N CA C O]
proc average {list} {
	    expr {[tcl::mathop::+ {*}$list 0.0] / max(1, [llength $list])}
    }
proc run {ofile id} {
	set outU [open $ofile.U.obs w]
	set outL [open $ofile.L.obs w]
	set Uchains [list A B C D E F]
	set Lchains [list G H I J K L]
	set DiHe [list N CA C O]
	align $id $id
	set numf [molinfo top get numframes]
	foreach chain $Uchains {
		puts $outU "Chain:\t$chain"
		set DiHeIn {}
		foreach name $DiHe {
			set sel [atomselect top "protein and chain $chain and resid 47 and name $name"]
			lappend DiHeIn [$sel get index]
			$sel delete
		}
		set sel [atomselect top "protein and chain $chain and ((resid 3 and name CG) or (resid 5 and name HG1))"]
		set D3S5In [$sel get index]
		$sel delete
		for {set i 0} {$i < $numf} {incr i} {
			animate goto $i
			set dihed [measure dihed $DiHeIn]
			set D3S5 [measure bond $D3S5In]
			puts $outU "$i\t$dihed\t$D3S5"
		}	
	}
	close $outU
	foreach chain $Lchains {
		puts $outL "Chain:\t$chain"
		set DiHeIn {}
		foreach name $DiHe {
			set sel [atomselect top "protein and chain $chain and resid 47 and name $name"]
			lappend DiHeIn [$sel get index]
			$sel delete
		}
		set sel [atomselect top "protein and chain $chain and ((resid 3 and name CG) or (resid 5 and name HG1))"]
		set D3S5In [$sel get index]
		$sel delete
		for {set i 0} {$i < $numf} {incr i} {
			animate goto $i
			set dihed [measure dihed $DiHeIn]
			set D3S5 [measure bond $D3S5In]
			puts $outL "$i\t$dihed\t$D3S5"
		}
	}
	close $outL
}
