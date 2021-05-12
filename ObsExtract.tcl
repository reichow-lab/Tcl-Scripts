source ~/Scripts/TCL/Tcl-Scripts/auto-ionz-IN.tcl
proc average {list} {
	    expr {[tcl::mathop::+ {*}$list 0.0] / max(1, [llength $list])}
    }
proc run {ofile id {from "CZ"}} {
	set outU [open $ofile.U.obs w]
	set outL [open $ofile.L.obs w]
	set Uchains [list A B C D E F]
	set Lchains [list G H I J K L]
	set UJchains [list B C D E F A]
	set LJchains [list H I J K L G]
	set DiHe [list N CA C O]
	align $id $id
	set numf [molinfo top get numframes]
	set counter 1
	foreach chainI $Uchains chainJ $UJchains {
		puts $outU "Chain:\t$chainI"
		set DiHeIn {}
		foreach name $DiHe {
			set sel [atomselect top "protein and chain $chainI and resid 47 and name $name"]
			lappend DiHeIn [$sel get index]
			$sel delete
		}
		set sel1 [atomselect top "protein and chain $chainI and ((resid 3 and name CG) or (resid 5 and name HG1))"]
		set D3S5In [$sel1 get index]
		set sel2 [atomselect top "(protein and chain $chainI and resid 9 and name $from) or (protein and chain $chainJ and resid 12 and name CG)"]
		set RN9E12In [$sel2 get index]
		$sel1 delete
		$sel2 delete
		for {set i 0} {$i < $numf} {incr i} {
			animate goto $i
			set dihed [measure dihed $DiHeIn]
			set D3S5 [measure bond $D3S5In]
			set RN9E12 [measure bond $RN9E12In]
			puts $outU "$i\t$dihed\t$D3S5\t$RN9E12"
		}
	}
	close $outU
	set counter 1
	foreach chainI $Lchains chainJ {
		puts $outL "Chain:\t$chainI"
		set DiHeIn {}
		foreach name $DiHe {
			set sel [atomselect top "protein and chain $chain and resid 47 and name $name"]
			lappend DiHeIn [$sel get index]
			$sel delete
		}
		set sel1 [atomselect top "protein and chain $chainI and ((resid 3 and name CG) or (resid 5 and name HG1))"]
		set D3S5In [$sel1 get index]
		set sel2 [atomselect top "(protein and chain $chainI and resid 9 and name $from) or (protein and chain $chainJ and resid 12 and name CG)"]
		set RN9E12In [$sel2 get index]
		$sel1 delete
		$sel2 delete
		for {set i 0} {$i < $numf} {incr i} {
			animate goto $i
			set dihed [measure dihed $DiHeIn]
			set D3S5 [measure bond $D3S5In]
			set RN9E12 [measure bond $RN9E12In]
			puts $outL "$i\t$dihed\t$D3S5\t$RN9E12"
		}
	}
	close $outL
}
