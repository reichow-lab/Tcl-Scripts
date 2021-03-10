puts "create a list of dcd's that you want processed."
puts "create a list of the min-mid & mid-max boundaries"
proc orderbycurrent {TextIN PSFin DCDinList peaklist outname st} {
  set INFILE [open $TextIN r]
  set INDATA [read $INFILE]
  close $INFILE
  set DATA [split $INDATA "\n"]
  # read through the data and find the maximum Current
  set max 0
  foreach line $DATA {
    if {[string is double -strict [lindex $line 1]]} {
      if {[expr abs([lindex $line 1])] > [expr abs($max)]} {set max [lindex $line 1]}
    }
  }
  # simply split into 3 groups: min mid max
  set minmid [lindex $peaklist 0]
  set midmax [lindex $peaklist 1]
  # generate list of frames and their assignments
  set Bin0 {}
  set Bin1 {}
  set Bin2 {}
  foreach line $DATA {
    if {[string is double -strict [lindex $line 1]]} {
      if {[lindex $line 1] <= [expr $minmid]} {
        lappend Bin0 [expr int([lindex $line 0])]
      } elseif {([lindex $line 1]] > $minmid) && ([lindex $line 1] <= $midmax)} {
        lappend Bin1 [expr int([lindex $line 0])]
      } elseif {[lindex $line 1] > $midmax} {
        lappend Bin2 [expr int([lindex $line 0])]
      }
    }
  }
  # rearrange the order of the assigned frame_list
  set BINS [list $Bin0 $Bin1 $Bin2]
  # Now that the frames are ordered, properly order the dcd
  set b 0
  foreach Bin $BINS {
    puts "$Bin"
    mol new $PSFin
    foreach dcd $DCDinList {
      mol addfile $dcd step $st waitfor all
    }
    set numframes [molinfo top get numframes]
    for {set i [expr $numframes -1]} {$i >= 0} {incr i -1} {
      animate goto $i
      if {[lsearch $Bin $i] >= 0} {continue
  		} else {animate delete beg $i end $i top}
    }
    set all [atomselect top all]
    animate write dcd $outname.Bin$b.dcd beg 0 end -1 waitfor all sel $all top
    mol delete top
    incr b
  }
}
