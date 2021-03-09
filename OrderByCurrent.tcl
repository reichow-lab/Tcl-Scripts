puts "create a list of dcd's that you want processed."

proc orderbycurrent {TextIN PSFin DCDinList outname st} {
  set INFILE [open $TextIN r]
  set INDATA [read $INFILE]
  close $INFILE
  set DATA [split $INDATA "\n"]
  # read through the data and find the maximum Current
  set max 0
  foreach line $DATA {
    if {[string is double -strict [lindex $line 1]]} {
      if {[expr abs([lindex $line 1])] > $max} {set max [expr abs([lindex $line 1])]}
    }
  }
  # simply split into 3 groups: min mid max
  set min [expr $max / 3.0]
  set mid [expr $max * (2.0 / 3.0)]
  puts "$min, $mid, $max"
  # generate list of frames and their assignments
  set Bin0 {}
  set Bin1 {}
  set Bin2 {}
  foreach line $DATA {
    if {[string is double -strict [lindex $line 1]]} {
      if {[expr abs([lindex $line 1])] <= $min} {
        lappend Bin0 [expr int([lindex $line 0])]
      } elseif {([expr abs([lindex $line 1])] > $min) && ([expr abs([lindex $line 1])] <= $mid)} {
        lappend Bin1 [expr int([lindex $line 0])]
      } elseif {[expr abs([lindex $line 1])] > $mid} {
        lappend Bin2 [expr int([lindex $line 0])]
      }
    }
  }
  puts "$Bin1"
  # rearrange the order of the assigned frame_list
  set BINS [list $Bin0 $Bin1 $Bin2]
  set b 0
  set ordered {}
  foreach bin $BINS {
    foreach line $bin {
      if {[lindex $line 2] == $b} {
        lappend ordered [lindex $line 0]
      }
    }
    incr b
  }
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
      if {[lsearch $Bin $i] >= 0} {puts "save"
  		} else {animate delete beg $i end $i top}
    }
    set all [atomselect top all]
    animate write dcd $outname.Bin$b.dcd beg 0 end -1 waitfor all sel $all top
    mol delete top
    incr b
  }
}
