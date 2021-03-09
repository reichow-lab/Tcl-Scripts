puts "create a list of dcd's that you want processed."

proc orderbycurrent {TextIN PSFin DCDinList outname st} {
  set INFILE [open $TextIN r]
  set INDATA [read $INFILE]
  close $INFILE
  set DATA [split $INDATA "\n"]
  # read through the data and find the maximum Current
  set max 0
  foreach line $DATA {
    puts "[lindex $line 1]"
    if {[expr abs([lindex $line 1])] >= $max} {set max [expr abs([lindex $line 1])]}
    puts "$max"
  }
  puts "debug-1"
  # simply split into 3 groups: min mid max
  set min [expr $max / 3.0]
  set mid [expr $max * (2.0 / 3.0)]
  # generate list of frames and their assignments
  set Bin0 {}
  set Bin1 {}
  set Bin2 {}
  foreach line $DATA {
    if {[expr abs([lindex $line 1])] <= $min} {
      lappend Bin0 [list [lindex $line 0] [lindex $line 1] 0]
    } elseif {([expr abs([lindex $line 1])] > $min) && ([expr abs([lindex $line 1])] <= $mid)} {
      lappend Bin1 [list [lindex $line 0] [lindex $line 1] 1]
    } elseif {[expr abs([lindex $line 1])] > $mid} {
      lappend Bin2 [list [lindex $line 0] [lindex $line 1] 2]
    }
  }
  puts "debug-2"
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
  puts "debug-3"
  # Now that the frames are ordered, properly order the dcd
  set b 0
  foreach Bin $BINS {
    mol new $PSFin.psf
    foreach dcd $DCDinList {
      mol addfile $dcd step $st waitfor all
    }
    set numframes [molinfo top get numframes]
    for {set i [expr $numframes -1]} {$i >= 0} {incr i -1} {
      animate goto $i
      if {[lsearch -exact $Bin $i] >= 0} {continue
  		} else {animate delete beg $i end $i top}
    }
    set all [atomselect top all]
    animate write dcd $outname.Bin$b.dcd beg 0 end -1 waitfor all sel $all top
    mol delete top
    incr b
  }
}
