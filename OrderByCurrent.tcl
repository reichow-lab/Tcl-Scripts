puts "create a list of dcd's that you want processed."
puts "create a list of the min-mid & mid-max boundaries"
proc orderbycurrent {TextIN PSFin DCDinList peaklist outname st} {
  set INFILE [open $TextIN r]
  set INDATA [read $INFILE]
  set last [expr [llength $peaklist] - 1]
  set peakFROM [lreplace $peaklist $last $last]
  set peakTO  [lreplace $peaklist 0 0]
  close $INFILE
  set DATA [split $INDATA "\n"]
  # peakFROM and peakTO are the boundaries between "current-states", and allow
  # for an arbitrary number of "current-states" than simply min mid & max
  # Generate list of Bins. BINS is initialized with the min and max bins, middle
  # bins are added in the next loop
  set BINS "{} {}"
  set BinN ""
  for {set i 1} {$i <= [expr [llength $peaklist] - 1]} {incr i} {
    lappend BinN $i
    lappend BINS ""
  }
  puts $BinN
  puts $BINS
  puts $peakFROM
  puts $peakTO
  # generate list of frames and their assignments
  foreach from $peakFROM to $peakTO bin $BinN {
    foreach line $DATA {
      if {[string is double -strict [lindex $line 1]]} {
        puts "debug-0"
        if {[lindex $line 1] <= [expr [lindex $peakFROM 0]]} {
          puts "debug-0.5"
          lset BINS 0 end+1 [expr int([lindex $line 0])]
          puts "debug-1"
        } elseif {([lindex $line 1] > $from) && ([lindex $line 1] <= $to)} {
          lset BINS $bin end+1 [expr int([lindex $line 0])]
          puts "debug-2"
        } elseif {[lindex $line 1] > [lindex $peakTO [expr [llength $peakTO] - 1]} {
          lset BINS [expr [llength $BINS] - 1] end+1 [expr int([lindex $line 0])]
          puts "debug-3"
        }
      }
    }
  }
  puts "debug-4"
  # sort bin lists to remove redundant frame-IDs
  for {set i 0} {$i < [llength $BINS]} {incr i} {
    lappend BINSort [lreplace $BINS 0 [expr [llength $BINS] - 1] [lsort -unique [lindex $BINS $i]]]
  }
  # Now that the frames are ordered, properly order the dcd
  set b 0
  foreach Bin $BINSort {
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
