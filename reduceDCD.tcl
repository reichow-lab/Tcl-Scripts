# Takes list of .dcds, associated .psf, their current frame-step, and the desired frame-step.
puts "create a list of dcd's that you want processed."
proc reduce {DCDlist PSFin InitFrameStep FinalFrameStep} {
  # calculate step-size needed to get to 10ps/frame
  set stride [expr $FinalFrameStep/$InitFrameStep]
  set i 0
  foreach dcd $DCDlist {
    mol new $PSFin
    mol addfile $dcd step $stride waitfor all
    animate write dcd $dcd.10ps.dcd beg 0 end -1 waitfor all sel [atomselect top all] $i
    incr i
    mol delete top
    file delete $dcd
  }
}
