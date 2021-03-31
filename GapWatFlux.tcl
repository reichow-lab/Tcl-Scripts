# Program: GapWatFlux.tcl
# Author: Bassam Haddad
#
# This program calculates the water flux through a gap junction pore at the constriction region (abs(z)=45), the parahelix (abs(z)=30), and the center (z=0)

puts "Ensure the trajectories are aligned before the calculation."
puts "To run: gapwatflux outfile dt (100ps = 0.1ns)"

proc gapwatflux {outfile dt} {
  set numframes [molinfo top get numframes]
  set Zbounds [list 45 30 0 -30 -45]
  set Zlists  [list Z1 Z2 Z3 Z4 Z5]
  #loop through boundaries and run the calculation
  foreach zl $Zlists zb $Zbounds {
    animate goto 0
    # reset all of the waters to beta = 0
    set AllWater [atomselect top water]
    $AllWater set beta 0
    # make initial 'chamber' assignments
    set C1 [atomselect top "name OH2 and ((z < 60 and z > $zb) and (x^2 + y^2 < 400))"]
    set C2 [atomselect top "name OH2 and ((z > -60 and z < $zb) and (x^2 + y^2 < 400))"]
    $C1 set beta 1
    $C2 set beta 2
    unset C1 C2
    # loop through frames and perform analysis
    for {set i 1} {$i < $numframes} {incr i} {
      animate goto $i
      set UpWat [atomselect top "name OH2 and ((z < 60 and z > $zb) and (x^2 + y^2 < 400)) and beta = 2"]
      set DnWat [atomselect top "name OH2 and ((z > -60 and z < $zb) and (x^2 + y^2 < 400)) and beta = 1"]
      set UpFlux  [$UpWat num]
      set DnFlux  [$DnWat num]
      set NetFlux [expr $UpFlux - $DnFlux]
      lappend $zl $NetFlux
      set C1 [atomselect top "name OH2 and ((z < 60 and z > $zb) and (x^2 + y^2 < 400))"]
      set C2 [atomselect top "name OH2 and ((z > -60 and z < $zb) and (x^2 + y^2 < 400))"]
      $C1 set beta 1
      $C2 set beta 2
      unset UpWat DnWat C1 C2 UpFlux DnFlux NetFlux
    }
  }
  # Write out the fluxes
  set out [open $outfile.txt w]
  set i 1
  puts $out "Time(ns)\tz = 45\tz = 30\tz = 0\tz = -30\tz = -45"
  foreach z1 $Z1 z2 $Z2 z3 $Z3 z4 $Z4 z5 $Z5 {
    set t [expr $i * $dt]
    puts $out "$t\t$z1\t$z2\t$z3\t$z4\t$z5"
    incr i
  }
}
