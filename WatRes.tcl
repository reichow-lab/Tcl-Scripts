puts "To Run: watres <outname> <cation name>"
proc watres {outname CAT} {
  set numf  [molinfo top get numframes]
  set out   [open $outname-WatRatio.txt w]
  puts $out "Frame\tPoreRatio\tBulkRatio"
  for {set i 0} {$i < $numf} {incr i} {
    animate goto $i
    set pore    [atomselect top "water and same residue as ((abs(z) < 70) and (x^2 + y^2 < 500))"]
    set poreIon [atomselect top "water and same residue as ((abs(z) < 70) and (x^2 + y^2 < 500) and within 3 of (name $CAT))"]
    set bulk    [atomselect top "water and same residue as (abs(z) > 70)"]
    set bulkIon [atomselect top "water and same residue as ((abs(z) > 70) and within 3 of (name $CAT))"]
    set numPore     [$pore num]
    set numPoreIon  [$poreIon num]
    set numBulk     [$bulk num]
    set numBulkIon  [$bulkIon num]
    set PoreRatio   [expr $numPoreIon / $numPore]
    set BulkRatio   [expr $numBulkIon / $numBulk]
    puts $out "$i\t$PoreRatio\t$BulkRatio"
    $pore delete
    $poreIon delete
    $bulk delete
    $bulkIon delete
    unset numPore, numPoreIon, numBulk, numBulkIon, PoreRatio, BulkRatio
  }
}
