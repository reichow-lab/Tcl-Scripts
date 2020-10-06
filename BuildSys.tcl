# This is the build script for the final ion-conduction project. The purpose of
# this script is to construct a gap-junction simulation system inside a membrane
# and solvated (hexagonal box) given a protein .pdb, a membrane .pdb/.psf and a
# system radius.
#
# In the meantime, you must hard-code in the resID of the cysteins participating in disulfide bonds
proc help {} {
  puts "To run this program type: run <protein-prefix> <membrane-prefix> <outname> <isoform(26/46/50)>"
}
help
proc run {prot mem outname iso {rad 70}} {

  mol new $prot.pdb
  # Find all terminal residues for patching
  set termID ""
  set CAlist [[atomselect top "name CA and chain A"] get resid]
  lappend termID [lindex $CAlist 0]
  for {set i 1} {$i < [llength $CAlist]} {incr i} {
    set a [lindex $CAlist $i]
    set b [lindex $CAlist [expr $i + 1]]
    set c [lindex $CAlist [expr $i - 1]]
    if {[expr $b - $a] != [expr $a - $c]} {lappend termID $a}
  }

  # Create protein psf (acetylated n-termini) and center it in the system.
  buildCx [lindex $termID 0] [lindex $termID 1] [lindex $termID 2] [lindex $termID 3] HOLD-temp $iso
  mol delete all
  mol new HOLD-temp.psf
  mol addfile HOLD-temp.pdb
  set prot [atomselect top all]
  set vec [vecinvert [measure center $prot]]
  $prot moveby $vec
  $prot writepdb HOLD-temp.pdb

  # Add membranes to the system
  mol new $mem.psf
  mol addfile $mem.pdb
  set memb [atomselect top "lipids"]
  set vec [vecinvert [measure center $memb]]
  $memb moveby $vec
  $memb moveby {0 0 40}
  $memb writepdb memA.pdb
  $memb writepsf memA.psf
  $memb moveby {0 0 -80}
  set segPlist [[atomselect top "lipids and name P"] get segid]
  set segList [lsort -unique $segPlist]
  foreach seg $segList {
    set sel [atomselect top "all and segid $seg"]
    set newseg [string tolower $seg]
    $sel set segid $newseg
    unset sel newseg
  }
  $memb writepdb memB.pdb
  $memb writepsf memB.psf

  # Merge the membranes together
  merge "memA" "memB" "mem-temp"

  # Merge membranes with protein
  merge "mem-temp" "HOLD-temp" HOLD-TEMP

  # Delete overlapping lipids
  mol delete all
  mol new HOLD-TEMP.psf
  mol addfile HOLD-TEMP.pdb
  set all [atomselect top all]
  set badlip [atomselect top "lipid and same residue as ((within 1.5 of protein) or (x^2 + y^2 < 256))"]
  $all set beta 0
  $badlip set beta 1
  set good [atomselect top "all and beta = 0"]
  $good writepsf HOLD-lipids.psf
  $good writepdb HOLD-lipids.pdb

  # Solvate the system in a waterbox that is ever-so-slightly smaller than that of the lipids
  # Here I am using the coordinates of the phosphates to help define the boundaries of the waterbox
  package require solvate
  mol delete all
  mol new HOLD-lipids.psf
  mol addfile HOLD-lipids.pdb
  set phosP [atomselect top "lipids and name P"]
  set vec [measure minmax $phosP]
  set negZ [vecinvert [list 0 0 [lindex $vec 0 2]]]
  set watvec1 [vecadd [lindex $vec 0] $negZ {0 0 -100}]
  set posZ [vecinvert [list 0 0 [lindex $vec 1 2]]]
  set watvec2 [vecadd [lindex $vec 1] $posZ {0 0 100}]
  solvate HOLD-lipids.psf HOLD-lipids.pdb -o HOLD-RAW -b 1.5 -minmax [list $watvec1 $watvec2]

  # Delete waters that are overlapping in the lipid region
  mol delete all
  mol new HOLD-RAW.psf
  mol addfile HOLD-RAW.pdb
  set phospD [atomselect top "lipids and name P and z < 0"]
  set VecD [measure minmax $phospD]
  set phospU [atomselect top "lipids and name P and z > 0"]
  set VecU [measure minmax $phospU]
  set all [atomselect top all]
  set badwat [atomselect top "water and same residue as (((z > [expr [lindex $VecD 0 2] + 5] and z < [expr [lindex $VecD 1 2] - 5]) or (z > [expr [lindex $VecU 0 2] + 5] and z < [expr [lindex $VecU 1 2] - 5])) and (x^2 + y^2 > 700))"]
  $all set beta 0
  $badwat set beta 1
  set good [atomselect top "all and beta = 0"]
  $good writepsf HOLD-solv.psf
  $good writepdb HOLD-solv.pdb

  # Ionize system, first by splitting into Intracellular & Extracellular compartments (as usual)
  package require autoionize
  mol delete all
  mol new HOLD-solv.psf
  mol addfile HOLD-solv.pdb
  set all [atomselect top all]
  set outwat [atomselect top "water and same residue as (abs(z) < 35 and (x^2 + y^2 > 700))"]
  $all set beta 0
  $outwat set beta 1
  set inwat [atomselect top "all and beta = 0"]
  $outwat writepsf outwat.psf
  $outwat writepdb outwat.pdb
  $inwat writepsf inwat.psf
  $inwat writepdb inwat.pdb
  # Now that they are split, ionize accordingly
  autoionize -psf outwat.psf -pdb outwat.pdb -o oution -seg EON -sc 0.15 -cation SOD -anion CLA -from 5 -between 5
  autoionize -psf inwat.psf -pdb inwat.pdb -o inion -seg ION -sc 0.15 -cation POT -anion CLA -from 5 -between 5

  # Merge ionized systems back together
  merge "oution" "inion" $outname-lwi

  # Cut into a Hexagon
  mol delete all
  mol new $outname-lwi.psf
  mol addfile $outname-lwi.pdb
  set sqrt3 [expr sqrt(3.0)]
  set bound1text "abs(x) < 0.5*$sqrt3*$rad"
  set bound2text "x < $sqrt3*(y+$rad) and x > $sqrt3*(y-$rad)"
  set bound3text "x < $sqrt3*($rad-y) and x > $sqrt3*(-y-$rad)"
  set all [atomselect top all]
  set InBounds [atomselect top "same residue as $bound1text and $bound2text and $bound3text"]
  $InBounds writepdb $outname-lwih.pdb
  $InBounds writepsf $outname-lwih.psf

  # Prepare NAMD input files
  mol delete all
  mol new $outname-lwih.psf
  mol addfile $outname-lwih.pdb
  set all [atomselect top all]
  set prot [atomselect top protein]
  set back [atomselect top "protein and backbone"]
  $all set beta 0
  $prot set beta 1
  $all writepdb $outname-PROTCONST.pdb
  $all set beta 0
  $back set beta 1
  $all writepdb $outname-BACKCONST.pdb
  set pbc [open "$outname-PBC.txt" w]
  set minmax [measure minmax [atomselect top water]]
  set a [expr [lindex $minmax 1 1] - [lindex $minmax 0 1]]
  set b [expr {$a / 2 * sqrt(3.0)}]
  set d [expr {$a / 2}]
  set c [expr [lindex $minmax 1 2] - [lindex $minmax 0 2]]
  puts $pbc "#Curtesy of CHARMM-GUI.."
  puts $pbc "#Hexagonal Periodic Boundary conditions."
  puts $pbc "cellBasisVector1       $a  0.0 0.0"
  puts $pbc "cellBasisVector2       $d  $b  0.0"
  puts $pbc "cellBasisVector3       0.0 0.0 $c"
  puts $pbc "cellOrigin             [measure center [atomselect top protein]]"
  close $pbc

  # Clean up directory
  lappend BADFILES [glob chain*]
  lappend BADFILES [glob mem*]
  lappend BADFILES [glob *wat*]
  lappend BADFILES [glob *ion*]
  lappend BADFILES [glob HOLD*]
  foreach buddy $BADFILES {
    foreach boy $buddy {
      file delete $boy
    }
  }
}
proc merge {pdb1 pdb2 outname} {
  mol delete all
  package require psfgen
  resetpsf
  readpsf $pdb1.psf
  coordpdb $pdb1.pdb
  readpsf $pdb2.psf
  coordpdb $pdb2.pdb
  writepsf $outname.psf
  writepdb $outname.pdb
}
proc buildCx {NT ICC ICN CT outname iso} {
  resetpsf
  if {$iso == "46"} {set disuList [list 54 61 65 189 183 178]
  } elseif {$iso == "50"} {set disuList [list 54 61 65 201 195 190]
  } elseif {$iso == "26"} {set disuList [list 53 60 64 180 174 169]}
  package require psfgen
  set i 1
  set j 2
  set chains [list A B C D E F G H I J K L]
  foreach chain $chains {
    set sel1 [atomselect top "chain $chain and resid $NT to $ICC"]
    set sel2 [atomselect top "chain $chain and resid $ICN to $CT"]
    $sel1 writepdb chain-$chain$i.pdb
    $sel2 writepdb chain-$chain$j.pdb
  }
  topology /media/bassam/Stuff/WORK/Topology\ and\ Parameters/top_all36_prot.rtf
  topology /media/bassam/Stuff/WORK/Topology\ and\ Parameters/toppar_water_ions_namd.str

  pdbalias residue HIS HSD
  set segN [list A1 B1 C1 D1 E1 F1 G1 H1 I1 J1 K1 L1]
  set segC [list A2 B2 C2 D2 E2 F2 G2 H2 I2 J2 K2 L2]
  foreach segn $segN chain $chains {
    resetpsf
    pdbalias residue HIS HSD
    segment $segn {
      first ACE
      last NONE
      pdb chain-$segn.pdb
    }
    coordpdb chain-$segn.pdb $segn
    guesscoord
    writepdb chain-$segn.pdb
    writepsf chain-$segn.psf
  }
  foreach segc $segC {
    resetpsf
    pdbalias residue HIS HSD
    segment $segc {
      first NONE
      last CTER
      pdb chain-$segc.pdb
    }
    coordpdb chain-$segc.pdb $segc
    guesscoord
    writepdb chain-$segc.pdb
    writepsf chain-$segc.psf
  }
  resetpsf
  foreach segn $segN segc $segC {
    readpsf chain-$segn.psf
    coordpdb chain-$segn.pdb
    readpsf chain-$segc.psf
    coordpdb chain-$segc.pdb
  }
  foreach segn $segN segc $segC {
    patch DISU $segn:[lindex $disuList 0] $segc:[lindex $disuList 3]
    patch DISU $segn:[lindex $disuList 1] $segc:[lindex $disuList 4]
    patch DISU $segn:[lindex $disuList 2] $segc:[lindex $disuList 5]
  }
  writepsf $outname.psf
  writepdb $outname.pdb
}
