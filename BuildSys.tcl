# This is the build script for the final ion-conduction project. The purpose of
# this script is to construct a gap-junction simulation system inside a membrane
# and solvated (hexagonal box) given a protein .pdb, a membrane .pdb/.psf and a
# system radius.
#
# In the meantime, you must hard-code in the resID of the cysteins participating in disulfide bonds
proc help {} {
  puts "To run this program type: build <protein-prefix> <membrane-prefix> -o <outname> -iso <26/\[46\]/50> -ion <ion name> -wat <\[0\]/1> -rad <minimum radius> -hex <0/\[1\]> -mut <\[0\]/1>"
}
help
proc build {prot mem args}  {
  # Set named arguments
  array set opt [concat {-o "TEST" -iso "46" -ion "POT" -wat 0 -rad 70 -hex 1 -mut 0} $args]
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
  buildCx [lindex $termID 0] [lindex $termID 1] [lindex $termID 2] [lindex $termID 3] HOLD-temp $opt(-iso) $opt(-wat) $opt(-mut)
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
  set badwat [atomselect top "water and same residue as (((z > [expr [lindex $VecD 0 2] + 5] and z < [expr [lindex $VecD 1 2] - 5]) or (z > [expr [lindex $VecU 0 2] + 5] and z < [expr [lindex $VecU 1 2] - 5])) and (x^2 + y^2 > 740))"]
  $all set beta 0
  $badwat set beta 1
  if {$opt(-wat) == 1} {
    set goodwat [atomselect top "water and segname SW"]
    $goodwat set beta 0
    }
  set good [atomselect top "all and beta = 0"]
  $good writepsf HOLD-solv.psf
  $good writepdb HOLD-solv.pdb

# Cut into a Hexagon
  if {$opt(-hex) == 1} {
    mol delete all
    mol new HOLD-solv.psf
    mol addfile HOLD-solv.pdb
    set sqrt3 [expr sqrt(3.0)]
    set bound1text "abs(x) < 0.5*$sqrt3*$opt(-rad)"
    set bound2text "x < $sqrt3*(y+$opt(-rad)) and x > $sqrt3*(y-$opt(-rad))"
    set bound3text "x < $sqrt3*($opt(-rad)-y) and x > $sqrt3*(-y-$opt(-rad))"
    set all [atomselect top all]
    set InBounds [atomselect top "same residue as $bound1text and $bound2text and $bound3text"]
    $InBounds writepdb $opt(-o)-lwh.pdb
    $InBounds writepsf $opt(-o)-lwh.psf
  } elseif {$opt(-hex) == 0} {
    $good writepsf $opt(-o)-lwh.psf
    $good writepdb $opt(-o)-lwh.pdb
  }

  # Ionize system, first by splitting into Intracellular & Extracellular compartments (as usual)
  package require autoionize
  mol delete all
  mol new $opt(-o)-lwh.psf
  mol addfile $opt(-o)-lwh.pdb
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
  autoionize -psf inwat.psf -pdb inwat.pdb -o inion -seg ION -sc 0.15 -cation $opt(-ion) -anion CLA -from 5 -between 5

  # Merge ionized systems back together
  merge "oution" "inion" $opt(-o)-lwih

  # Prepare NAMD input files
  mol delete all
  mol new $opt(-o)-lwih.psf
  mol addfile $opt(-o)-lwih.pdb
  set all [atomselect top all]
  set prot [atomselect top protein]
  set back [atomselect top "protein and backbone"]
  set lip [atomselect top "lipids"]
  $all set beta 1
  $lip set beta 0
  $all writepdb $opt(-o)-LIPMELT.pdb
  $all set beta 0
  $prot set beta 1
  $all writepdb $opt(-o)-PROTCONST.pdb
  $all set beta 0
  $back set beta 1
  $all writepdb $opt(-o)-BACKCONST.pdb
  set pbc [open "$opt(-o)-PBC.txt" w]
  set a [expr $opt(-rad) * 2]
  set b [expr {$a / 2 * sqrt(3.0)}]
  set d [expr {$a / 2}]
  set c 200
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
proc buildCx {NT ICC ICN CT outname iso strucwat mut} {
  if {$mut == 1} {
    puts "How many mutations would you like to make? "
    flush stdout
    set mutnum [gets stdin]
    for {set m 0} {$m < $mutnum} {incr m} {
      puts "RESID? "
      flush stdout
      lappend mres [gets stdin]
      puts "segn or segc (type 'n' or 'c')"
      flush stdout
      lappend mseg [gets stdin]
      puts "New AA? "
      flush stdout
      lappend mcod [gets stdin]
    }
  }
  if {$iso == "46"} {set disuList [list 54 61 65 189 183 178]
  } elseif {$iso == "50"} {set disuList [list 54 61 65 201 195 190]
  } elseif {$iso == "26"} {set disuList [list 53 60 64 180 174 169]}
  package require psfgen
  resetpsf
  set i 1
  set j 2
  set chains [list A B C D E F G H I J K L]
  foreach chain $chains {
    set sel1 [atomselect top "chain $chain and resid $NT to $ICC"]
    set sel2 [atomselect top "chain $chain and resid $ICN to $CT"]
    $sel1 writepdb chain-$chain$i.pdb
    $sel2 writepdb chain-$chain$j.pdb
  }

  topology /home/bassam/Topology\ and\ Parameters/top_all36_prot.rtf
  topology /home/bassam/Topology\ and\ Parameters/toppar_water_ions_namd.str

  pdbalias residue HIS HSD
  set segN [list A1 B1 C1 D1 E1 F1 G1 H1 I1 J1 K1 L1]
  set segC [list A2 B2 C2 D2 E2 F2 G2 H2 I2 J2 K2 L2]
  foreach segn $segN chain $chains {
    resetpsf
    pdbalias residue HIS HSD
    segment $segn {
      first ACE
      last NONE
      if {[llen $mutnum] > 0} {
        foreach seg $mseg res $mres code $mcod {
          if {$seg == 'n'} {
          mutate $res $code
          }
        }
      }
      pdb chain-$segn.pdb
      #if {[llen $mutnum] > 0} {
      #  foreach seg $mseg res $mres code $mcod {
      #    if {$seg == 'n'} {
      #    mutate $res $code
      #    }
      #  }
      #}
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
      #if {[llen $mutnum] > 0} {
      #  foreach seg $mseg res $mres code $mcod {
      #    if {$seg == 'c'} {
      #      mutate $res $code
      #    }
      #  }
      #}
    }
    coordpdb chain-$segc.pdb $segc
    guesscoord
    writepdb chain-$segc.pdb
    writepsf chain-$segc.psf
  }
  if {$strucwat == 1} {
    set wat [atomselect top water]
    $wat writepdb EMwat.pdb
    resetpsf
    pdbalias residue HOH TIP3
    segment SW {
      first NONE
      last NONE
      pdb EMwat.pdb
    }
    pdbalias atom HOH O OH2
    coordpdb EMwat.pdb SW
    guesscoord
    writepdb strucwat.pdb
    writepsf strucwat.psf
  }
  resetpsf
  foreach segn $segN segc $segC {
    readpsf chain-$segn.psf
    coordpdb chain-$segn.pdb
    readpsf chain-$segc.psf
    coordpdb chain-$segc.pdb
  }
  if {$strucwat == 1} {
    readpsf strucwat.psf
    coordpdb strucwat.pdb
  }
  foreach segn $segN segc $segC {
    patch DISU $segn:[lindex $disuList 0] $segc:[lindex $disuList 3]
    patch DISU $segn:[lindex $disuList 1] $segc:[lindex $disuList 4]
    patch DISU $segn:[lindex $disuList 2] $segc:[lindex $disuList 5]
  }
  writepsf $outname.psf
  writepdb $outname.pdb
}
