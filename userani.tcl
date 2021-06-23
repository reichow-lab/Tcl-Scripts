# UserAni 0.92 - VMD Easy User Defined Movie Animation System (Norman Geist 2020; norman.geist@uni-greifswald.de)

# args: material startvalue endvalue
proc userani_FadeTransparency { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 3} {error "$me:improper number of arguments:should be material startval endval!"}
  set val [expr [lindex $funcargs 1] + ($prog * ([lindex $funcargs 2] - [lindex $funcargs 1]))]
  set mat [lindex $funcargs 0]
  eval "material change opacity $mat $val"
}

# args: material property startvalue endvalue
proc userani_FadeMaterialProperty { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 4} {error "$me:improper number of arguments:should be material property startval endval!"}
  set val [expr [lindex $funcargs 2] + ($prog * ([lindex $funcargs 3] - [lindex $funcargs 2]))]
  set mat [lindex $funcargs 0]
  set prop [lindex $funcargs 1]
  eval "material change $prop $mat $val"
}

#args: from to
proc userani_FadeZoom { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 2} {error "$me:improper number of arguments:should be from to!"}
  set val [expr [lindex $funcargs 0] + ($prog * ([lindex $funcargs 1] - [lindex $funcargs 0]))]
  eval "scale to $val"
# }

#args: axis degrees
proc userani_rotate { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 2} {error "$me:improper number of arguments:should be \[x|y|z\] degrees!"}
  set val [expr $part * [lindex $funcargs 1] ]
  eval "rotate [lindex $funcargs 0] by $val"
}

#args: molid axis degrees
#rotate whole molecule, rotating only parts is
#too dangerous and hard to restore
proc userani_rotatemol { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 3} {error "$me:improper number of arguments:should be molid \[x|y|z\] degrees!"}

  set molid [lindex $funcargs 0]
  set axis [lindex $funcargs 1]
  set degrees [lindex $funcargs 2]
  set sel [atomselect $molid "all"]

  if {[$sel num] > 0} {
    set val [expr $play == 1 ? $prog * $degrees : $part * $degrees ]
    set center [measure center $sel]
    set tran [trans center $center axis $axis $val ]
    $sel move $tran
    $sel delete
  } else {error "$me: ups, selection was empty, does molid exist!?"}
}

#args: molid {x y z}
proc userani_movemol { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 2} {error "$me:improper number of arguments:should be molid {x y z}"}

  set molid [lindex $funcargs 0]
  set shift [lindex $funcargs 1]
  set sel [atomselect $molid "all"]

  if {[$sel num] > 0} {
    if {$play == 1} {
      set shift [vecscale $shift $prog]
    } else {
      set shift [vecscale $shift $part]]
    }
    $sel moveby $shift
    $sel delete
  } else {error "$me: ups, selection was empty, does molid exist!?"}
}

#args: factor
proc userani_FadeZoomFromCurrent { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 1} {error "$me:improper number of arguments:should be factor!"}

  #do stuff
  set val [expr pow(10.,log10([lindex $funcargs 0]) * $part)]
  eval "scale by $val"
}

#fade the global translation to 0
proc userani_FadeOutTranslate { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  set global_matrix [lindex [molinfo top get global_matrix] 0]
  set ctrans [list [lindex [lindex $global_matrix 0] 3] [lindex [lindex $global_matrix 1] 3] [lindex [lindex $global_matrix 2] 3]]
  set ntrans [vecscale $ctrans [expr pow(abs($prog-1.),$prog)]]
  eval "translate to $ntrans"
}

#auto scale to bounding box of all thats visible
#regarding DrawMolecule.C a fitting scale value can be
#computed as 1.5/max(range(x),range(y)), Z is neglectable !?
#args [scalefactor]
proc userani_AutoScaleAllVisible { args } {
  global userani_auto_maxprog
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  set minX "false"
  set minY "false"
  set maxX "false"
  set maxY "false"
  set zoom 1.8
  if {[llength $funcargs] > 0 && [string is double [lindex $funcargs 0]]} {
    set zoom [expr $zoom * [lindex $funcargs 0]]
  }
  #compute bb for all visible stuff
  set mols [molinfo list]
  foreach molid $mols {
    if {[molinfo $molid get displayed] && [molinfo $molid get numframes] > 0} {
      set num_reps [molinfo $molid get numreps]
      if {$num_reps > 0} {
	set seltext ""
	for {set i 0} {$i<$num_reps} {incr i} {
	  if {[mol showrep $molid $i]} {
	    if {[string length $seltext] > 0} { set seltext "$seltext or " }
	    set temp [molinfo $molid get "{selection $i}"]
	    set seltext "${seltext}(${temp})"
	  }
	}
      }
      if {[string length $seltext] > 0} {
	set sel [atomselect $molid ($seltext)]
	set mm [measure minmax $sel]
	$sel delete
	set minXtemp [lindex [lindex $mm 0] 0]
	set minYtemp [lindex [lindex $mm 0] 1]
	set maxXtemp [lindex [lindex $mm 1] 0]
	set maxYtemp [lindex [lindex $mm 1] 1]
	set minX [expr $minXtemp < $minX || $minX == "false" ? $minXtemp : $minX]
	set minY [expr $minYtemp < $minY || $minY == "false" ? $minYtemp : $minY]
	set maxX [expr $maxXtemp > $maxX || $maxX == "false" ? $maxXtemp : $maxX]
	set maxY [expr $maxYtemp > $maxY || $maxY == "false" ? $maxYtemp : $maxY]
      }
    }
  }
  if {$minX != "false"} {#true for 1 true for all
    set rangeX [expr $maxX - $minX]
    set rangeY [expr $maxY - $minY]
    set maxrange [expr max($rangeX,$rangeY)]
    set target [expr $zoom/$maxrange]
    set cscale [lindex [lindex [lindex [molinfo top get scale_matrix] 0] 0] 0]
    set nscale [expr $target + (($cscale - $target) * (1.-$userani_auto_maxprog))]
    eval "scale to $nscale"
  } else {
      puts "$me: nothing seems to be visible!"
  }
}

#auto scale to bounding box of visibles in passed molids
#args {molid molid ...} [scalefactor]
proc userani_AutoScaleMolidsVisible { args } {
  global userani_auto_maxprog
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] < 1} {error "$me:improper number of arguments:should be {molid [molid...]}!"}
  set minX "false"
  set minY "false"
  set maxX "false"
  set maxY "false"
  set zoom 1.8
  if {[llength $funcargs] > 1 && [string is double [lindex $funcargs 1]]} {
    set zoom [expr $zoom * [lindex $funcargs 1]]
  }
  #compute bb for all visible stuff
  set mols [lindex $funcargs 0]
  foreach molid $mols  {
    if {[molinfo $molid get displayed] && [molinfo $molid get numframes] > 0} {
      set num_reps [molinfo $molid get numreps]
      if {$num_reps > 0} {
	set seltext ""
	for {set i 0} {$i<$num_reps} {incr i} {
	  if {[mol showrep $molid $i]} {
	    if {[string length $seltext] > 0} { set seltext "$seltext or " }
	    set temp [molinfo $molid get "{selection $i}"]
	    set seltext "${seltext}(${temp})"
	  }
	}
      }
      if {[string length $seltext] > 0} {
	set sel [atomselect $molid ($seltext)]
	set mm [measure minmax $sel]
	$sel delete
	set minXtemp [lindex [lindex $mm 0] 0]
	set minYtemp [lindex [lindex $mm 0] 1]
	set maxXtemp [lindex [lindex $mm 1] 0]
	set maxYtemp [lindex [lindex $mm 1] 1]
	set minX [expr $minXtemp < $minX || $minX == "false" ? $minXtemp : $minX]
	set minY [expr $minYtemp < $minY || $minY == "false" ? $minYtemp : $minY]
	set maxX [expr $maxXtemp > $maxX || $maxX == "false" ? $maxXtemp : $maxX]
	set maxY [expr $maxYtemp > $maxY || $maxY == "false" ? $maxYtemp : $maxY]
      }
    }
  }
  if {$minX != "false"} {#true for 1 true for all
    set rangeX [expr $maxX - $minX]
    set rangeY [expr $maxY - $minY]
    set maxrange [expr max($rangeX,$rangeY)]
    set target [expr $zoom/$maxrange]
    set cscale [lindex [lindex [lindex [molinfo top get scale_matrix] 0] 0] 0]
    set nscale [expr $target + (($cscale - $target) * (1.-$userani_auto_maxprog))]
    eval "scale to $nscale"
  } else {
      puts "$me: nothing seems to be visible!"
  }
}

#auto scale to bounding box of repid in passed molid
#args molid repid [scalefactor]
proc userani_AutoScaleMolRepVisible { args } {
  global userani_auto_maxprog
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] < 2} {error "$me:improper number of arguments:should be molid repid!"}
  set zoom 1.8
  if {[llength $funcargs] > 2 && [string is double [lindex $funcargs 2]]} {
    set zoom [expr $zoom * [lindex $funcargs 2]]
  }
  #compute bb for all visible stuff
  set molid [lindex $funcargs 0]
  set repid [lindex $funcargs 1]
  if {[molinfo $molid get displayed] && [molinfo $molid get numframes] > 0} {
    if {[mol showrep $molid $repid]} {
      set seltext [molinfo $molid get "{selection $repid}"]
      set sel [atomselect $molid ($seltext)]
      set mm [measure minmax $sel]
      $sel delete
      set minX [lindex [lindex $mm 0] 0]
      set minY [lindex [lindex $mm 0] 1]
      set maxX [lindex [lindex $mm 1] 0]
      set maxY [lindex [lindex $mm 1] 1]
      set rangeX [expr $maxX - $minX]
      set rangeY [expr $maxY - $minY]
      set maxrange [expr max($rangeX,$rangeY)]
      set target [expr $zoom/$maxrange]
      set cscale [lindex [lindex [lindex [molinfo top get scale_matrix] 0] 0] 0]
      set nscale [expr $target + (($cscale - $target) * (1.-$userani_auto_maxprog))]
      eval "scale to $nscale"
      return
    }
  }
  puts "$me: nothing seems to be visible!"
}

#auto fit to com of everthing visible
proc userani_AutoFocusAllVisible { args } {
  global userani_auto_maxprog
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  set avg_center {0 0 0}
  set avg_count 0
  #compute center for all visible stuff
  set mols [molinfo list]
  foreach molid $mols {
    if {[molinfo $molid get displayed] && [molinfo $molid get numframes] > 0} {
      set num_reps [molinfo $molid get numreps]
      if {$num_reps > 0} {
	set seltext ""
	for {set i 0} {$i<$num_reps} {incr i} {
	  if {[mol showrep $molid $i]} {
	    if {[string length $seltext] > 0} { set seltext "$seltext or " }
	    set temp [molinfo $molid get "{selection $i}"]
	    set seltext "${seltext}(${temp})"
	  }
	}
      }
      if {[string length $seltext] > 0} {
	set sel [atomselect $molid ($seltext)]
	set center [measure center $sel weight mass]
	$sel delete
	set avg_center [vecadd $avg_center $center]
	incr avg_count
      }
    }
  }
  if {$avg_count > 0} {
    set avg_center [vecscale $avg_center [expr 1./$avg_count]]
    foreach molid $mols {
      eval "molinfo $molid set center [list [list $avg_center]]"
    }
    set global_matrix [lindex [molinfo top get global_matrix] 0]
    set ctrans [list [lindex [lindex $global_matrix 0] 3] [lindex [lindex $global_matrix 1] 3] [lindex [lindex $global_matrix 2] 3]]
    set ntrans [vecscale $ctrans [expr (1.-$userani_auto_maxprog)]]
    eval "translate to $ntrans"
  } else {
      puts "$me: nothing seems to be visible!"
  }
}

#auto fit to com of visibles in passed molids
proc userani_AutoFocusMolidsVisible { args } {
  global userani_auto_maxprog
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 1} {error "$me:improper number of arguments:should be {molid [molid...]}!"}
  set avg_center {0 0 0}
  set avg_count 0
  #compute center for all visible stuff
  set mols [molinfo list]
  foreach molid [lindex $funcargs 0]  {
    if {[molinfo $molid get displayed] && [molinfo $molid get numframes] > 0} {
      set num_reps [molinfo $molid get numreps]
      if {$num_reps > 0} {
	set seltext ""
	for {set i 0} {$i<$num_reps} {incr i} {
	  if {[mol showrep $molid $i]} {
	    if {[string length $seltext] > 0} { set seltext "$seltext or " }
	    set temp [molinfo $molid get "{selection $i}"]
	    set seltext "${seltext}(${temp})"
	  }
	}
      }
      if {[string length $seltext] > 0} {
	set sel [atomselect $molid ($seltext)]
	set center [measure center $sel weight mass]
	$sel delete
	set avg_center [vecadd $avg_center $center]
	incr avg_count
      }
    }
  }
  if {$avg_count > 0} {
    set avg_center [vecscale $avg_center [expr 1./$avg_count]]
    foreach molid $mols {
      eval "molinfo $molid set center [list [list $avg_center]]"
    }
    set global_matrix [lindex [molinfo top get global_matrix] 0]
    set ctrans [list [lindex [lindex $global_matrix 0] 3] [lindex [lindex $global_matrix 1] 3] [lindex [lindex $global_matrix 2] 3]]
    set ntrans [vecscale $ctrans [expr (1.-$userani_auto_maxprog)]]
    eval "translate to $ntrans"
  } else {
      puts "$me: nothing seems to be visible"
  }
}

#auto fit to com of visibles of molid and repid
#args: molid repid
proc userani_AutoFocusMolRepVisible { args } {
  global userani_auto_maxprog
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 2} {error "$me:improper number of arguments:should be molid repid!"}
  #compute center for all visible stuff
  set molid [lindex $funcargs 0] 
  set repid [lindex $funcargs 1]

  if {[molinfo $molid get displayed] && [molinfo $molid get numframes] > 0} {
    if {[mol showrep $molid $repid]} {
      set seltext [molinfo $molid get "{selection $repid}"]
      set sel [atomselect $molid ($seltext)]
      set center [measure center $sel weight mass]
      $sel delete
      foreach molid [molinfo list] {
	eval "molinfo $molid set center [list [list $center]]"
      }
      set global_matrix [lindex [molinfo top get global_matrix] 0]
      set ctrans [list [lindex [lindex $global_matrix 0] 3] [lindex [lindex $global_matrix 1] 3] [lindex [lindex $global_matrix 2] 3]]
      set ntrans [vecscale $ctrans [expr (1.-$userani_auto_maxprog)]]
      eval "translate to $ntrans"
      return
    }
  }
  puts "$me: nothing seems to be visible"
}

# animated rastering
#args: shiftX shiftY
proc userani_raster_replicas { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 2} {error "$me:improper number of arguments:should be shiftx shifty!"}
  global num_replicas

  set shiftX [expr $play == 1 ? [lindex $funcargs 0]*$prog : [lindex $funcargs 0]*$part]
  set shiftY [expr $play == 1 ? [lindex $funcargs 1]*$prog : [lindex $funcargs 1]*$part]

  set replicas_per_row [expr ceil(sqrt($num_replicas))]
  set num_rows [expr ceil($num_replicas/$replicas_per_row)]
  set rowc 0
  set colc 0

  for {set replica_id 0} {$replica_id < $num_replicas} {incr replica_id} {
    set everyone [atomselect $replica_id "all"]
    #moving with center shift <- looks like expanding in all directions
    $everyone moveby [list [expr ($rowc * $shiftX) - (($shiftX/2.) * ($replicas_per_row-1))] [expr ($colc * $shiftY) - (($shiftY/2.) * ($num_rows-1))] 0]
    $everyone delete 

    incr rowc
    if {$rowc == $replicas_per_row} {
      set rowc 0
      incr colc
    }
  }
}


# This is a skeleton functions for userani_ animations
# You will automatically find your arguments in funcargs.
# These functions are especially needed for procedures
# that require any progressive threatment inside frameranges
# like fading etc.
# Therefore use the automatically generated $prog scaling as 0 <  $prog <= 1 progressively
# or the value of $part which gives just a constant factor based on 1/frames2do
# and so can easily scale any parameter you want to fade. No need
# to deal with frames etc from here.
# To invert the order of $prog simply do [expr abs($prog-1)]
# See FadeZoomFromCurrent for complex example
proc userani_skeleton { args } {
  set id [lindex $args 0]
  set prog [lindex $args 1]
  set part [lindex $args 2]
  set play [lindex $args 3]
  set funcargs [lindex $args 4]
  set me [lindex [info level [info level]] 0]
  #----
  if {[llength $funcargs] != 2} {error "$me:improper number of arguments:should be whatever!"}
}

##################################################################################################

### NON userani_* PROCEDURES #####################################################################

#required to keep the rastering static after animation
proc raster_replicas {shiftX shiftY} {
  global num_replicas

  set replicas_per_row [expr ceil(sqrt($num_replicas))]
  set num_rows [expr ceil($num_replicas/$replicas_per_row)]
  set rowc 0
  set colc 0

  set gcenter {0. 0. 0.}
  for {set replica_id 0} {$replica_id < $num_replicas} {incr replica_id} {
    set everyone [atomselect $replica_id "all"]
    #moving with center shift <- looks like expanding in all directions
    $everyone moveby [list [expr ($rowc * $shiftX) - (($shiftX/2.) * ($replicas_per_row-1))] [expr ($colc * $shiftY) - (($shiftY/2.) * ($num_rows-1))] 0]
    $everyone delete 

    incr rowc
    if {$rowc == $replicas_per_row} {
      set rowc 0
      incr colc
    }
  }
}

proc dummy {} {
  sleep 0
  return
}

#reads through all mols representations
#and disables it if the material it uses has
#opacity 0 or enables it if opacity becomes
#bigger 0
#Why? Saves a lot of rendering time to
#not render invisible stuff
proc TransparencyStalker {} {
  set mols [molinfo list]
  foreach molid $mols {
    if {[molinfo $molid get displayed] && [molinfo $molid get numframes] > 0} {
      set num_reps [molinfo $molid get numreps]
      if {$num_reps > 0} {
	for {set i 0} {$i<$num_reps} {incr i} {
	  set material [molinfo $molid get "{material $i}"]
	  set opacity [lindex [material settings $material] 5]
	  if {$opacity > 0} {
	    mol showrep $molid $i 1
	  } else {
	    mol showrep $molid $i 0
	  }
	}
      }
    }
  }
  return
}

##################################################################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
######  DO NOT CHANGE BELOW  ######
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#################################### ANIMATION SYSTEM ############################################

#helpers ------------------------------------------------------------
#is float a whole number?
proc float_whole { float } {
  return [expr abs($float - int($float)) > 0 ? 0 : 1]
}

#get biggest possibe multiple of divisor below divident
proc max_multiple_before {divident divisor} {
  return [expr int($divident/$divisor)*$divisor]
}

#get smallest multiple of divisor above divident
proc first_multiple_from {divident divisor} {
  return [expr ceil(1. * $divident/$divisor) * $divisor];# + (ceil(1. * $divident/$divisor) - int(1. * $divident/$divisor))]
}

#number of multiples in range for step in frame
#this also returns if frame fits the step if the return is a whole
proc frames_in_range {framenum step} {
  return [expr $framenum / (1. * $step)]
}

#helper to calc with percent values inside animation script
proc expr% { args } {
  set result [expr [string map {"%" ""} $args]]
  return "${result}%"
}

proc smooth {{len 3} {molid -1}} {
    set molid [expr $molid < 0 ? [molinfo top] : $molid]
    set nr [molinfo $molid get numreps]
    for {set i 0} {$i < $nr} {incr i} {
        mol smoothrep $molid $i $len
    }
}

#configure max progress for autoscale/focus functions
set userani_auto_maxprog 0.05
proc userani_maxprog { {maxprog 0.05} } {
	global userani_auto_maxprog
	set userani_auto_maxprog $maxprog
}

#prepend mystring with spaces to reach maxlength (align right)
proc spacefill {mystring maxlength} {
  set pre [string repeat " " [expr $maxlength-[string length $mystring]]]
  return "${pre}${mystring}"
}

#calc number of signs correspoding to frame
proc spaceshift {numframes framesperspace} {
    return [expr round(1. * $numframes / $framesperspace)]
}

#env setup ----------------------------------------------------------
array unset useranireg
array set useranireg {}
array unset useranimark
array set useranimark {}
set useranicount 0
#defaults
menu vmdmovie on; #makes sure $::MovieMaker:numframes is set
set ::MovieMaker::movietype userdefined
set ::MovieMaker::renderer dryrun
display update
global useruseranicount useranireg useranimark
#--------------------------------------------------------------------

#args: name [-start start] [-end end] [-step step] \
#           [-shiftstart shift] [-shiftend shift] [-shiftall shift]
#	    [-like sectionname] || [-after sectionname] [-noplay]
#defaulting keywords: first,last, n% where n can be any decimal number
#CAUTION: order of options is relevant.
# -like = copies start,step,end
# -after = copies end as own start
# -shiftall = adds shift to start and end
# -noplay = Trajectory will not play (if section has an action; use dummy for fake action)
proc add_section { args } {
  global useranimark useranidefstart useranidefstep
  if {[llength $args] < 1} {error "Error) Improper number of arguments, need at least a name!"}
  
  #get args
  set name [lindex $args 0]
  set options [lrange $args 1 end]
  set userstep [expr $::MovieMaker::trjstep ]
  set useranidefend [expr $::MovieMaker::numframes * $userstep]
  set useranidefstart 0
  set useranidefstep $userstep

  #setup vars
  set start -1
  set step -1
  set end -1
  set shiftstart 0
  set shiftend 0
  set aftername -1
  set likename -1
  set play 1
  
  #replace defaulting keywords
  set options [string map [list "first" $useranidefstart] $options]
  set options [string map [list "last" $useranidefend] $options]

  set matches [lsort -unique [regexp -all -inline {([0-9]*\.?[0-9]*\%)} $options]]
  foreach match [lsort -unique $matches] {
    set perc [string map {"%" ""} $match]
    set perc [expr $perc/100. * $useranidefend]
    set perc [max_multiple_before $perc $userstep]
    set options [string map [list $match $perc] $options]
  }

  #check args
  for {set i 0} {$i < [llength $options]} {incr i} {
    if {[lindex $options $i] == "-start"} {
      set start [lindex $options [expr $i+1]]
      if {![string is integer $start]} {error "Error) -start needs an integer value!"}
      incr i
    } elseif {[lindex $options $i] == "-end"} {
      set end [lindex $options [expr $i+1]]
      if {![string is integer $end]} {error "Error) -end needs an integer value!"}
      incr i
    } elseif {[lindex $options $i] == "-step"} {
      set step [expr [lindex $options [expr $i+1]] * $userstep]
      if {![string is integer $step]} {error "Error) -step needs an integer value!"}
      incr i
    } elseif {[lindex $options $i] == "-after"} {
      set aftername [lindex $options [expr $i+1]]
      if {[string length $aftername] < 1} {error "Error) -after needs a value!"}
      incr i
      #copy start value from aftername
      if {[info exists useranimark($aftername)]} {
	set start [lindex $useranimark($aftername) 2]
      } else { error "Error) Section name for after $aftername undefined!" }
    } elseif {[lindex $options $i] == "-shiftstart"} {
      set shiftstart [lindex $options [expr $i+1]]
      if {![string is integer $shiftstart]} {error "Error) -shiftstart needs an integer value!"}
      incr i
    } elseif {[lindex $options $i] == "-shiftend"} {
      set shiftend [lindex $options [expr $i+1]]
      if {![string is integer $shiftend]} {error "Error) -shiftend needs an integer value!"}
      incr i
    } elseif {[lindex $options $i] == "-shiftall"} {
      set shiftstart [lindex $options [expr $i+1]]
      set shiftend [lindex $options [expr $i+1]]
      if {![string is integer $shiftend]} {error "Error) -shiftall needs an integer value!"}
      incr i
    } elseif {[lindex $options $i] == "-noplay"} {
      set play -1
    } elseif {[lindex $options $i] == "-like"} {
      set likename [lindex $options [expr $i+1]]
      if {[string length $likename] < 1} {error "Error) -like needs a value!"}
      incr i
      #copy values from likename
      if {[info exists useranimark($likename)]} {
	set start [lindex $useranimark($likename) 0]
	set step [lindex $useranimark($likename) 1]
	set end [lindex $useranimark($likename) 2]
      } else { error "Error) Section name for like $likename undefined!" }
    }
  }

  #insert whats not already theres
  if {$start == -1} { set start $useranidefstart }
  if {$step == -1} { set step $useranidefstep }
  if {$end == -1} { set end [expr $start+$step] }

  #apply shifts
  set start [expr $start+$shiftstart]
  set end [expr $end+$shiftend]

  #check values for being integer
  foreach val [list $start $step $end] {
    if {![string is integer $val]} {error "Error) Invalid settings in section $name \{$start $step $end\}!"}
  }

  #check range settings
  #check if section can actually happen due step size
  #should also avoid cased like start > end or start + step > end etc.
  set realstart [first_multiple_from $start $step]
  if {[frames_in_range [expr $end-$realstart] $step] < 1} {
    puts "Warning) Section '$name' to short for step size $step, would be overjumped."
    puts "Extending to minimum length (step+1)!"
    set end [expr int($realstart + $step)]
  }
  #negative values
  if {$start < 0 || $step < 0 || $end < 0} {
    error "Error) Section $name negative values for start,step or end!"
  }
  #bounds
  if {$start > $useranidefend || $step > $useranidefend || $end > $useranidefend} {
    error "Error) Section $name values out of bounds \{$start $step $end\} (frame max = $useranidefend)!"
  }
  
  #setup range
  set range [list $start $step $end $play]

  #check if we overwrite
  if {[info exists useranimark($name)]} {
    puts "Warning) Section name $name overwritten!"
  }

  set useranimark($name) $range
}

proc add_userani {sectionname function args} {
  global useranicount useranireg useranimark
  if {![info exists useranimark($sectionname)]} {
    puts "Warning) Section name $sectionname not yet defined!"
  }
  set useranireg($useranicount) [list $sectionname $function $args 0]
  incr useranicount
}

proc add_eval {sectionname function args} {
  global useranicount useranireg useranimark
  if {![info exists useranimark($sectionname)]} {
    puts "Warning) Section name $sectionname not yet defined!"
  }
  set useranireg($useranicount) [list $sectionname $function $args 1]
  incr useranicount
}

proc print_story {} {
  global useranireg 
  puts "\nInfo) Storyboard: (ordered like added)"
  foreach ani [lsort -integer -increasing [array names useranireg]] {
    puts "[lindex $useranireg($ani) 0]\t\t[lrange $useranireg($ani) 1 2]"
  }
  puts ""
}

proc print_sections {} {
  global useranimark
  array set temp {}
  puts "\nInfo) Time sections: (ordered by start time)"
  foreach sect [array names useranimark] {
    set temp([list [lindex $useranimark($sect) 0] [lindex $useranimark($sect) 2] $sect]) $sect
  }
  foreach sect [lsort -increasing -dictionary [array names temp]] {
    puts "$temp($sect) \t\t \[$useranimark($temp($sect))\]"
  }
  array unset temp
  puts ""
}

proc print_timeline { {scale 1} } {
  global useranimark useranireg
  set userstep [expr $::MovieMaker::trjstep ]
  set useranidefend [expr $::MovieMaker::numframes * $userstep]
  set lengthperval [string length $useranidefend]
  set scale [expr int($scale*10.)/10.]; #restrict to one digit behind comma;
  set area [expr int($scale * 70.)]; #fixed width for timeline
  set numticks [expr int(10. * $scale)]; #num ticks
  set ticksspacing [expr 1. * $area / $numticks];
  set framesperspace [expr 1. * $useranidefend / $area]; #num frames each char represents
  set framerate $MovieMaker::framerate

  puts "\nInfo) Timeline: (* = freeze)"

  #reorder sections and find largest section name
  set maxsectname 0
  array set temp {}
  foreach sect [array names useranimark] {
    set temp([list [lindex $useranimark($sect) 0] [lindex $useranimark($sect) 2] $sect]) $sect
    set maxsectname [expr [string length $sect] > $maxsectname ? [string length $sect] : $maxsectname]
  }

  #print axis and labels
  set axis "|"
  set framebar "0"
  set percbar "0%"
  set secbar "0s"
  for {set i 0} {$i<$numticks} {incr i} {
    set tick [string repeat "_" [expr round($ticksspacing-1)]]
    set axis "$axis$tick|"

    set framelabel [expr round(($i+1.)*$ticksspacing*$framesperspace)]
    set tick [string repeat " " [expr round($ticksspacing-[string length $framelabel])]]
    set framebar "$framebar$tick$framelabel"

    set perclabel [expr round(1. * $framelabel/$useranidefend * 100)]
    set perclabel "$perclabel%"
    set tick [string repeat " " [expr round($ticksspacing-[string length $perclabel])]]
    set percbar "$percbar$tick$perclabel"
    
    set seclabel [expr round( 1. * $framelabel / ($framerate * $userstep))]
    set seclabel "${seclabel}s"
    set tick [string repeat " " [expr round($ticksspacing-[string length $seclabel])]]
    set secbar "$secbar$tick$seclabel"
  }
  set shift [string repeat " " [expr $maxsectname + 1 + (($lengthperval * 3) +2) +3]]; #the +n are just the spaces and braces
  puts "$shift$framebar"
  puts "$shift$axis"

  #actually output time sections
  foreach sect [lsort -increasing -dictionary [array names temp]] {
    set start [lindex $useranimark($temp($sect)) 0]
    set step [lindex $useranimark($temp($sect)) 1]
    set end [lindex $useranimark($temp($sect)) 2]
    set play [lindex $useranimark($temp($sect)) 3]
    set playchar "-"
    if { $play < 0 } { set playchar "*" }
    set name [spacefill $temp($sect) $maxsectname]
    set range "\[[spacefill $start $lengthperval] [spacefill $step $lengthperval] [spacefill $end $lengthperval]\]"

    set from [spaceshift $start $framesperspace]
    set to [spaceshift $end $framesperspace]
    
    set align [string repeat " " $from]
    set length [string repeat $playchar [expr $to-$from-1]]

    puts "$name $range $align|$length>"
  }
  #print lower % bar
  puts "$shift$axis"
  puts "$shift$percbar"
  puts "$shift$secbar"

  array unset temp
  puts ""
}

proc followregister { args } {
  global useranicount useranireg useranimark
  set userstep [expr $::MovieMaker::trjstep ]
  set userframe [expr $::MovieMaker::userframe * $userstep]
  set useranidefend [expr $::MovieMaker::numframes * $userstep]
  set restoreframe -1; #no play

  #goto 0 on start
  if {$userframe == 0} { animate goto 0 }

  #look for noplay sections around current frame
  foreach ani [lsort -integer -increasing [array names useranireg]] {
    set framerange $useranimark([lindex $useranireg($ani) 0])
    set start [lindex $framerange 0]
    set step [lindex $framerange 1]
    set end [lindex $framerange 2]
    set play [lindex $framerange 3]
    set stepcond [frames_in_range $userframe $step]
    if {$play < 0 && $userframe >= $start && $userframe < $end} {#noplay?
      set restoreframe [expr $start < $restoreframe || $restoreframe == -1 ? $start : $restoreframe]
    } elseif {$play > 0 && $userframe >= $start && $userframe < $end && [float_whole $stepcond]} {#slow play
      set restoreframe -1
    }
  }
  if {$restoreframe == -1} {#we may play
    animate goto $userframe
  }

  foreach ani [lsort -integer -increasing [array names useranireg]] {
    set framerange $useranimark([lindex $useranireg($ani) 0])
    set start [lindex $framerange 0]
    set step [lindex $framerange 1]
    set end [lindex $framerange 2]
    set play [lindex $framerange 3]
    set funcargs [lindex $useranireg($ani) 2]
    set is_eval [lindex $useranireg($ani) 3]
    if {$is_eval == 1} {#call any tcl stuff
      set function [lindex $useranireg($ani) 1]
      set stepcond [frames_in_range $userframe $step]
      if {$userframe >= $start && $userframe < $end && [float_whole $stepcond]} {
	eval "$function $funcargs"
      }
    } else {#call userani_* procedures
      set function "userani_[lindex $useranireg($ani) 1]"
      set stepcond [frames_in_range $userframe $step]
      if {$userframe >= $start && $userframe < $end && [float_whole $stepcond]} {
	set done [frames_in_range [expr $userframe-$start] $step]
	set todo [frames_in_range [expr $end-$start] $step]
	set prog [expr ($done + 1) / $todo]
	set part [expr 1. / $todo]
	eval "$function $ani $prog $part $play {$funcargs}"
      }
    }
  }
}

proc updatesettings { args } {
  set ::MovieMaker::numframes [max_multiple_before [expr int([molinfo top get numframes] / $::MovieMaker::trjstep)] $::MovieMaker::framerate]
  set ::MovieMaker::movieduration [expr int($::MovieMaker::numframes / $::MovieMaker::framerate)]
  if {$::MovieMaker::trjstep > $::MovieMaker::numframes} {
    error "Error) Trajectory step size exceeds number of frames!"
  } elseif {$::MovieMaker::numframes > [molinfo top get numframes]} {
    error "Error) Not enough frames, check your top mol!"
  }
}

proc onchangewarn { args } {
  updatesettings
  puts "Warning) Some settings changed, please re-source your animation procedure to update the time sections!"
}

proc start_userani { }  {
  trace add variable ::MovieMaker::userframe write followregister
  trace add variable ::MovieMaker::trjstep write onchangewarn
}

proc stop_userani { } {
  trace remove variable ::MovieMaker::userframe write followregister
  trace remove variable ::MovieMaker::trjstep write onchangewarn
}
#---------------------------------------------------------------------

# Safe to re-source --
stop_userani
start_userani
updatesettings; #update numframes f.i. due top mol change
#-----------------
