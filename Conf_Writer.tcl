##########################################################################
##	Preparation of MD Configuration Files				##
##									##
##	Portland State University					##
##	P.I.	: Steve Reichow						##
##	Author	: Matthew Veter						##
##########################################################################

set outname "OUTNAME"
set all [atomselect top all]
set box [atomselect top "all and not lipids" ]
set center [measure center $box]
set minmax [measure minmax $box]

proc bounds {arg1 arg2} {
expr ([lindex $arg1 1 $arg2] - [lindex $arg1 0 $arg2])
}

set suffix [list \
_LipMelt \
_ProtConst \
_BackConst \
_Equil \
_Produc \
]

set configs [list \
Lipid-Melt.conf \
Protein-Constrained.conf \
Backbone-Constrained.conf \
Equilibration.conf \
Production.conf \
]

foreach ending $suffix cfg $configs {

	set ofile [open $cfg w]

	puts $ofile [format "structure\t\t$outname\_popcwi.psf"]
	puts $ofile [format "coordinates\t\t$outname\_popcwi.pdb"]
	puts $ofile [format "outputName\t\t$outname$ending\n"]
	if {$ending == "_LipMelt"} {
		puts $ofile [format "temperature\t\t300\n"]
		puts $ofile "if \{0\} \{"
	} else {
		puts $ofile [format "temperature\t\t310\n"]
		puts $ofile "if \{1\} \{"
	}
	puts $ofile [format "set inputname\t\t$outname$ending"]
	puts $ofile [format "binCoordinates\t\t\$inputname.restart.coor"]
	puts $ofile [format "binVelocities\t\t\$inputname.restart.vel"]
	puts $ofile [format "extendedSystem\t\t\$inputname.restart.xsc\n\}\n"]
	puts $ofile [format "firsttimestep\t\t0\n"]
	puts $ofile [format "paraTypeCharmm\t\ton"]
	puts $ofile [format "parameters\t\t/home/bassam/Topology\ and\ Parameters/par_all36m_prot.prm"]
	puts $ofile [format "mergeCrossterms\t\tyes"]
	puts $ofile [format "parameters\t\t/home/bassamh/Topology\ and\ Parameters/par_all36_lipid.prm"]
	puts $ofile [format "parameters\t\t/home/bassamh/Topology\ and\ Parameters/toppar_water_ions_namd.str"]
	puts $ofile [format "parameters\t\t/home/bassamh/Topology\ and\ Parameters/par_all36_cgenff.prm"]
	puts $ofile [format "parameters\t\t/home/bassamh/Topology\ and\ Parameters/par_all36_carb.prm"]
	puts $ofile [format "parameters\t\t/home/bassamh/Topology\ and\ Parameters/cAMP.prm\n"]
	puts $ofile [format "temperature\t\t \$temperature\n"]
	if {$ending == "_LipMelt"} {
		puts $ofile "if \{1\} \{"
	} else {
		puts $ofile "if \{0\} \{"
	}
	puts $ofile [format "cellBasisVector1\t[bounds $minmax 0]\t0.\t0."]
	puts $ofile [format "cellBasisVector2\t0.\t[bounds $minmax 1]\t0."]
	puts $ofile [format "cellBasisVector3\t0.\t0.\t[bounds $minmax 2]"]
	puts $ofile [format "cellOrigin\t\t[lindex $center 0]\t[lindex $center 1]\t[lindex $center 2]\n\}\n"]
	puts $ofile [format "wrapWater\t\ton"]
	puts $ofile [format "wrapAll\t\t\ton\n"]
	puts $ofile "# Force-Field Parameters"
	puts $ofile [format "exclude\t\t\tscaled1-4"]
	puts $ofile [format "1-4scaling\t\t1.0"]
	puts $ofile [format "cutoff\t\t\t12.0"]
	puts $ofile [format "switching\t\ton"]
	puts $ofile [format "switchdist\t\t10.0"]
	puts $ofile [format "pairlistdist\t\t13.5\n"]
	puts $ofile "# Integrator Parameter"
	if {$ending == "_LipMelt" || $ending == "_ProtConst" || $ending == "_BackConst"} {
		puts $ofile [format "timestep\t\t1.0"]
	} else {puts $ofile [format "timestep\t\t2.0"]
	}
	puts $ofile [format "rigidBonds\t\tall"]
	puts $ofile [format "nonbondedFreq\t\t1"]
	puts $ofile [format "fullElectFrequency\t4\t;# product of this and timestep must not exceed 8 (leave this @ 4)"]
	puts $ofile [format "stepspercycle\t\t20\n"]
	puts $ofile [format "PME\t\t\ton"]
	puts $ofile [format "PMEGridSpacing\t\t1.0\n"]
	puts $ofile "# Constant Temperature"
	puts $ofile [format "langevin\t\ton"]
	puts $ofile [format "langevinDamping\t\t1"]
	puts $ofile [format "langevinTemp\t\t\$temperature\n"]
	puts $ofile "# Constant Pressure"
	if {$ending == "_LipMelt"} {
		puts $ofile "if \{0\} \{"
	} else {
		puts $ofile "if \{1\} \{"
	}
	puts $ofile [format "useGroupPressure\tyes ;# needed for 2fs steps"]
	puts $ofile [format "useFlexibleCell\t\tyes ;# no for water box, yes for membrane"]
	puts $ofile [format "useConstantArea\t\tyes ;# no for water box, yes for membrane\n"]
	puts $ofile [format "langevinPiston\t\ton"]
	puts $ofile [format "langevinPistonTarget\t1.01325 ;# in bar --> 1 atm"]
	puts $ofile [format "langevinPistonPeriod\t200.0"]
	puts $ofile [format "langevinPistonDecay\t50.0"]
	puts $ofile [format "langevinPistonTemp\t\$temperature\n\}\n"]
	puts $ofile [format "restartfreq\t\t5000"]
	puts $ofile [format "dcdfreq\t\t\t1000"]
	puts $ofile [format "xstfreq\t\t\t1000"]
	puts $ofile [format "outputEnergies\t\t500"]
	puts $ofile [format "outputPressure\t\t500\n"]
	puts $ofile "# Steered MD"
	puts $ofile [format "SMD\t\t\toff"]
	puts $ofile [format "SMDFile\t\t\tlip-ref.pdb"]
	puts $ofile [format "SMDVel\t\t\t0"]
	puts $ofile [format "SMDk\t\t\t5"]
	puts $ofile [format "SMDDir\t\t\t0 0 1"]
	puts $ofile [format "SMDOutputFreq\t\t10000\n"]
	puts $ofile "# Steered restraint for cAMP"
	puts $ofile [format "colvars\t\t\toff"]
	puts $ofile [format "colvarsConfig\t\t.in\n"]
	puts $ofile "# Constrained Atoms"
	if {$ending == "_ProtConst" || $ending == "_BackConst"} {
		puts $ofile "if \{1\} \{"
	} else {
		puts $ofile "if \{0\} \{"
	}
	puts $ofile [format "constraints\t\toff"]
	puts $ofile [format "consexp\t\t\t2"]
	puts $ofile [format "consref\t\t\t.pdb"]
	puts $ofile [format "conskfile\t\t.pdb"]
	puts $ofile [format "conskcol\t\tB\n\}\n"]
	puts $ofile "# Minimization"
	puts $ofile "if \{1\} \{"
	puts $ofile [format "minBabyStep\t\t1.0e-3"]
	puts $ofile [format "minimize\t\t1000"]
	puts $ofile [format "reinitvels\t\t\$temperature\n\}\n"]
	if {$ending == "_LipMelt" || $ending == "_ProtConst" || $ending == "_BackConst"} {
		puts $ofile "run 999000 ;# 1 ns"
	} elseif {$ending == "_Equil"} {
		puts $ofile "run 13499000 ;# 27 ns"
	} else {puts $ofile "run 25000000 ;# 50 ns"
	}
	close $ofile
}
