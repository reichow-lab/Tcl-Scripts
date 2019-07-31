## 	For this script it will need to have a psf, and a pdb loaded...
##
##	Creates the .fix and both .cnst files

proc run {outname} {

	set all [atomselect top all]

	set fix [atomselect top "protein or water or name CLA SOD POT or (resname POPC and name N P HA HB HS HX HY C1 C2 C3 C11 H11A H11B O11 C12 H12A H12B O12 C13 H13A H13B H13C O13 C14 H14A H14B H14C O14 C15 H15A H15B H15C C21 O21 C31 O31 O22 O32)" ]

	set prot [atomselect top protein]

	set back [atomselect top "protein and backbone" ]

	$all set beta 0
	$fix set beta 1
	set end1 "_popcwi.fix"
	$all writepdb $outname$end1

	$all set beta 0
	$prot set beta 1
	set end2 "_PROT.cnst"
	$all writepdb $outname$end2

	$all set beta 0
	$back set beta 1
	set end3 "_BACK.cnst"
	$all writepdb $outname$end3


##
## 	Creates the config files

	set box [atomselect top "all and not lipids"]
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

	set newending ""

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
		puts $ofile [format "set inputname\t\t$outname$newending"]

		set newending $ending

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
		puts $ofile [format "temperature\t\t\$temperature\n"]
		if {$ending == "_LipMelt"} {
			puts $ofile "if \{1\} \{"
		} else {
			puts $ofile "if \{0\} \{"
		}
		puts $ofile [format "cellBasisVector1\t[bounds $minmax 0]\t0.\t0."]
		puts $ofile [format "cellBasisVector2\t0.\t[bounds $minmax 1]\t0."]
		puts $ofile [format "cellBasisVector3\t0.\t0.\t[bounds $minmax 2]"]
		puts $ofile [format "cellOrigin\t\t[lindex $center 0]\t[lindex $center 1]\t[lindex $center 2]\n\}\n"]
		puts $ofile [format "wrapWater\t\toff"]
		puts $ofile [format "wrapAll\t\t\toff\n"]

		puts $ofile "# Force-Field Parameters"
		puts $ofile [format "exclude\t\t\tscaled1-4"]
		puts $ofile [format "1-4scaling\t\t1.0"]
		puts $ofile [format "cutoff\t\t\t12.0"]
		puts $ofile [format "switching\t\ton"]
		puts $ofile [format "switchdist\t\t10.0"]
		puts $ofile [format "pairlistdist\t\t13.5\n"]

		puts $ofile "# Integrator Parameter"
		puts $ofile [format "timestep\t\t2.0"]
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
		puts $ofile [format "SMD\t\t\ton"]
		puts $ofile [format "SMDFile\t\t\tlip-ref.pdb"]
		puts $ofile [format "SMDVel\t\t\t0"]
		puts $ofile [format "SMDk\t\t\t5"]
		puts $ofile [format "SMDDir\t\t\t0 0 -1"]
		puts $ofile [format "SMDOutputFreq\t\t10000\n"]

		puts $ofile "# Steered restraint for cAMP"
		puts $ofile [format "colvars\t\t\ton"]
		puts $ofile [format "colvarsConfig\t\tcAMP-SMD.in\n"]

		puts $ofile "# Fixed Atoms Constraint"
		if {$ending == "_LipMelt"} {
			puts $ofile "if \{1\} \{"
		} else {
			puts $ofile "if \{0\} \{"
		}
		puts $ofile [format "fixedAtoms\t\ton"]
		puts $ofile [format "fixedAtomesFile\t\t$outname$end1"]
		puts $ofile [format "fixedAtomsCol\t\tB"]
		puts $ofile [format "fixedAtomsForces\ton\n\}\n"]

		puts $ofile "# Constrained Atoms (prevent drift from SMD...only restraining CA from resid 208)"
		if {$ending == "_ProtConst" || $ending == "_BackConst"} {
			puts $ofile "if \{1\} \{"
		} else {
			puts $ofile "if \{0\} \{"
		}
		puts $ofile [format "constraints\t\toff"]
		puts $ofile [format "consexp\t\t\t2"]
		puts $ofile [format "consref\t\t\tSMD.pdb"]
		puts $ofile [format "conskfile\t\tSMD.pdb"]
		puts $ofile [format "conskcol\t\tB\n\}\n"]

		puts $ofile "# Minimization"
		puts $ofile "if \{1\} \{"
		puts $ofile [format "minBabyStep\t\t1.0e-3"]
		puts $ofile [format "minimize\t\t1000"]
		puts $ofile [format "reinitvels\t\t\$temperature\n\}\n"]
		if {$ending == "_LipMelt"} {
			puts $ofile "run 1000000 ;# 2 ns"
		} else {
			puts $ofile "run 5000000 ;# 10 ns"
		}
		close $ofile
	}
}
