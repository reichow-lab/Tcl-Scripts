puts "To run namdenergy type the following: Calc_Energy 'outfile name' 'time-step'"

proc Calc_Energy {ofile temp timemult sel1 sel2} {

	package require namdenergy
	namdenergy -elec -vdw -nonb -sel $sel1 $sel2 -T 310 -ofile $ofile -tempname $temp -diel 1 -timemult $timemult -stride 1 -par /media/bassam/Lab_Files/Topology\ and\ Parameters/par_all36_carb.prm -par /media/bassam/Lab_Files/Topology\ and\ Parameters/par_all36_cgenff.prm -par /media/bassam/Lab_Files/Topology\ and\ Parameters/par_all36_lipid.prm -par /media/bassam/Lab_Files/Topology\ and\ Parameters/par_all36m_prot.prm -par /media/bassam/Lab_Files/Topology\ and\ Parameters/toppar_water_ions_namd.str
}
