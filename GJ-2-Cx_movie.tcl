# This program will take a .dcd trajectory, then superimpose the individual chains of a 	
# Gapjunction onto one another, and allow you to watch a movie of the superimposed chains.	
# 												
# First it will intake an atomselect ion which it will align to. With this information		
# given the program will have to go through a set of embedded for-loops...			
#												
#	The first loop will choose the chain, and the second loop will iterate through every	
#	frame of the trajectory applying the the transformation (alignment) to each frame.	
#												
proc gj2cx {seltext} {
	set chainsel [atomselect top "protein and resid 9 and name CA"]
	set chainlist [$chainsel get chain]
	set ref_mol	[atomselect top "protein and chain A and $seltext" frame 0]
	set numF	[molinfo top get numframes]
	foreach chain $chainlist {
		for {set f 0} {$f < $numF} {incr f} {
			animate goto $f
			set align_mol	 	[atomselect top "protein and chain $chain and $seltext"]
			set chainsel		[atomselect top "protein and chain $chain"]
			set trans_matrix	[measure fit $align_mol $ref_mol]
			$chainsel move $trans_matrix
		}
	}
}


