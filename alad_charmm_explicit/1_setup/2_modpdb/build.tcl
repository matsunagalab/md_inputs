# Load the original PDB
mol load pdb ../1_oripdb/ala2.pdb

# Rename "PDB general atom name" to "CHARMM-specific atom name"
#   HIS => HSD (but not included in this protein)
#   CD1 atom of ILE => CD
#   C-terminal carboxyl oxygen O and OXT => OT1 and OT2
#[atomselect top "resname HIS"                      ] set resname HSD
#[atomselect top "resname ILE and name CD1"         ] set name CD
#[atomselect top "resid 378 and name O"  ] set name OT1
#[atomselect top "resid 378 and name OXT"] set name OT2

# Measure the center of mass (com) of the selected atoms (protein)
# Shift the com of the protein to the origin
set sel [atomselect top "all"]
set com [measure center $sel weight mass]
$sel moveby [vecscale -1.0 $com]

# Write the modified PDB of the selected atoms
set sel [atomselect top "all"]
$sel writepdb prot.pdb

exit
