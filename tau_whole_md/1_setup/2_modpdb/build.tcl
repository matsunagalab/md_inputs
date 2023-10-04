# Load the original PDB
mol load pdb ../1_oripdb/5O3L.pdb

# Rename "PDB general atom name" to "CHARMM-specific atom name"
#   HIS => HSD (but not included in this protein)
#   CD1 atom of ILE => CD
#   C-terminal carboxyl oxygen O and OXT => OT1 and OT2
[atomselect top "resname HIS"                      ] set resname HSD
[atomselect top "resname ILE and name CD1"         ] set name CD
[atomselect top "resid 378 and name O"  ] set name OT1
[atomselect top "resid 378 and name OXT"] set name OT2

# Measure the center of mass (com) of the selected atoms (protein)
# Shift the com of the protein to the origin
set sel [atomselect top "protein"]
set com [measure center $sel weight mass]
$sel moveby [vecscale -1.0 $com]

# Write the modified PDB of the selected atoms
set sel [atomselect top "chain A"]
$sel writepdb proa.pdb

set sel [atomselect top "chain B"]
$sel writepdb prob.pdb

set sel [atomselect top "chain C"]
$sel writepdb proc.pdb

set sel [atomselect top "chain D"]
$sel writepdb prod.pdb

set sel [atomselect top "chain E"]
$sel writepdb proe.pdb

set sel [atomselect top "chain F"]
$sel writepdb prof.pdb

set sel [atomselect top "chain G"]
$sel writepdb prog.pdb

set sel [atomselect top "chain H"]
$sel writepdb proh.pdb

set sel [atomselect top "chain I"]
$sel writepdb proi.pdb

set sel [atomselect top "chain J"]
$sel writepdb proj.pdb

exit
