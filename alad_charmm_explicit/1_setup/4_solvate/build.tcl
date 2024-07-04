# Solvate the protein in 140x140x140 A^3 water box
package require solvate
set psffile "../3_psfgen/protein.psf"
set pdbfile "../3_psfgen/protein.pdb"
solvate $psffile $pdbfile -minmax {{-20.0 -20.0 -20.0} {20.0 20.0 20.0}} -o wbox 
exit
