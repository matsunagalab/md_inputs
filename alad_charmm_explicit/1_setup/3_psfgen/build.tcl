# Load psfgen-plugin and CHARMM topology file
package require psfgen
resetpsf
topology ../toppar/top_all36_prot.rtf

# Define the segment name as "PROA"
segment PROA {pdb ../2_modpdb/prot.pdb
              first NONE
              last NONE} 

# Assign the coordinates of atoms to PROA using the PDB data
coordpdb ../2_modpdb/prot.pdb PROA

# Guess the coordinates of missing atoms (mainly hydrogen)
guesscoord

regenerate angles dihedrals

# Generate PDB and PSF files
writepdb protein.pdb
writepsf protein.psf
exit

