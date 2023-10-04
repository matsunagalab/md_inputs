# Load psfgen-plugin and CHARMM topology file
package require psfgen
resetpsf
topology ../toppar/top_all36_prot.rtf

# Define the segment name as "PROA"
segment PROA {pdb ../2_modpdb/proa.pdb
              first ACE
              last CT3} 
segment PROB {pdb ../2_modpdb/prob.pdb
              first ACE
              last CT3} 
segment PROC {pdb ../2_modpdb/proc.pdb
              first ACE
              last CT3} 
segment PROD {pdb ../2_modpdb/prod.pdb
              first ACE
              last CT3} 
segment PROE {pdb ../2_modpdb/proe.pdb
              first ACE
              last CT3} 

# Assign the coordinates of atoms to PROA using the PDB data
coordpdb ../2_modpdb/proa.pdb PROA
coordpdb ../2_modpdb/prob.pdb PROB
coordpdb ../2_modpdb/proc.pdb PROC
coordpdb ../2_modpdb/prod.pdb PROD
coordpdb ../2_modpdb/proe.pdb PROE

# Guess the coordinates of missing atoms (mainly hydrogen)
guesscoord

regenerate angles dihedrals

# Generate PDB and PSF files
writepdb protein.pdb
writepsf protein.psf
exit

