# Load psfgen-plugin and CHARMM topology file
package require psfgen
resetpsf
topology ../toppar/top_all36_prot.rtf

# Define the segment name as "PROA"
segment PROA {pdb ../2_modpdb/proa.pdb
              first ACE
              last CT3} 
#segment PROB {pdb ../2_modpdb/prob.pdb
#              first ACE
#              last CT3} 
segment PROC {pdb ../2_modpdb/proc.pdb
              first ACE
              last CT3} 
#segment PROD {pdb ../2_modpdb/prod.pdb
#              first ACE
#              last CT3} 
segment PROE {pdb ../2_modpdb/proe.pdb
              first ACE
              last CT3} 
#segment PROF {pdb ../2_modpdb/prof.pdb
#              first ACE
#              last CT3} 
segment PROG {pdb ../2_modpdb/prog.pdb
              first ACE
              last CT3} 
#segment PROH {pdb ../2_modpdb/proh.pdb
#              first ACE
#              last CT3} 
segment PROI {pdb ../2_modpdb/proi.pdb
              first ACE
              last CT3} 
#segment PROJ {pdb ../2_modpdb/proj.pdb
#              first ACE
#              last CT3} 

# Assign the coordinates of atoms to PROA using the PDB data
coordpdb ../2_modpdb/proa.pdb PROA
#coordpdb ../2_modpdb/prob.pdb PROB
coordpdb ../2_modpdb/proc.pdb PROC
#coordpdb ../2_modpdb/prod.pdb PROD
coordpdb ../2_modpdb/proe.pdb PROE
#coordpdb ../2_modpdb/prof.pdb PROF
coordpdb ../2_modpdb/prog.pdb PROG
#coordpdb ../2_modpdb/proh.pdb PROH
coordpdb ../2_modpdb/proi.pdb PROI
#coordpdb ../2_modpdb/proj.pdb PROJ

# Guess the coordinates of missing atoms (mainly hydrogen)
guesscoord

regenerate angles dihedrals

# Generate PDB and PSF files
writepdb protein.pdb
writepsf protein.psf
exit

