# julia> using MDToolbox
# julia> psf = mdload("../3_psfgen/protein.psf")
# julia> pdb = mdload("../3_psfgen/protein.pdb, top=psf")
# 
# julia> minimum(pdb.xyz[1, 1:3:end])
# -57.367
# 
# julia> maximum(pdb.xyz[1, 1:3:end])
# 57.701
# 
# julia> 
# 
# julia> minimum(pdb.xyz[1, 2:3:end])
# -45.45
# 
# julia> maximum(pdb.xyz[1, 2:3:end])
# 46.818
# 
# julia> 
# 
# julia> minimum(pdb.xyz[1, 3:3:end])
# -19.729
# 
# julia> maximum(pdb.xyz[1, 3:3:end])
# 23.927

# Solvate the protein in 140x140x140 A^3 water box
package require solvate
set psffile "../3_psfgen/protein.psf"
set pdbfile "../3_psfgen/protein.pdb"
solvate $psffile $pdbfile -minmax {{-70.0 -70.0 -70.0} {70.0 70.0 70.0}} -o wbox 
exit
