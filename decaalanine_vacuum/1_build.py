import pdbfixer
import openmm as mm
import openmm.app as app
from openmm import unit
import sys
import numpy as np

pdbfile = 'ala9.pdb'

def prepare_protein(pdb_file, ignore_missing_residues=True, ignore_terminal_missing_residues=True, ph=7.0):
    """
    Use pdbfixer to prepare the protein from a PDB file. Hetero atoms such as ligands are
    removed and non-standard residues replaced. Missing atoms to existing residues are added.
    Missing residues are ignored by default, but can be included.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ignore_missing_residues: bool, optional
        If missing residues should be ignored or built.
    ignore_terminal_missing_residues: bool, optional
        If missing residues at the beginning and the end of a chain should be ignored or built.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.

    CC-BY 4.0, TeachOpenCADD T019,  Molecular dynamics simulation
    (c) Pietro Gerletti, Mareike Leja, Jeffrey R Wagner, David Schaller, Andrea Volkamer, 2020.
    https://projects.volkamerlab.org/teachopencadd/talktorials/T019_md_simulation.html
    """
    fixer = pdbfixer.PDBFixer(str(pdb_file))
    fixer.removeHeterogens()  # co-crystallized ligands are unknown to PDBFixer
    fixer.findMissingResidues()  # identify missing residues, needed for identification of missing atoms

    # if missing terminal residues shall be ignored, remove them from the dictionary
    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

    # if all missing residues shall be ignored ignored, clear the dictionary
    if ignore_missing_residues:
        fixer.missingResidues = {}

    # Add ACE or NME cap  
    for chain in fixer.topology.chains():
        lastIndexInChain = [i for i, res in enumerate(chain.residues())][-1]
        fixer.missingResidues[(chain.index, lastIndexInChain + 1)] = ["NME"]
        fixer.missingResidues[(chain.index, 0)] = ["ACE"]

    fixer.findNonstandardResidues()  # find non-standard residue
    fixer.replaceNonstandardResidues()  # replace non-standard residues with standard one
    fixer.findMissingAtoms()  # find missing heavy atoms
    fixer.addMissingAtoms()  # add missing atoms and residues
    fixer.addMissingHydrogens(ph)  # add missing hydrogens
    return fixer

def orient_molecule(modeller):
    # Extract the positions as a numpy array
    positions = np.array([[atom.x, atom.y, atom.z] for atom in modeller.positions])

    # Find the indices of the N-terminal C atom and C-terminal N atom
    n_terminal_c_atom_index = 14
    c_terminal_n_atom_index = 96

    # Get the positions of the N-terminal C and C-terminal N atoms
    c_atom_position = positions[n_terminal_c_atom_index]
    n_atom_position = positions[c_terminal_n_atom_index]

    # Translate the molecule to place the N-terminal C at the origin
    translation_vector = -c_atom_position
    positions += translation_vector

    # Align the C-terminal N atom along the z-axis
    n_atom_position_translated = positions[c_terminal_n_atom_index]
    z_axis = np.array([0.0, 0.0, 1.0])
    norm_vector = n_atom_position_translated / np.linalg.norm(n_atom_position_translated)
    axis = np.cross(norm_vector, z_axis)
    angle = np.arccos(np.clip(np.dot(norm_vector, z_axis), -1.0, 1.0))

    # Rodrigues' rotation formula
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    rotation_matrix = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

    # Apply rotation
    rotated_positions = np.dot(positions, rotation_matrix.T) * 10.0

    # Convert numpy array back to a list of Vec3 for OpenMM
    modeller.positions = [mm.Vec3(x, y, z) for x, y, z in rotated_positions]

# Load the PDB structure
fixer = prepare_protein(pdbfile, ignore_missing_residues=False, ph=7.0)

# Now you can continue with saving the system
#with open('oriented_system.pdb', 'w') as f:
#    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

# And also save the system configuration as before
#with open('system.xml', 'w') as f:
#    f.write(mm.openmm.XmlSerializer.serialize(system))

# Solvate
modeller = app.Modeller(fixer.topology, fixer.positions)
forcefield = app.ForceField('amber14-all.xml')
#forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
#modeller.addSolvent(forcefield, padding=1.0 * unit.nanometers, ionicStrength=0.15 * unit.molar)

# Apply the orientation function
orient_molecule(modeller)

# Save topology and positions
with open('system.pdb', 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

# System Configuration
nonbondedMethod = app.NoCutoff
nonbondedCutoff = 3.0*unit.nanometers
constraints = None
rigidWater = True
hydrogenMass = 1.0*unit.amu

# Create system
system = forcefield.createSystem(modeller.topology,
                                 nonbondedMethod=nonbondedMethod,
                                 nonbondedCutoff=nonbondedCutoff,
                                 constraints=constraints)

# Save system
with open('system.xml', 'w') as f:
    f.write(mm.openmm.XmlSerializer.serialize(system))

