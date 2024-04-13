import pdbfixer
import openmm as mm
import openmm.app as app
from openmm import unit
import sys

# Load the topology and positions
with open('system.pdb', 'r') as f:
    pdb = app.PDBFile(f)

# Load the system
with open('system.xml', 'r') as f:
    system = mm.openmm.XmlSerializer.deserialize(f.read())

# Integrator options
dt = 0.001*unit.picoseconds
temperature = 300*unit.kelvin
friction = 1.0/unit.picosecond

# Simulation options
steps = 10000

# Linux with GPU
#platform = mm.openmm.Platform.getPlatformByName('CUDA')
#platformProperties = {'Precision': 'mixed'}

# Mac
platform = mm.openmm.Platform.getPlatformByName('OpenCL')
platformProperties = {'Precision': 'single'}

dcdReporter = app.DCDReporter('2_equilibration.dcd', 10)
dataReporter = app.StateDataReporter(sys.stdout, 10, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
#checkpointReporter = app.CheckpointReporter('2_equilibration.chk', 10000)

# Prepare simulation object
print('Building system...')
integrator = mm.openmm.LangevinMiddleIntegrator(temperature, friction, dt)
#integrator.setConstraintTolerance(constraintTolerance)

# Positional restraints
restraint = mm.openmm.CustomExternalForce('k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
system.addForce(restraint)
#restraint.addGlobalParameter('k', 10.0*unit.kilojoules_per_mole/(unit.nanometer*unit.nanometer))
restraint.addPerParticleParameter('k')
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')
for atom in pdb.topology.atoms():
    if atom.name == 'CA':
        restraint.addParticle(atom.index, [100 * unit.kilojoule_per_mole/unit.angstroms**2, 
              pdb.positions[atom.index][0], pdb.positions[atom.index][1], pdb.positions[atom.index][2]] )

simulation = app.Simulation(pdb.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(pdb.positions)

# Minimize
print('Performing energy minimization...')
simulation.minimizeEnergy()

# Equilibrate
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.step(steps)

#system.removeForce(system.getNumForces() - 1)
simulation.saveState("2_equilibration.xml")

