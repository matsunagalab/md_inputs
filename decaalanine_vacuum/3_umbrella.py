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
steps = 100000000 # 100 ns

# Linux with GPU
#platform = mm.openmm.Platform.getPlatformByName('CUDA')
#platformProperties = {'Precision': 'mixed'}

# Mac
platform = mm.openmm.Platform.getPlatformByName('OpenCL')
platformProperties = {'Precision': 'single'}

dcdReporter = app.DCDReporter('3_umbrella.dcd', 1000)
dataReporter = app.StateDataReporter(sys.stdout, 1000, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
#checkpointReporter = app.CheckpointReporter('3_umbrella.chk', 10000)

# Prepare simulation object
print('Building system...')

## Positional restraints
k = 50.0 * unit.kilojoule_per_mole/unit.angstroms**2  # force constant
force = mm.CustomExternalForce("k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
force.addGlobalParameter("k", k)
force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")

c_atom_index = 14
n_atom_index = 96

force.addParticle(c_atom_index, [0.0 * unit.nanometers, 0.0 * unit.nanometers, 0.0 * unit.nanometers])
force.addParticle(n_atom_index, [0.0 * unit.nanometers, 0.0 * unit.nanometers, 1.4 * unit.nanometers])

system.addForce(force)

## integrator
integrator = mm.openmm.LangevinMiddleIntegrator(temperature, friction, dt)
#integrator.setConstraintTolerance(constraintTolerance)

## create simulation
simulation = app.Simulation(pdb.topology, system, integrator, platform, platformProperties)
#simulation.context.setPositions(pdb.positions)

# Load the state
simulation.loadState('2_equilibration.xml')

# Production
print('Production...')
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.currentStep = 0
simulation.step(steps)

simulation.saveState("3_umbrella.xml")

