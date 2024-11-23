import pdbfixer
import openmm as mm
import openmm.app as app
from openmm import unit
import sys
from math import ceil

# Load the topology and positions
with open('system.pdb', 'r') as f:
    pdb = app.PDBFile(f)

# Load the system
with open('system.xml', 'r') as f:
    system = mm.openmm.XmlSerializer.deserialize(f.read())

# Integrator options
constraintTolerance = 0.000001
dt = 0.004*unit.picoseconds
temperature = 300*unit.kelvin
friction = 1.0/unit.picosecond
pressure = 1.0*unit.atmospheres
barostatInterval = 25

# Simulation options
steps = 10*int(ceil(1e6 / dt.value_in_unit(unit.picoseconds))) # 10 microseconds
#steps = 100000 # 400 picoseconds
n_steps_save  = int(ceil(1e2 / dt.value_in_unit(unit.picoseconds))) # 0.1 ns
# Linux with GPU
platform = mm.openmm.Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}
# Mac
#platform = mm.openmm.Platform.getPlatformByName('OpenCL')
#platformProperties = {'Precision': 'single'}
dcdReporter = app.DCDReporter('3_production.dcd', n_steps_save)
dataReporter = app.StateDataReporter(sys.stdout, n_steps_save, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, elapsedTime=True, separator='\t')
#checkpointReporter = CheckpointReporter('run1_checkpoint.chk', 10000)

# Prepare simulation object
print('Building system...')
system.addForce(mm.openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = mm.openmm.LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
integrator.setRandomNumberSeed(283)
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

simulation.saveState("3_production.xml")

