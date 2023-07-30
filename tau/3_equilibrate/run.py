from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np
import sys

# Load the system

psf = CharmmPsfFile('../1_setup/5_ionize/ionized.psf')
pdb = PDBFile('../1_setup/5_ionize/ionized.pdb')

xyz = np.array(pdb.positions/nanometer)
xyz[:,0] -= np.amin(xyz[:,0])
xyz[:,1] -= np.amin(xyz[:,1])
xyz[:,2] -= np.amin(xyz[:,2])
pdb.positions = xyz*nanometer

with open('../2_equilibrate/system.xml', 'r') as f:
    system = openmm.XmlSerializer.deserialize(f.read())

# Integration Options

constraintTolerance = 0.000001
dt = 0.004*picoseconds
temperature = 300*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres
barostatInterval = 25

# Simulation Options

steps = 250000
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}
dcdReporter = DCDReporter('run.dcd', 2500)
dataReporter = StateDataReporter(sys.stdout, 2500, totalSteps=steps,
    step=True, time=True, speed=True, progress=True, elapsedTime=True,
    remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, volume=True, density=True, separator='\t')
#checkpointReporter = CheckpointReporter('checkpoint.chk', 10000)

# Prepare the Simulation

system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

# Positional restraints

restraint = openmm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
system.addForce(restraint)
restraint.addPerParticleParameter('k')
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')
for atom in psf.topology.atoms():
    if atom.name == 'CA' and atom.residue.chain.id != 'PROI':
        restraint.addParticle(atom.index, [1.0*kilocalories_per_mole/(angstrom**2), 
              pdb.positions[atom.index][0], pdb.positions[atom.index][1], pdb.positions[atom.index][2]] )

# COM restraint

indices_PROI = [atom.index for atom in psf.topology.atoms() if atom.name == 'CA' and atom.residue.chain.id == 'PROI']
z_coords = pdb.positions[indices_PROI]._value[:, 2]
z0 = np.mean(z_coords) * nanometer
com = CustomCentroidBondForce(1, 'k*(z1-z0)^2')
k = 73.0 * kilocalories_per_mole / (angstrom**2)
com.addGlobalParameter('k', k)
com.addGlobalParameter('z0', z0)
group_index = com.addGroup(indices_PROI)
com.addBond([group_index], [])
system.addForce(com)

# Create simulation object

integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(pdb.topology, system, integrator, platform, platformProperties)
#simulation.context.setPositions(pdb.positions)

# Load the state

simulation.loadState('../2_equilibrate/run.xml')

# Equilibrate

simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)

print('Continued equilibrating...')
#simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(steps)

# Simulate

#print('Simulating...')
#simulation.reporters.append(dcdReporter)
#simulation.reporters.append(dataReporter)
#simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0
#simulation.step(steps)

# Write file with final simulation state

simulation.saveState("run.xml")

