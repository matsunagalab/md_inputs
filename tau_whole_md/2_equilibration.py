from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np
import sys

# Input Files

psf = CharmmPsfFile('./1_setup/5_ionize/ionized.psf')
pdb = PDBFile('./1_setup/5_ionize/ionized.pdb')
#crd = CharmmPdbFile('./1_setup/5_ionize/ionized.crd')
params = CharmmParameterSet('./toppar/top_all36_prot.rtf', './toppar/par_all36m_prot.prm', './toppar/toppar_water_ions.str')
#params = CharmmParameterSet('par_all36m_prot.prm')

xyz = np.array(pdb.positions/nanometer)
xyz[:,0] -= np.amin(xyz[:,0])
xyz[:,1] -= np.amin(xyz[:,1])
xyz[:,2] -= np.amin(xyz[:,2])
pdb.positions = xyz*nanometer

#psf.topology.setPeriodicBoxVectors(pdb.topology.getPeriodicBoxVectors())
psf.setBox(14.0*nanometer, 14.0*nanometer,  14.0*nanometer)

# System Configuration

nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometer
ewaldErrorTolerance = 0.0005
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001
hydrogenMass = 1.5*amu

# Integration Options

dt = 0.004*picoseconds
temperature = 300*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres
barostatInterval = 25

# Simulation Options

steps = 250000
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}
dcdReporter = DCDReporter('2_equilibration.dcd', 2500)
dataReporter = StateDataReporter(sys.stdout, 2500, totalSteps=steps,
    step=True, time=True, speed=True, progress=True, elapsedTime=True,
    remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, volume=True, density=True, separator='\t')
#checkpointReporter = CheckpointReporter('checkpoint.chk', 10000)

# Prepare the Simulation

print('Building system...')
topology = pdb.topology
positions = pdb.positions
system = psf.createSystem(params, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
    constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance, hydrogenMass=hydrogenMass)
with open('system.pdb', 'w') as f:
    PDBFile.writeFile(pdb.topology, pdb.positions, f)
with open("system.xml", mode="w") as file:
    file.write(XmlSerializer.serialize(system))
system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

# Positional restraints

restraint = openmm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
system.addForce(restraint)
#restraint.addGlobalParameter('k', 10.0*kilojoules_per_mole/(nanometer**2))
restraint.addPerParticleParameter('k')
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')
for atom in pdb.topology.atoms():
    if atom.name == 'CA':
        restraint.addParticle(atom.index, [10.0*kilocalories_per_mole/(angstrom**2), 
              pdb.positions[atom.index][0], pdb.positions[atom.index][1], pdb.positions[atom.index][2]] )

# Create simulation object

integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
integrator.setRandomNumberSeed(283)
simulation = Simulation(topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)
with open("integrator.xml", mode="w") as file:
    file.write(XmlSerializer.serialize(integrator))

# Minimize and Equilibrate

print('Performing energy minimization...')
simulation.minimizeEnergy()

print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.step(steps)

#system.removeForce(system.getNumForces() - 1)
simulation.saveState("2_equilibration.xml")

