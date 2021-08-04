from ase import Atoms
from ase.constraints import FixBondLengths
from ase.calculators.emt import EMT
from ase.md import Langevin
import ase.units as units
from ase.io.trajectory import Trajectory
import numpy as np
from ase.visualize import view
from ase.build import molecule


# Set up water box at 20 deg C density
atoms = molecule('O2')

vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
atoms.set_cell((vol, vol, vol))
atoms.center()

atoms = atoms.repeat((2, 2, 1))
atoms.set_pbc(True)
#print(int(len(atoms)/2))

view(atoms)

atoms.constraints = FixBondLengths([(2*i , 2*i+1) for i in range(int(len(atoms)/2))])

tag = 'O2_EMT_4mol_equil'
atoms.calc = EMT()
md = Langevin(atoms, 1 * units.fs, temperature_K=300,
              friction=0.01, logfile=tag + '.log')

traj = Trajectory(tag + '.traj', 'w', atoms)
md.attach(traj.write, interval=1)
md.run(5000)
