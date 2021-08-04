from ase import Atoms
from ase.calculators.emt import EMT
from ase.visualize import view
from ase.constraints import FixBondLengths
from ase.calculators.tip3p import TIP3P, rOH, angleHOH
from ase.md import Langevin
import ase.units as units
from ase.io.trajectory import Trajectory
from ase.build import molecule
import numpy as np


traj = Trajectory('tip3p_27mol_equil.traj')
atoms = traj[-1]

# Repeat box and equilibrate further.
tag = 'tip3p_216mol_equil'
atoms.set_constraint()  # repeat not compatible with FixBondLengths currently.
atoms = atoms.repeat((2, 2, 2))
atoms.constraints = FixBondLengths([(3 * i + j, 3 * i + (j + 1) % 3)
    for i in range(int(len(atoms) / 3)) 
    for j in [0, 1, 2]])

OO = molecule('O2')
vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
OO.set_cell((vol, vol, vol))
OO = OO.repeat((2,2,2))
#atoms = atoms + OO


view(OO)

#atoms.calc = TIP3P(rc=7.)
#md = Langevin(atoms, 1 * units.fs, temperature_K=300,
#              friction=0.01, logfile=tag + '.log')

#traj = Trajectory(tag + '.traj', 'w', atoms)
#md.attach(traj.write, interval=1)
#md.run(2000)