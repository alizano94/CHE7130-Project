from ase import Atoms
from ase.calculators.emt import EMT
from ase.visualize import view
from ase.constraints import FixBondLengths
from ase.calculators.tip3p import TIP3P, rOH, angleHOH
from ase.md import Langevin
import ase.units as units
from ase.io.trajectory import Trajectory
from ase.build import molecule
from ase.build.attach import attach
import numpy as np

vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
WB_traj = Trajectory('tip3p_27mol_equil.traj')
water_box = WB_traj[-1]
water_box.set_constraint()  # repeat not compatible with FixBondLengths currently.
#water_box = water_box.repeat((2, 2, 2))
water_box.constraints = FixBondLengths([(3 * i + j, 3 * i + (j + 1) % 3)
                                        for i in range(int(len(water_box)/3))
                                        for j in [0, 1, 2]])
#view(water_box)


Oxygen_box = molecule('O2')

vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
Oxygen_box.set_cell((vol, vol, vol))
Oxygen_box.center()

#atoms = atoms.repeat((2, 2, 2))
Oxygen_box.set_pbc(True)
#print(int(len(atoms)/2))

Oxygen_box.constraints = FixBondLengths([(2*i , 2*i+1) 
                                        for i in range(int(len(Oxygen_box)/2))])

atoms = attach(water_box,Oxygen_box,vol,maxiter=100000)
#atoms = water_box + Oxygen_box 
view(atoms)