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



O2_traj = Trajectory('O2_EMT_4mol_equil.traj')
Oxygen_box = O2_traj[-1]
Oxygen_box.set_constraint()  # repeat not compatible with FixBondLengths currently.



Oxygen_box.translate((0,-2.0*vol,0))
Oxygen_box.rotate(90, 'x')
Oxygen_box.translate((0,2.0*vol,0.5*vol))

atoms = water_box + Oxygen_box 
atoms.center()
atoms.set_pbc(True)
fa = FixBondLengths([(3 * i + j, 3 * i + (j + 1) % 3)
                    for i in range(int(len(water_box)/3))
                    for j in [0, 1, 2]])

fb = FixBondLengths([(2*i+81 , 2*i+82) for i in range(int(len(Oxygen_box)/2))])

atoms.set_constraint([fa, fb])

view(atoms)
tag = 'O2_WB_EMT_4mol_27mol_equil-2'
atoms.calc = EMT()
md = Langevin(atoms, 1 * units.fs, temperature_K=300,
              friction=0.01, logfile=tag + '.log')

traj = Trajectory(tag + '.traj', 'w', atoms)
md.attach(traj.write, interval=1)
md.run(5000)






