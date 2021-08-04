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
from ase.md.analysis import DiffusionCoefficient


traj = Trajectory('O2_WB_EMT_4mol_27mol_equil-2.traj')
dif = DiffusionCoefficient(traj,1 * units.fs,molecule=True)

dif.calculate(ignore_n_images=2000, number_of_segments=10)

D,std = dif.get_diffusion_coefficients()

print(D)
print(std)

d = D[0] *units.fs /10

print(d)

dif.plot(show=True)
