from ase import Atoms
from ase.io.trajectory import Trajectory
import numpy as np
import matplotlib.pyplot as plt

MSD = np.zeros((5000,2))

traj = Trajectory('O2_WB_EMT_4mol_27mol_equil-2.traj')

init_atoms = traj[0]
init_pos = init_atoms.get_positions()

for i in range(1,5000):   #each step
	atoms = traj[i]
	pos = atoms.get_positions()
	MSD[i][0] = i
	for j in range(81,89):  #each atom
		r_k = 0.0
		r_0 = 0.0
		for k in range(0,3): #each dimension
			r_k += pos[j:][0][k]**2
			r_0 += init_pos[j:][0][k]**2
		MSD[i][1] += (r_k**0.5-r_0**0.5)**2
	MSD[i][1] = MSD[i][1]/8/800

fit = np.polyfit(MSD[4200:4350,0],MSD[4200:4350,1],1)
print(fit)

plot_flag = True
if plot_flag:
	plt.plot(MSD[4200:,0],MSD[4200:,1])
	plt.xlabel("time (fs)")
	plt.ylabel("MSD ($\AA^2$)")
	plt.show()

