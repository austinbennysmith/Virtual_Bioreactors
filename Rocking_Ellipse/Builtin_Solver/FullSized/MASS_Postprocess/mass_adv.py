import numpy as np
import matplotlib.pyplot as plt 

mass_total = np.genfromtxt("FullSized_advectionh/mass.dat")
mass_IN = np.genfromtxt("FullSized_advectionh/mass_IN.dat")
mass_OUT = np.genfromtxt("FullSized_advectionh/mass_OUT.dat")
mass_water = np.genfromtxt("FullSized_advectionh/mass_water.dat")
mass_air = np.genfromtxt("FullSized_advectionh/mass_air.dat")
# print(mass_total)
# print("hello world")
mass_total = mass_total[:, [1, 4]] # Python automatically splits up the columns where the are spaces in the table. I'm picking the columns with numbers in them.
mass_IN = mass_IN[:, [1, 4]]
mass_OUT = mass_OUT[:, [1, 4]]
mass_water = mass_water[:, [1, 4]]
mass_air = mass_air[:, [1, 4]]

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(20, 20))
axlist = ax.flatten()

# plt.plot(mass_total[:, 0], mass_total[:, 1])

axlist[0].plot(mass_IN[:, 0], mass_IN[:, 1])
axlist[0].set_xlabel('Dimensionless time')
axlist[0].set_ylabel('Tracer Mass')
axlist[0].set_title('Tracer Mass Inside the Reactor')

axlist[1].plot(mass_OUT[:, 0], mass_OUT[:, 1])
axlist[1].set_xlabel('Dimensionless time')
axlist[1].set_ylabel('Tracer Mass')
axlist[1].set_title('Tracer Mass Outside the Reactor')

axlist[2].plot(mass_water[:, 0], mass_water[:, 1])
axlist[2].set_xlabel('Dimensionless time')
axlist[2].set_ylabel('Tracer Mass')
axlist[2].set_title('Tracer Mass in the Water')

axlist[3].plot(mass_air[:, 0], mass_air[:, 1])
axlist[3].set_xlabel('Dimensionless time')
axlist[3].set_ylabel('Tracer Mass')
axlist[3].set_title('Tracer Mass in the Air')

plt.savefig('tracer_masses_adv.png')
plt.show()