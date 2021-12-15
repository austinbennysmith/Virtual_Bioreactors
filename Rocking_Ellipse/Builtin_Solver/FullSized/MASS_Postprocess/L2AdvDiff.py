import numpy as np
import matplotlib.pyplot as plt 

# With just advection:
mass_total1 = np.genfromtxt("FullSized_advectionh/mass.dat")
mass_IN1 = np.genfromtxt("FullSized_advectionh/mass_IN.dat")
mass_OUT1 = np.genfromtxt("FullSized_advectionh/mass_OUT.dat")
mass_water1 = np.genfromtxt("FullSized_advectionh/mass_water.dat")
mass_air1 = np.genfromtxt("FullSized_advectionh/mass_air.dat")

# With diffusion and advection both:
mass_total2 = np.genfromtxt("FullSized_diffusionh/mass.dat")
mass_IN2 = np.genfromtxt("FullSized_diffusionh/mass_IN.dat")
mass_OUT2 = np.genfromtxt("FullSized_diffusionh/mass_OUT.dat")
mass_water2 = np.genfromtxt("FullSized_diffusionh/mass_water.dat")
mass_air2 = np.genfromtxt("FullSized_diffusionh/mass_air.dat")

mass_total1 = mass_total1[:, [1, 4]] # Python automatically splits up the columns where the are spaces in the table. I'm picking the columns with numbers in them.
mass_IN1 = mass_IN1[:, [1, 4]]
mass_OUT1 = mass_OUT1[:, [1, 4]]
mass_water1 = mass_water1[:, [1, 4]]
mass_air1 = mass_air1[:, [1, 4]]

mass_total2 = mass_total2[:, [1, 4]] # Python automatically splits up the columns where the are spaces in the table. I'm picking the columns with numbers in them.
mass_IN2 = mass_IN2[:, [1, 4]]
mass_OUT2 = mass_OUT2[:, [1, 4]]
mass_water2 = mass_water2[:, [1, 4]]
mass_air2 = mass_air2[:, [1, 4]]

fig, ax0 = plt.subplots(nrows=2, ncols=2, figsize=(20, 20))
axlist0 = ax0.flatten()

axlist0[0].plot(mass_IN1[:, 0], mass_IN1[:, 1], label="No diffusion")
axlist0[0].plot(mass_IN2[:, 0], mass_IN2[:, 1], label="Diffusion")
axlist0[0].set_xlabel('Dimensionless time')
axlist0[0].set_ylabel('Tracer Mass')
axlist0[0].set_title('Mass in Ellipse, advection.h vs diffusion.h')
axlist0[0].legend()

axlist0[1].plot(mass_OUT1[:, 0], mass_OUT1[:, 1], label="No diffusion")
axlist0[1].plot(mass_OUT2[:, 0], mass_OUT2[:, 1], label="Diffusion")
axlist0[1].set_xlabel('Dimensionless time')
axlist0[1].set_ylabel('Tracer Mass')
axlist0[1].set_title('Mass outside Ellipse, advection.h vs diffusion.h')
axlist0[1].legend()

axlist0[2].plot(mass_water1[:, 0], mass_water1[:, 1], label="No diffusion")
axlist0[2].plot(mass_water2[:, 0], mass_water2[:, 1], label="Diffusion")
axlist0[2].set_xlabel('Dimensionless time')
axlist0[2].set_ylabel('Tracer Mass')
axlist0[2].set_title('Mass in Water, advection.h vs diffusion.h')
axlist0[2].legend()

axlist0[3].plot(mass_air1[:, 0], mass_air1[:, 1], label="No diffusion")
axlist0[3].plot(mass_air2[:, 0], mass_air2[:, 1], label="Diffusion")
axlist0[3].set_xlabel('Dimensionless time')
axlist0[3].set_ylabel('Tracer Mass')
axlist0[3].set_title('Mass in Air, advection.h vs diffusion.h')
axlist0[3].legend()

plt.savefig('Comparison.png')
plt.show()

# print(mass_total1[:,1])
# firstarr = np.stack((mass_total1[:,1], mass_total2[:,1]), axis=0)

massList = []
inList = []
outList = []
waterList = []
airList = []

# print(len(mass_total1[:, 0]))
# print(mass_total1.shape)
# print(mass_total1[1,1])
for i in range(len(mass_total1)):
	massList.append(np.linalg.norm(mass_total1[i, 1]-mass_total2[i, 1]))
	inList.append(np.linalg.norm(mass_IN1[i, 1]-mass_IN2[i, 1]))
	outList.append(np.linalg.norm(mass_OUT1[i, 1]-mass_OUT2[i, 1]))
	waterList.append(np.linalg.norm(mass_water1[i, 1]-mass_water2[i, 1]))
	airList.append(np.linalg.norm(mass_air1[i, 1]-mass_air2[i, 1]))

massArr = np.array(massList)
inArr = np.array(inList)
outArr = np.array(outList)
waterArr = np.array(waterList)
airArr = np.array(airList)

t = mass_total1[:, 0]

massArr = np.stack((t, massArr), axis=1)
inArr = np.stack((t, inArr), axis=1)
outArr = np.stack((t, outArr), axis=1)
waterArr = np.stack((t, waterArr), axis=1)
airArr = np.stack((t, airArr), axis=1)

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(20, 20))
axlist = ax.flatten()

axlist[0].plot(inArr[:, 0], inArr[:, 1])
axlist[0].set_xlabel('Dimensionless time')
axlist[0].set_ylabel('L2 Norm')
axlist[0].set_title('Mass in Ellipse, advection.h vs diffusion.h')

axlist[1].plot(outArr[:, 0], outArr[:, 1])
axlist[1].set_xlabel('Dimensionless time')
axlist[1].set_ylabel('L2 Norm')
axlist[1].set_title('Mass outside Ellipse, advection.h vs diffusion.h')

axlist[2].plot(waterArr[:, 0], waterArr[:, 1])
axlist[2].set_xlabel('Dimensionless time')
axlist[2].set_ylabel('L2 Norm')
axlist[2].set_title('Mass in Water, advection.h vs diffusion.h')

axlist[3].plot(airArr[:, 0], airArr[:, 1])
axlist[3].set_xlabel('Dimensionless time')
axlist[3].set_ylabel('L2 Norm')
axlist[3].set_title('Mass in Air, advection.h vs diffusion.h')

plt.savefig('L2_norms.png')
plt.show()