import numpy as np 
import matplotlib.pyplot as plt 

FILE = input("FILENAME?")
myarray = np.genfromtxt(FILE, delimiter=' ')

# Below I figure out how many cells in the domain:
mask = myarray[:, 2]==0
domainSize = len(myarray[mask, :])
dt = myarray[domainSize+1, 3]
print(dt)

# fig, axarr = plt.subplots(nrows=2, ncols=5, figsize=(40, 40))
# axlist = axarr.flatten()
for i in range(0, int(myarray.shape[0]/domainSize), int(1/dt)): # If doing fig, ax, either reorganize somehow or subtract 1 from the second argument
	# interf=interf[interf[:, 0].argsort()] # When doing parallel processing, the output gets jumbled. This line sorts the output file so that the x values are in order, which makes the plots look correct.
	myarrayNOW = myarray[i*domainSize:(i+1)*domainSize, :]
	# print(interf)
	

	# axlist[int(i)].plot(myarrayNOW[:,  0], myarrayNOW[:, 1], label="t="+str(myarrayNOW[0, 2]))
	# axlist[int(i)].legend()
	# axlist[int(i)].set_xlabel("x")
	# axlist[int(i)].set_ylabel("T")

	plt.plot(myarrayNOW[:,  0], myarrayNOW[:, 1], label="t="+str(myarrayNOW[0, 3]))
	plt.legend()
	plt.xlabel("x")
	plt.ylabel("T")
	# if i==100:
	# 	break
plt.show()