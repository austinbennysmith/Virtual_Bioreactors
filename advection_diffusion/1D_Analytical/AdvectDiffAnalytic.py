## Plotting analytical solutions to the advection-diffusion equation.
## This combines the two codes for advection & diffusion analytical solutions.

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import quad

U = 2.0

X = np.linspace(-5, 5, 1000)
def A0integrand(x):
	return (1/10)*np.exp(-np.square(x))
def Aintegrand(x, n):
	return (1/5)*np.exp(-np.square(x))*np.cos((n*np.pi*x)/5)
def Bintegrand(x, n):
	return (1/5)*np.exp(-np.square(x))*np.sin((n*np.pi*x)/5)

for t in range(10):
	# Getting A0 (has to be found separately from the other An)
	I0 = quad(A0integrand, -5, 5)
	A0 = I0[0]
	# Initializing the array that will be iteratively added to:
	Asum = A0*np.cos((0*np.pi*X)/5)*np.exp(-t*np.square((0*np.pi)/5))
	Bsum = np.zeros(1000)
	# print(Asum)

	# AnLIST = []

	for n in range(1, 20): # Ideally this would be done for infinite steps. However, the sum converges I think (An decreases monotonically, can be checked by plotting it) so it should be a dcent approximation for a sufficiently small number of steps
		IAn = quad(Aintegrand, -5, 5, args=(n))
		An = IAn[0]
		# print("An: ", An)

		IBn = quad(Bintegrand, -5, 5, args=(n))
		Bn = IBn[0]
		# print("Bn: ", Bn)
		# AnLIST.append(An)
		Asummand = An*np.cos((n*np.pi*X)/5)*np.exp(-t*np.square((n*np.pi)/5))
		Bsummand = Bn*np.sin((n*np.pi*X)/5)*np.exp(-t*np.square((n*np.pi)/5))
		
		Asum = Asum+Asummand
		Bsum = Bsum+Bsummand

		SUM = Asum+Bsum
	
	# Shiting the graph by the veolcity * time elapsed:
	distance = U*t
	Yleft = SUM[X>=max(X)-(10 - (distance % 10))]
	Yright = SUM[X<max(X)-(10 - (distance % 10))]
	Yall = np.concatenate((Yleft, Yright))

	plt.plot(X, Yall, label="t="+str(t))
plt.legend()
plt.title("Heat Equation Analytical Solution")
plt.xlabel("x")
plt.ylabel("T")
plt.show()