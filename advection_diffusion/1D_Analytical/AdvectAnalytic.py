import numpy as np 
import matplotlib.pyplot as plt

U = 2.0

X = np.linspace(-5, 5, 1000)
Yanalytical = np.exp(-X*X)
for t in range(10):
	distance = U*t
	Yleft = Yanalytical[X>=max(X)-(10 - (distance % 10))]
	Yright = Yanalytical[X<max(X)-(10 - (distance % 10))]
	Yall = np.concatenate((Yleft, Yright))
	plt.plot(X, Yall, label="t="+str(t))
plt.legend()
plt.title("Advection Analytical Solution")
plt.xlabel("x")
plt.ylabel("T")
plt.show()