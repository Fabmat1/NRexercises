import numpy as np
from matplotlib import pyplot as plt


A = np.genfromtxt("../output/A.csv", delimiter=",")
r = np.genfromtxt("../output/r.csv", delimiter=",")


# plt.ylim(0, 10e-6)
plt.plot(r[0, :], A[2, :])
plt.tight_layout()
plt.show()


