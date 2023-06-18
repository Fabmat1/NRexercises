import numpy as np
from matplotlib import pyplot as plt


A = np.genfromtxt("../output/A.csv", delimiter=",")
r = np.genfromtxt("../output/r.csv", delimiter=",")


# plt.ylim(0, 10e-6)
plt.yscale("log")
plt.plot(r[0, :], A[1, :])
plt.plot(r[0, :], A[2, :])
plt.plot(r[0, :], A[3, :])
plt.plot(r[0, :], A[4, :])
plt.tight_layout()
plt.show()


