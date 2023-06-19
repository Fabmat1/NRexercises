import numpy as np
from matplotlib import pyplot as plt


A = np.genfromtxt("../output/A.csv", delimiter=",")
r = np.genfromtxt("../output/r.csv", delimiter=",")


# plt.yscale("log")
plt.ylim(0, 30)
plt.xlim(0, 10)

plt.plot(r[0, :], A[1, :])
plt.plot(r[0, :], A[2, :])
plt.plot(r[0, :], A[3, :])
plt.plot(r[0, :], A[4, :])
plt.plot(r[0, :], A[5, :])
plt.plot(r[0, :], A[6, :])
plt.plot(r[0, :], A[7, :])
plt.plot(r[0, :], A[8, :])
plt.plot(r[0, :], A[9, :])
plt.tight_layout()
plt.show()


