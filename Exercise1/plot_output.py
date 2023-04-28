import numpy as np
from matplotlib import pyplot as plt

data = np.genfromtxt("output.txt", delimiter="\t")

x = data[:, 0]
y = data[:, 1]

print(2*x-y)

plt.plot(x, 2*x)
plt.scatter(x, y)
plt.xlim((0, 1))
plt.ylim((0, 1))
plt.tight_layout()
plt.show()