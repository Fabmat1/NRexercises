import numpy as np
from matplotlib import pyplot as plt

data = np.genfromtxt("output.txt", delimiter="\t")

x = data[:, 0]
y = data[:, 1]

plt.scatter(x, y)
plt.tight_layout()
plt.show()
