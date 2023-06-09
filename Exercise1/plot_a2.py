import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":
    data = np.genfromtxt("../output/second_derivative.txt", delimiter="\t")

    x = data[:, 0]
    y = data[:, 1]

    # Plotting stuff...
    plt.plot(x, np.full(y.shape, 2), color="grey", zorder=1)
    plt.scatter(x, y, color="darkred", marker="x", zorder=5)
    plt.title("Second derivative approximation of $x^2$")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(["$2$ (exact)", "Approximation"])
    plt.tight_layout()
    plt.savefig("first_derivative_plot.pdf")
    plt.show()
