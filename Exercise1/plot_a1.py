import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":
    data = np.genfromtxt("../output/first_derivative.txt", delimiter="\t")

    # Plot only every 20th datapoint, to not overcrowd the plot
    x = data[::20, 0::20]
    y = data[::20, 1::20]

    # Plotting stuff...
    plt.plot(x, 2 * x, color="grey", zorder=1)
    plt.scatter(x, y, color="darkred", marker="x", zorder=5)
    plt.xlim((1e-10, 1e32))
    plt.ylim((1e-10, 2e32))
    plt.title("First derivative approximation of $x^2$")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(["$2x$ (exact)", "Approximation"])
    plt.tight_layout()
    plt.savefig("first_derivative_plot.pdf")
    plt.show()
