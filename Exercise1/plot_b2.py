import numpy as np
from matplotlib import pyplot as plt


if __name__ == "__main__":
    data = np.genfromtxt("../output/test.txt", delimiter="\t")

    # Plot only every 20th datapoint, to not overcrowd the plot
    x = data[::20, 0::20]
    y = data[::20, 1::20]

    # Plotting stuff...
    plt.scatter(x, y, color="darkred", marker="x", zorder=5)
    plt.xlim((0, 1))
    # plt.title("")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("b2_plot.pdf")
    plt.show()