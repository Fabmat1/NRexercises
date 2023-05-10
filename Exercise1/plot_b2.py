import numpy as np
from matplotlib import pyplot as plt


def u(x):
    return np.exp(-2 * np.cos(2 * np.pi * x))


if __name__ == "__main__":
    data = np.genfromtxt("../output/b2_data.txt", delimiter="\t")

    # Plot only every 20th datapoint, to not overcrowd the plot
    x = data[:, 0]
    y = data[:, 1]

    # Plotting stuff...
    plt.plot(x, y, color="darkred", zorder=5)
    plt.plot(x, u(x + 15000 / 100000), color="lightgrey", zorder=1)
    plt.xlim((0, 1))
    # plt.title("")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("b2_plot.pdf")
    plt.show()
