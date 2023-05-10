import numpy as np
from matplotlib import pyplot as plt


def u_prime(x):
    return 4*np.pi*np.sin(2*np.pi*x)*np.exp(-2 * np.cos(2 * np.pi * x))


if __name__ == "__main__":
    data = np.genfromtxt("../output/b2_data.txt", delimiter="\t")

    x = data[:, 0]
    y = data[:, 1]

    # Calculate chi-squared values
    chi_squared = ((u_prime(x) - y) ** 2)/y

    # Plotting
    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    # Main plot
    ax1.scatter(x[::20], y[::20], color="darkred", marker="x", zorder=5)
    ax1.plot(x, u_prime(x), color="grey", zorder=1)
    ax1.set_xlim((0, 1))
    ax1.set_title("Derivative approximation")
    ax1.set_ylabel("y")
    ax1.legend(["Approximation", "Analytic solution"])

    # Chi-squared scatter plot
    ax2.scatter(x, chi_squared, s=1, color="navy", marker="o", zorder=5)
    ax2.set_xlabel("x")
    ax2.set_ylabel("$\chi^2$")
    ax2.grid(True)

    # Adjust spacing between subplots
    plt.subplots_adjust(hspace=0)

    plt.tight_layout()
    plt.savefig("b2_plot.pdf")
    plt.show()
