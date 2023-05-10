import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":
    # Plot the example RMSE growing for the Euler method as time increases
    data = np.genfromtxt("../output/RMSE_RK4.txt", delimiter="\t")

    x = data[:, 0]
    y = data[:, 1]

    # Plotting stuff...
    plt.plot(x, y, color="darkred", zorder=5)
    plt.xlim((0, 1))
    plt.title("RMSE over time with the RK4 method")
    plt.xlabel("t")
    plt.ylabel("RMSE")
    plt.tight_layout()
    plt.savefig("b3_RMSE.pdf")
    plt.show()

    # Plot heatmap of RMSE for different delta_t and delta_x
    data = np.genfromtxt("../output/RMSEspace_RK4.txt", delimiter="\t")
    params = np.genfromtxt("../output/RMSEspace_RK4.txt_params.txt", delimiter="\t")

    # Plotting stuff...
    plt.imshow(data, cmap="viridis", interpolation="nearest", aspect="auto", extent=[10, params[1], params[0], 10])
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('Average RSME', rotation=270)

    plt.gca().invert_yaxis()
    plt.title("RMSE heatmap for the RK4 method")
    plt.xlabel("$N_x$")
    plt.ylabel("$N_t$")
    plt.tight_layout()
    plt.savefig("b3_RMSE_heatmap_RK4.pdf")
    plt.show()
