import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


def create_animation(file, r_array, outname, dt, y_limits=None, x_limits=None, logarithmic=False, title=None, yaxis_name=None, *args, **kwargs):
    data = np.genfromtxt(file, delimiter=",")
    fig, ax = plt.subplots(figsize=(4.8 * 16 / 9, 4.8))

    def update(frame):
        ax.clear()

        if logarithmic:
            ax.set_yscale("symlog")
        y = data[frame, :]
        if title is not None:
            ax.set_title(title)

        timestamp = f"$t$ = {round(frame*dt, 2)}"

        while len(timestamp) < 10:
            timestamp += "0"

        ax.annotate(timestamp, xy=(1, 1), xycoords='axes fraction', xytext=(-10, -10),
                    textcoords='offset points', ha='right', va='top', fontname='monospace')

        if yaxis_name is not None:
            ax.set_ylabel(yaxis_name)

        ax.set_xlabel("Radius $r$")

        if y_limits is not None:
            plt.ylim(y_limits)
        if x_limits is not None:
            plt.xlim(x_limits)
        ax.plot(r_array, y, **kwargs)
        plt.tight_layout()

    extra_args = ['-s', '3840x2160', '-b:v', '2M']
    animation = FuncAnimation(fig, update, frames=data.shape[0], interval=40)

    plt.rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ffmpeg\bin\ffmpeg.exe'
    animation.save(outname, writer="ffmpeg", dpi=250, extra_args=extra_args)

    plt.close(fig)


A = np.genfromtxt("../output/A.csv", delimiter=",")
r = np.genfromtxt("../output/r.csv", delimiter=",")


# plt.yscale("log")
plt.ylim(0.9, 7)
plt.xlim(0, 10)

plt.plot(r[0, :], A[1, :])
plt.plot(r[0, :], A[10, :])
plt.plot(r[0, :], A[20, :])
plt.plot(r[0, :], A[30, :])
plt.plot(r[0, :], A[40, :])
plt.plot(r[0, :], A[50, :])
plt.plot(r[0, :], A[60, :])
plt.plot(r[0, :], A[70, :])
plt.plot(r[0, :], A[80, :])
plt.tight_layout()
plt.show()


r_bh = np.loadtxt("../output/r_BH_scalar.csv", delimiter=",")
t = np.linspace(0, 10, len(r_bh))
plt.plot(t, r_bh, color="crimson")
plt.title("Black hole apparent radius as a function of time.")
plt.ylabel("Radius $r_h$")
plt.xlabel("Time $t$")
plt.ylim(0.5, 2.5)
plt.xlim(0, 10)
plt.tight_layout()
plt.savefig("r_bh_over_time.pdf")
plt.show()

S_bh = np.loadtxt("../output/S_BH_apparent.csv", delimiter=",")
plt.plot(t, S_bh, color="crimson")
plt.title("Black hole apparent surface area as a function of time.")
plt.ylabel("Radius $r_h$")
plt.xlabel("Time $t$")
plt.xlim(0, 10)
plt.tight_layout()
plt.savefig("S_bh_over_time.pdf")
plt.show()


create_animation("../output/A.csv", r[0, :], "A_animation.mp4", 0.005*25, y_limits=(0, 10), x_limits=(0, 10), color="navy",
                 yaxis_name="Parameter $A$", title="Parameter $A$ plotted over $r$")
create_animation("../output/r_BH_apparent.csv", r[0, :], "r_BH_apparent_animation.mp4", 0.005*25, y_limits=(-3, 3),
                 x_limits=(0, 10), color="navy", yaxis_name="Black hole radius $r_h$",
                 title="Black hole radius plotted over $r$")
create_animation("../output/B.csv", r[0, :], "B_animation.mp4", 0.005*25, y_limits=(0, 1),
                 x_limits=(0, 10), color="navy", yaxis_name="Parameter $B$", title="Parameter $B$ plotted over $r$")
create_animation("../output/alpha.csv", r[0, :], "alpha_animation.mp4", 0.005*25, y_limits=(0, 5),
                 x_limits=(0, 10), color="navy", yaxis_name=r"Parameter $\alpha$", title=r"Parameter $\alpha$ plotted over $r$")


