import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

t = np.linspace(-1e10, 1e10, 100)
r = np.linspace(0, np.pi / 2, 30)


# Calculate penrose coords from r and t (infinity works too)
def penrose_coords(r, t):
    R = 1. / np.pi * (np.arctan(t + r) - np.arctan(t - r))
    T = 1. / np.pi * (np.arctan(t + r) + np.arctan(t - r))
    # T = T[R >= 0]
    # R = R[R >= 0]
    return R, T


fig, ax = plt.subplots(1, 1)

skeleton = mpl.patches.Polygon([(0, 1), (1, 0), (0, -1), (-1, 0)], closed=True, label="_nolegend_", fc="white",
                               ec="black", linewidth=1)

ax.axis("off")

ax.add_patch(skeleton)
ax.plot(*penrose_coords(0, t), color="orange")  # Line along which r=0
ax.plot(*penrose_coords(t, t), color="red")  # Line with t to infinity, r to infinity and constant t-r
ax.plot(*penrose_coords(-t, t), color="green")  # Line with t to negative infinity, r to infinity and constant t+r
ax.scatter(*penrose_coords(np.inf, 0), marker="x", color="crimson", zorder=5)  # t = const, r = inf
ax.scatter(*penrose_coords(0, -np.inf), marker="x", color="lime", zorder=5)  # t = -inf, r = const
ax.scatter(*penrose_coords(0, np.inf), marker="x", color="navy", zorder=5)  # t = inf, r = const
plt.legend(["$r=0$",
            r"$t\rightarrow\infty, r\rightarrow\infty, t-r=$const.",
            r"$t\rightarrow-\infty, r\rightarrow\infty, t+r=$const.",
            "$t=$const., $r=\infty$",
            "$t=-\infty$, $r=$const.",
            "$t=\infty$, $r=$const."], loc="upper left")
plt.tight_layout()
plt.show()
