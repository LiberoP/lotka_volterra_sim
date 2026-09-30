# python script that only visualizes the Lotka-Volterra simulation output (results.txt).

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "figure.dpi": 120,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "legend.fontsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

data = np.loadtxt("./results.txt", skiprows=1)
t, x, y, H = data.T

fig, ax = plt.subplots(3, 1, figsize=(10, 14))

# 1. Time series: prey & predator
ax[0].plot(t, x, label="prey  x(t)", color="tab:blue", lw=1.5)
ax[0].plot(t, y, label="predator  y(t)", color="tab:red", lw=1.5)
ax[0].set(xlabel="t", ylabel="density", title="Populations over time")
ax[0].legend(loc="upper left", frameon=True, facecolor="white", framealpha=1.0,
             edgecolor="0.8")

# 2. Phase space
ax[1].plot(x, y, color="tab:purple", lw=1.5)
ax[1].plot(x[0], y[0], "o", color="k", ms=6, label="start")
ax[1].set(xlabel="prey  x(t)", ylabel="predator  y(t)",
          title="Phase space")
ax[1].legend(loc="upper left", frameon=True, facecolor="white", framealpha=1.0,
             edgecolor="0.8")

# 3. Relative drift of H
drift = 100 * (H - H[0]) / H[0]
ax[2].plot(t, drift, color="tab:orange", lw=1.5)
ax[2].set(xlabel="t", ylabel="(H − H₀) / H₀  (%)",
          title="Relative energy drift")

plt.suptitle("Lotka–Volterra, explicit Euler", fontsize=18)
plt.tight_layout()
plt.savefig("./results.png", dpi=150, bbox_inches="tight")
plt.close()