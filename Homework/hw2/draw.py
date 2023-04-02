import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
plt.switch_backend('agg')
# Define the interval and the discretization step size
Omega = (0, 1)
h = 0.1

# Define n and k values
n = 6
k = 3/2 * n
k_alias = 1/2 * n

# Create a range of x values to plot
x = np.linspace(Omega[0], Omega[1], 100)

# Compute the Fourier modes with and without aliasing
w = np.sin(k * x * np.pi)
w_alias = np.sin(k_alias * x * np.pi)

# Plot the results
fig, axs = plt.subplots(figsize=(8, 4))
axs.plot(x, w, label=f'k = {k}')
axs.plot(x, -w_alias, label=f'k\' = {k_alias}')
axs.set_xlim(Omega)
axs.set_xlabel('x')
axs.set_ylabel('w(x)')
axs.set_title(f'Aliased Fourier modes for n = {n}')
axs.legend()
plt.savefig("./9_14.eps")