import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
plt.switch_backend('agg')
Omega = (0, 1)
h = 0.1
n = 6
k = 3/2 * n
k_alias = 1/2 * n
x = np.linspace(Omega[0], Omega[1], 100)
w = np.sin(k * x * np.pi)
w_alias = np.sin(k_alias * x * np.pi)
fig, axs = plt.subplots(figsize=(8, 4))
axs.plot(x, w, label=f'k = {k}')
axs.plot(x, -w_alias, label=f'k\' = {k_alias}')
axs.set_xlim(Omega)
axs.set_xlabel('x')
axs.set_ylabel('w(x)')
axs.set_title('the case of n = 6 for Example 9.13')
axs.legend()
plt.savefig("./9_14.eps")