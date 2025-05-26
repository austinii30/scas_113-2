import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import numpy as np

x = np.linspace(-3.0, 3.0, 100)
y = np.linspace(-2.0, 2.0, 100)
X, Y = np.meshgrid(x, y)

Z1 = np.exp(-(X)**2 - (Y)**2)
Z2 = np.exp(-(X * 10)**2 - (Y * 10)**2)
z = Z1 + 50 * Z2

fig, ax = plt.subplots()

n_levels = 50
cs = ax.contourf(X, Y, z, 
                 np.logspace(np.log10(z.min()),np.log10(z.max()), n_levels), 
                 locator=ticker.LogLocator(), 
                 cmap=cm.jet
                 )

cbar = fig.colorbar(cs)
cbar.locator = ticker.LogLocator(10)
cbar.set_ticks(cbar.locator.tick_values(z.min(), z.max()))
cbar.minorticks_off()

plt.savefig("my_density_plot.pdf", bbox_inches='tight')
