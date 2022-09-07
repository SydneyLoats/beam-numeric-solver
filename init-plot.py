#Sydney Loats SML4267

import matplotlib.pyplot as plt
import numpy as np
import math


x = np. linspace(0, 3, 100)
y = -0.1041666667*pow(x, 4) + 0.625*pow(x, 3) -2.8125 * x
plt.ylabel('displacement (m)')

plt.plot(x, y)

plt.xlabel('length (m)')
plt.title("Beam Displacement Analytic Solution")
plt.grid()
#plt.legend(['y = -37500x^4 + 225000x^3 - 1012500x'])
plt.legend(['y = -0.1042x^4 + 0.625x^3 - 2.8125x'])

plt.savefig('init-plot.png')
