#Sydney Loats SML4267

import matplotlib.pyplot as plt
import numpy as np

point1 = [10, 0.001064]
point2 = [20, 0.008169]
point3 = [50, 0.123086]
point4 = [100, 1.37639]
point5 = [200, 19.4254]


x = [point1[0], point2[0], point3[0], point4[0], point5[0]]

y = [point1[1], point2[1], point3[1], point4[1], point5[1]]

plt.plot(x, y)

plt.plot(10, 0.001064, marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')
plt.plot(20, 0.008169, marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')
plt.plot(50, 0.123086, marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')
plt.plot(100, 1.37639, marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')
plt.plot(200, 19.4254, marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')

plt.xlabel('n')
plt.ylabel('seconds')
plt.title('Runtime Performance as a Function of n')
plt.grid()

plt.savefig('runtime-plot.png')
