#Sydney Loats SML4267

import matplotlib.pyplot as plt
import numpy as np
import math

point1 = [math.log(10), math.log(0.0101603)]
point2 = [math.log(20), math.log(0.00227957)]
point3 = [math.log(50), math.log(0.000340377)]
point4 = [math.log(100), math.log(0.0000740612)]
point5 = [math.log(200), math.log(0.000019358)]

x = [point1[0], point2[0], point3[0], point4[0], point5[0]]

y = [point1[1], point2[1], point3[1], point4[1], point5[1]]


x2 = np.linspace(1, 5, 5) 
y2 = -2*x2

plt.plot(x2, y2)

plt.plot(x, y)

plt.legend(['log(l2)xlog(n)', 'y = -2x'])
plt.xlabel('log(n)')
plt.ylabel('log(l2)')
plt.title('log(l2) x log(n)')
plt.grid()

#plt.plot(math.log(10), math.log(0.0101603), marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')
#plt.plot(math.log(20), math.log(0.00227957), marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')
#plt.plot(math.log(50), math.log(0.000340377), marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')
#plt.plot(math.log(100), math.log(0.0000740612), marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')
#plt.plot(math.log(200), math.log(0.000019358), marker = 'o', markersize = 5, markeredgecolor = 'red', markerfacecolor = 'red')

plt.savefig('log-plot.png')
