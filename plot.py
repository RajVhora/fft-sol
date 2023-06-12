import numpy as np
import matplotlib.pyplot as plt

#Print to console that graphing has began
print("Graphing has begun")

# Import data from numerical.dat file
data = np.loadtxt('numerical.dat')

radius = 50

# Plot data
plt.plot(data[:,0]/radius, data[:,1], label='sigma_11')
plt.plot(data[:,0]/radius, data[:,2], label='sigma_22')
plt.plot(data[:,0]/radius, data[:,3], label='sigma_12')
plt.xlabel('x/a')
plt.ylabel('sigma')
plt.legend()

plt.savefig('graph.png')

#Print to console that graphing has ended
print("Graphing has ended")