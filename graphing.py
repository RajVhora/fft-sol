import numpy as np
import matplotlib.pyplot as plt

# Import data from numerical.dat file
data = np.loadtxt('numerical.dat')

# Plot data
plt.plot(data[:,0], data[:,1], label='Numerical')