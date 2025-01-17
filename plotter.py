'''
plotter.py --- Primary tester file.

This takes the data generated by main.py and turns it into our plots.

Attribution - Jimmy
'''

import numpy as np
import matplotlib.pyplot as plt


X, izzo = np.loadtxt('izzo2015.csv', delimiter=',', dtype=float)
X, arora = np.loadtxt('arora2013.csv', delimiter=',', dtype=float)
X, avanzini = np.loadtxt('avanzini2008.csv', delimiter=',', dtype=float)
X, gooding = np.loadtxt('gooding1990.csv', delimiter=',', dtype=float)

step_size = np.array([...])  # Array with step sizes (in days)
runtime1 = np.array([...])  # Array with runtime for function 1 (in seconds)
runtime2 = np.array([...])  # Array with runtime for function 2 (in seconds)
runtime3 = np.array([...])  # Array with runtime for function 3 (in seconds)
runtime4 = np.array([...])  # Array with runtime for function 4 (in seconds)

# Create two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot the data for the 4 functions on the first subplot (normal axes)
ax1.plot(X, izzo, label='Izzo')
ax1.plot(X, arora, label='Arora and Russel')
ax1.plot(X, avanzini, label='Avanzini')
ax1.plot(X, gooding, label='Gooding')

# Add labels and title for the first subplot
ax1.set_xlabel('Step Size (days)')
ax1.set_ylabel('Runtime (seconds)')
ax1.set_title('Runtime vs Step Size')
ax1.legend()

# Plot the data for the 4 functions on the second subplot (log scale on the y-axis)
ax2.plot(X, izzo, label='Izzo')
ax2.plot(X, arora, label='Arora and Russel')
ax2.plot(X, avanzini, label='Avanzini')
ax2.plot(X, gooding, label='Gooding')

# Set log scale for the y-axis, and add labels and title for the second subplot
ax2.set_yscale('log')
ax2.set_xlabel('Step Size (days)')
ax2.set_ylabel('Runtime (seconds, log scale)')
ax2.set_title('Runtime vs Step Size (Log Scale)')
ax2.legend()

# Display the plots
plt.tight_layout()
plt.savefig('resultsfig', dpi=400)
plt.show()
