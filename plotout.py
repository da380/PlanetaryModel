import numpy as np
import matplotlib.pyplot as plt

# d = np.loadtxt("Lagrange.out", delimiter=";")
d = np.loadtxt("SpecOutput.out", delimiter=";")
# d2 = np.loadtxt('time_test.STAT1.Z')

fig, ax = plt.subplots(2)
# fig  =plt.figure()
# ax = fig.add_subplot(111)
ax[0].plot(d[:, 0], d[:, 1], "k")
# Initialise the subplot function using number of rows and columns
# figure, axis = plt.subplots(1, 2)
# plt.plot(d[:, 0], d[:, 1], "k")

ax[0].plot(d[:,0],d[:,2],'r--')
ax[0].plot(d[:, 0], d[:, 1], "ko")
ax[0].legend(['Exact', 'Spectral element result', 'Nodes'])

ax[1].plot(d[:,0], abs(d[:,1] - d[:,2]))
ax[1].legend(['Absolute error'])

# plt.ylim(-0.00001, 0.00001)
plt.show()
