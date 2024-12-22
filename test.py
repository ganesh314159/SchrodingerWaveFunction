import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("cppsol/psi-osc.dat")

# print(data[0][1])

x_data = []
y_data = []

for i in range(len(data)):
    x_data.append(data[i][0])
    y_data.append(data[i][1])

plt.plot(x_data, y_data)
plt.title("Schrodinger's Wavefunction")
plt.xlabel("X")
plt.ylabel("Psi")
plt.savefig(f"psi-osc.jpg")
# plt.show()


