
import numpy as np
import matplotlib.pyplot as plt
import math

m = 940
k = 100
hbar = 197.3

N = 400
rmax = 5
dr = rmax / (N + 1)

r = np.linspace(dr, rmax-dr, N)

for b in [0, 10, 30]:

    V = 0.5*k*r**2 + (1/3)*b*r**3

    d = hbar**2/(m*dr**2) + V
    o = -hbar**2/(2*m*dr**2)

    H = np.diag(d) + np.diag(np.full(N-1, o), 1)
    H += np.diag(np.full(N-1, o), -1)

    E, U = np.linalg.eigh(H)

    energy = E[0]
    psi = U[:, 0]

    psi = psi / np.sqrt(np.sum(psi**2)*dr)

    print("b =", b, "Ground State Energy =", round(energy, 3), "MeV")

    plt.plot(r, psi, label="b = " + str(b))

plt.xlabel("r (fm)")
plt.ylabel("Wave Function")
plt.title("Anharmonic Oscillator Ground State")
plt.legend()
plt.grid()
plt.show()
