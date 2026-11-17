
import numpy as np
import matplotlib.pyplot as plt
import math

m = float(input("Enter reduced mass (MeV/c^2): "))
D = float(input("Enter dissociation energy D (eV): "))
a = float(input("Enter Morse parameter a (1/Angstrom): "))
r0 = float(input("Enter equilibrium distance r0 (Angstrom): "))
N = int(input("Enter number of grid points: "))
rmax = float(input("Enter maximum distance (Angstrom): "))

hbar = 1973
m = m * 1e6

alpha = hbar**2 / (2*m)

r = np.linspace(0.1, rmax, N+2)
dr = r[1] - r[0]
x = r[1:-1]

V = D * (np.exp(-2*a*(x-r0)) - 2*np.exp(-a*(x-r0)))

H = np.zeros((N, N))

for i in range(N):
    H[i, i] = 2*alpha/dr**2 + V[i]

    if i > 0:
        H[i, i-1] = -alpha/dr**2

    if i < N-1:
        H[i, i+1] = -alpha/dr**2

energy, wave = np.linalg.eigh(H)

psi = np.zeros(N+2)
psi[1:-1] = wave[:, 0]

psi = psi / math.sqrt(np.trapezoid(psi**2, r))

print("\nMorse Potential")
print("Ground State Energy =", round(energy[0], 5), "eV")

plt.figure(figsize=(9, 5))
plt.plot(r, np.r_[V[0], V, V[-1]])
plt.xlabel("Distance (Angstrom)")
plt.ylabel("Potential Energy (eV)")
plt.title("Morse Potential")
plt.grid()
plt.show()

plt.figure(figsize=(9, 5))
plt.plot(r, psi)
plt.xlabel("Distance (Angstrom)")
plt.ylabel("Wavefunction")
plt.title("Hydrogen Molecule Ground State")
plt.grid()
plt.show()
