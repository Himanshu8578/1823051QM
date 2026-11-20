
import numpy as np
import matplotlib.pyplot as plt
import math

l = int(input("Enter l: "))
m = int(input("Enter m: "))

if l < 0 or abs(m) > l:
    print("Invalid quantum numbers")
    exit()

theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)

TH, PH = np.meshgrid(theta, phi)

def P(l, m, x):
    p = np.ones_like(x)

    for i in range(1, m+1):
        p *= -(2*i-1)*np.sqrt(1-x**2)

    if l == m:
        return p

    q = x*(2*m+1)*p

    if l == m+1:
        return q

    for k in range(m+2, l+1):
        r = ((2*k-1)*x*q-(k+m-1)*p)/(k-m)
        p = q
        q = r

    return q

p = P(l, abs(m), np.cos(TH))

C = math.sqrt((2*l+1)/(4*math.pi) *
              math.factorial(l-abs(m)) /
              math.factorial(l+abs(m)))

Y = C*p*np.exp(1j*m*PH)
D = np.abs(Y)**2

R = np.sqrt(D/D.max())

X = R*np.sin(TH)*np.cos(PH)
Y1 = R*np.sin(TH)*np.sin(PH)
Z = R*np.cos(TH)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, Y1, Z, cmap='jet')

ax.set_title("Spherical Harmonics")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()

plt.contourf(PH, TH, D, 50, cmap='jet')
plt.colorbar()
plt.xlabel("Phi")
plt.ylabel("Theta")
plt.title("Probability Density")
plt.show()
