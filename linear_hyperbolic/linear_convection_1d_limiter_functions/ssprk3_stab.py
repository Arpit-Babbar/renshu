import numpy as np
import matplotlib.pyplot as plt


N = 500

zr = np.linspace(-5.0,1.0,N)
zi = np.linspace(-4.0,4.0,N)

zr,zi = np.meshgrid(zr,zi)
z = zr + 1j * zi

G = 1.0 + z + z**2/2.0 + z**3/6
sigma = 0.6282
x = np.linspace(-4.0*np.pi,4*np.pi,100)
z = sigma*(-1.5 + 2.0 * np.cos(x) - 0.5 * np.cos(2.0*x) + (-2.0*np.sin(x) + 0.5*np.sin(2.0*x))*1j)
X = [x.real for x in z]
Y = [x.imag for x in z]

plt.figure(figsize=(6,5))
cs = plt.contourf(zr,zi,np.abs(G),levels=30,cmap='jet')
plt.colorbar(cs)
plt.contour(zr,zi,np.abs(G),levels=[0.9999,1,1.0001],linewidths=2)
plt.axhline(0, color='white')
plt.axvline(0, color='white')
plt.text(-4, 0, 'Unstable region', verticalalignment='center',
         color='white', rotation=90,bbox=dict(facecolor='red'))
plt.title('RK4: Magnitude of amplification factor')
plt.xlabel('Real z')
plt.ylabel('Imaginary z')
plt.scatter(X,Y, color='red')
plt.show()
