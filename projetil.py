'''
Documentação:
Programa que calcula a trajetória de um projétil pelo método de Euler 
'''

# Biblioteca
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Funções
def f(r, t):
    # Documentação:
    """
    Derivadas da EDO de primeira ordem

    Parâmetros: x (array), y (array), z (array)
    Retorno: (array)
    """ 
    B2m = 4*10**(-5)
    g = 9.78
    return np.array([r[1], -(B2m)*(((r[1])**2 + (r[3])**2)**(1/2))*r[1], r[3], -(B2m)*(((r[1])**2 + (r[3])**2)**(1/2))*r[3] - g])

def Euler(f, ri, ti, tf, N):
    # Documentação:
    """
    Método de Euler para r' = f(x,y,z), r(0) = ri, com N passos de ti até tf

    Parâmetros: funcaoderiv, ri (float), ti (float), tf (float), N (int)
    Retorno: r (array), t (array)
    """
    t = np.zeros(N)

    if isinstance(ri, (float, int)):
        r = np.zeros(N)
    else:
        neq = len(ri)
        r = np.zeros((N, neq))

    r[0] = ri
    t[0] = ti
    h = (tf - ti)/float(N)

    for i in range(N-1):
        print(r[i, 0], r[i, 1], r[i, 2], r[i, 3], t[i])
        r[i+1] = r[i] + (f(r[i], t[i]))*h
        t[i+1] = t[i] + h

    return r, t

# Parâmetros
ti = 0.
tf = 95.
v = 700.
theta = 45
ri = np.array([0, v*np.cos(theta), 0, v*np.sin(theta)])
N = 50

# Método de Euler
r, t = Euler(f, ri, ti, tf, N)

# Gráfico X x t
plt.plot(t, r[:, 0], 'r-', label = 'numérico')
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("Gráfico X x t")
plt.legend()
plt.grid(True)
plt.show()

# Gráfico Y x t
plt.plot(t, r[:, 2], 'r-', label = 'numérico')
plt.xlabel("t")
plt.ylabel("y(t)")
plt.title("Gráfico Y x t")
plt.legend()
plt.grid(True)
plt.show()

# Gráfico Vx x t
plt.plot(t, r[:, 1], 'r-', label = 'numérico')
plt.xlabel("t")
plt.ylabel("Vx(t)")
plt.title("Gráfico Vx x t")
plt.legend()
plt.grid(True)
plt.show()

# Gráfico Vy x t
plt.plot(t, r[:, 3], 'r-', label = 'numérico')
plt.xlabel("t")
plt.ylabel("Vy(t)")
plt.title("Gráfico Vy x t")
plt.legend()
plt.grid(True)
plt.show()

# Gráfico 3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(t, r[:, 0], r[:, 2], lw=0.5)
ax.set_xlabel("t")
ax.set_ylabel("X")
ax.set_zlabel("Y")
ax.set_title("Gráfico XY x t")
plt.show()