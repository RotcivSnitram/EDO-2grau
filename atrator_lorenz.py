'''
Documentação:
Programa que calcula EDO de primeira ordem pelo método de Euler 
'''

# Biblioteca
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Funções
def funcaoderiv(x, y, z):
    # Documentação:
    """
    Derivadas da EDO de primeira ordem

    Parâmetros: x (array), y (array), z (array)
    Retorno: (array)
    """ 
    sigma = 10.
    pho = 28.
    beta = 8/3
    return np.array([sigma*(y - x), pho*x - y - x*z, x*y - beta*z])

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
        #print(r[i, 0], r[i, 1], r[i, 2], t[i])
        r[i+1] = r[i] + (funcaoderiv(r[i, 0], r[i, 1], r[i, 2]))*h
        t[i+1] = t[i] + h

    return r, t

# Parâmetros
ti = 0.
tf = 50.
ri = np.array([0, 1, 0])
N = 5000

# Método de Euler
r, t = Euler(funcaoderiv, ri, ti, tf, N)

# Gráfico y x t
plt.plot(t, r[:, 1], 'r-', label = 'numérico')
plt.xlabel("t")
plt.ylabel("y(t)")
plt.title("Atrator de Lorenz Y x t")
plt.legend()
plt.grid(True)
plt.show()

# Gráfico z x x
plt.plot(r[:, 0], r[:, 2], 'r-', label = 'numérico')
plt.xlabel("x")
plt.ylabel("z")
plt.title("Atrator de Lorenz Z x X")
plt.legend()
plt.grid(True)
plt.show()

# Gráfico 3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(r[:, 0], r[:, 1], r[:, 2], lw=0.5)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("Atrator de Lorenz")
plt.show()