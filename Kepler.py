'''
Documentação:
Programa que calcula a translação da Terra pelo método de Euler 
'''

# Biblioteca
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from rk4 import RungeKutta4

# Funções
def f(x, t):
    # Documentação:
    """
    Derivadas da EDO

    Parâmetros: x (array), y (array)
    Retorno: (array)
    """
    G = 6.67*10**(-11)
    Ms = 1.98*10**(30)
    return np.array([x[1], -G*Ms*(x[0]/(((x[0]**2 + x[2]**2)**(1/2))**3)), x[3], -G*Ms*(x[2]/(((x[0]**2 + x[2]**2)**(1/2))**3))])

def Euler(f, xi, ti, tf, N):
    # Documentação:
    """
    Método de Euler para r' = f(x,y), r(0) = ri, com N passos de ti até tf

    Parâmetros: funcaoderiv, xi (float), ti (float), tf (float), N (int)
    Retorno: r (array), t (array)
    """
    t = np.zeros(N)

    if isinstance(xi, (float, int)):
        x = np.zeros(N)
    else:
        neq = len(xi)
        x = np.zeros((N, neq))

    x[0] = xi
    t[0] = ti
    h = (tf - ti)/float(N)

    for i in range(N-1):
        #print(x[i, 0], x[i, 1], x[i, 2], x[i, 3], t[i])
        x[i+1] = x[i] + (f(x[i], t[i]))*h
        t[i+1] = t[i] + h

    return x, t

# Parâmetros
ti = 0.
tf = 31536000.
xi = np.array([1.496*10**(11), 0, 0, 2.97*10**(4)])
N = 365

# Método de Euler
x, t = Euler(f, xi, ti, tf, N)

# Método de Runge Kutta 4° ordem
x1, t1 = RungeKutta4(f, xi, ti, tf, N)

# Gráfico Y x t
plt.plot(t, x[:, 2], 'r-', label = 'numérico')
plt.xlabel("t")
plt.ylabel("y(t)")
plt.title("Gráfico y(t) x t (Método de Euler)")
plt.legend()
plt.grid(True)
plt.show()

# Gráfico Y x X
plt.plot(x[:, 0], x[:, 2], 'r-', label = 'h = 1 dia')
plt.xlabel("x")
plt.ylabel("y")
plt.title("Trajetória da Terra ao redor do Sol")
plt.legend()
plt.axes().set_aspect('equal')
plt.grid(True)
plt.show()

# Gráfico 3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(t, x[:, 0], x[:, 2], lw=0.5)
ax.set_xlabel("t")
ax.set_ylabel("X")
ax.set_zlabel("Y")
ax.set_title("Gráfico XY x t (Método de Euler)")
plt.show()

# Gráfico Y x t
plt.plot(t1, x1[:, 2], 'r-', label = 'numérico')
plt.xlabel("t")
plt.ylabel("y(t)")
plt.title("Gráfico y(t) x t (Método de Runge Kutta de 4° Ordem)")
plt.legend()
plt.grid(True)
plt.show()

# Gráfico Y x X
plt.plot(x1[:, 0], x1[:, 2], 'r-', label = 'h = 1 dia')
plt.xlabel("x")
plt.ylabel("y")
plt.title("Trajetória da Terra ao redor do Sol (Runge-Kutta)")
plt.legend()
plt.axes().set_aspect('equal')
plt.grid(True)
plt.show()

# Gráfico 3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(t1, x1[:, 0], x1[:, 2], lw=0.5)
ax.set_xlabel("t")
ax.set_ylabel("X")
ax.set_zlabel("Y")
ax.set_title("Gráfico XY x t (Método de Runge Kutta de 4° Ordem)")
plt.show()