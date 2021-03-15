'''
Documentação:
Programa de propagação de doenças pelo método de Euler 
'''

# Biblioteca
import numpy as np
import math
import matplotlib.pyplot as plt

# Funções
def f(x, t):
    # Documentação:
    """
    Sistema de equações diferenciais para simular a propagação de doenças:
    
    dS/dt = -beta*S*I
    dI/dt = beta*S*I - ni*I
    dR/dt = ni*I

    Parâmetros: x (array), t (array)
    Retorno: (array)
    """
    beta = 0.0005           # Constante relacionada à probabilidade de transmissão
    ni = 0.1                # Constante relacionada à duração da epidemia
    return np.array([-beta*x[0]*x[1] , beta*x[0]*x[1] - ni*x[1], ni*x[1]])

def Euler(f, xi, ti, tf, N):
    # Documentação:
    """
    Método de Euler para x' = f(x,t), x(0) = xi, com N passos de ti até tf

    Parâmetros: f(x,t), xi (float), ti (float), tf (float), N (int)
    Retorno: x (array), t (array)
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
        print(x[i, 0], x[i, 1], x[i, 2], t[i], x[i, 0] + x[i, 1] + x[i, 2])             # S(t), I(t), R(t), t, S(t) + I(t) + R(t)
        x[i+1] = x[i] + (f(x[i], t[i]))*h
        t[i+1] = t[i] + h

    return x, t

# Parâmetros
ti = 0.
tf = 90.
xi = np.array([1500, 1, 0])         # Condições iniciais [S = indivíduos suscetíveis, I = indivíduos infectados, R = indivíduos imunes]
N = 120

# Método de Euler
x, t = Euler(f, xi, ti, tf, N)

# Gráfico da evolução temporal das curvas S, I e R
plt.plot(t, x[:, 0], 'r-', label = 'S = suscetíveis')
plt.plot(t, x[:, 1], 'g-', label = 'I = infectados')
plt.plot(t, x[:, 2], 'b-', label = 'R = imunes')
plt.xlabel("t")
plt.ylabel("S(t), I(t), R(t)")
plt.title("Evolução temporal das curvas S, I e R")
plt.legend()
plt.grid(True)
plt.show()