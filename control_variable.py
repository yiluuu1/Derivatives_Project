from main import binary_tree_Y, delta_y, gamma_y, theta_y
from monte_carlo import y_premium
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

r = math.log(1.01496)
q = math.log(1.0058 * 1.0064 * 1.006 * 1.0058)
weeks = pd.read_csv('AAPL.csv')['Adj Close']
sigma = (weeks.apply(math.log) - weeks.shift(1).apply(math.log))[1:].std() * math.sqrt(52)
s0 = 161.94
n=30

t1=0.75
t2=0.25
k1=23
k2=150
delta=np.array([])
gamma=np.array([])
theta=np.array([])
simu=np.linspace(5, 1000, 200).astype(int)
ys=np.array([])
for n in simu:
    y=y_premium(s0, sigma, r, k1, k2, t1, t2, q, n)
    print(y)
    ys=np.append(ys, y)

plt.plot(simu, ys,'-b', label='Monte Carlo')
plt.axhline(y=7.237064885, xmin=0, xmax=1000, label='Binary Tree', color='red',linestyle="--")
plt.legend(loc="upper left")
plt.show()