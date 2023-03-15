from numpy.random import randn
import math
import pandas as pd

def cal_S_T(s0, sigma, r, t, q):
    """Simulation ST follows Geometric Brownean Motion """
    return s0 * math.exp((r - q - 0.5 * sigma ** 2) * t + sigma * math.sqrt(t) * randn())

def monte_carlo(s0, sigma, r, k1, k2, t1, t2, q, simulation):
    ysum=0
    for i in range(simulation):
        csum=0
        s1 = cal_S_T(s0, sigma, r, t1, q)  # Simulate S at t1 from t0
        for j in range(simulation):
            s2=cal_S_T(s1, sigma, r, t2, q)  # Simulate S at t2 from t1
            csum+=max(s2-k2, 0)
        c=math.exp(-r *t2)*csum/simulation
        ysum+=max(c-k1, 0)
    y=math.exp(-r *t1)*ysum/simulation
    return y

if __name__ == '__main__':
    r = math.log(1.01496)
    q = math.log(1.0058 * 1.0064 * 1.006 * 1.0058)
    weeks = pd.read_csv('AAPL.csv')['Adj Close']
    sigma = (weeks.apply(math.log) - weeks.shift(1).apply(math.log))[1:].std() * math.sqrt(52)
    s0 = 161.94
    simulation = 1000

    for t1 in [0.25, 0.5, 0.75, 1]:
        for t2 in [0.25, 0.5, 0.75, 1]:
            if t1 + t2 == 0.5:  # March 2022
                for k1, k2 in zip([7, 11, 17], [170, 160, 150]):
                    y = monte_carlo(s0, sigma, r, k1, k2, t1, t2, q, simulation)
                    print(y)

            if t1 + t2 == 0.75:  # June 2022
                for k1, k2 in zip([11, 15, 21], [170, 160, 150]):
                    y = monte_carlo(s0, sigma, r, k1, k2, t1, t2, q, simulation)
                    print(y)
            if t1 + t2 == 1:  # Sep 2022
                for k1, k2 in zip([14, 18, 23], [170, 160, 150]):
                    y = monte_carlo(s0, sigma, r, k1, k2, t1, t2, q, simulation)
                    print(y)

            if t1 + t2 == 1.25:  # Dec 2022
                for k1, k2 in zip([18, 22, 28], [170, 160, 150]):
                    y = monte_carlo(s0, sigma, r, k1, k2, t1, t2, q, simulation)
                    print(y)
            if t1 + t2 == 1.5:  # March 2023
                for k1, k2 in zip([19, 23, 29], [170, 160, 150]):
                    y = monte_carlo(s0, sigma, r, k1, k2, t1, t2, q, simulation)
                    print(y)
            if t1 + t2 == 1.75:  # June 2023
                for k1, k2 in zip([21, 26, 31], [170, 160, 150]):
                    y = monte_carlo(s0, sigma, r, k1, k2, t1, t2, q, simulation)
                    print(y)
            if t1 + t2 == 2:  # Sep 2023
                for k1, k2 in zip([24, 28, 33], [170, 160, 150]):
                    y = monte_carlo(s0, sigma, r, k1, k2, t1, t2, q, simulation)
                    print(y)
