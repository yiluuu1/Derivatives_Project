import math
import pandas as pd


def binary_s(s0, u, d, node):
    """Calculate s at node(step, up, down)"""
    s = s0 * (d ** node[1]) * (u ** (node[0] - node[1]))
    return s


def forward_tree_O(n, t, r, q, sigma, s0, k, types='E', direct='C'):
    option = {}
    h = t / n
    cn = {}
    u = math.exp((r - q) * h + sigma * math.sqrt(h))
    d = math.exp((r - q) * h - sigma * math.sqrt(h))
    p = (math.exp((r - q) * h) - d) / (u - d)
    for i in range(n + 1):
        if direct == 'P':
            cn[i] = max(k - binary_s(s0, u, d, [n, i]), 0)
        else:
            cn[i] = max(binary_s(s0, u, d, [n, i]) - k, 0)
    option[n] = cn
    for i in reversed(range(n)):
        cx = {}
        for j in range(i + 1):
            cu = option[i + 1][j]
            cd = option[i + 1][j + 1]
            sx = binary_s(s0, u, d, [i, j])
            cr = math.exp(-r * h) * (p * cu + (1 - p) * cd)
            if types == 'A':
                if direct == 'P':
                    cx[j] = max(k - sx, cr)
                else:
                    cx[j] = max(sx - k, cr)
            else:
                cx[j] = cr
        option[i] = cx
    return option


def binary_tree_Y(n1, n2, t1, t2, r, q, sigma, s0, k1, k2, types='E', direct='C'):
    y = {}
    h1 = t1 / n1
    yn = {}
    u1 = math.exp((r - q) * h1 + sigma * math.sqrt(h1))
    d1 = math.exp((r - q) * h1 - sigma * math.sqrt(h1))
    p1 = (math.exp((r - q) * h1) - d1) / (u1 - d1)
    for i in range(n1 + 1):
        s = binary_s(s0, u1, d1, [n1, i])
        c = forward_tree_O(n2, t2, r, q, sigma, s, k2)[0][0]
        if direct == 'P':
            yn[i] = max(k1 - c, 0)
        else:
            yn[i] = max(c - k1, 0)
    y[n1] = yn
    for i in reversed(range(n1)):
        yx = {}
        for j in range(i + 1):
            cu = y[i + 1][j]
            cd = y[i + 1][j + 1]
            cr = math.exp(-r * h1) * (p1 * cu + (1 - p1) * cd)
            s = binary_s(s0, u1, d1, [i, j])
            c = forward_tree_O(n2, t2, r, q, sigma, s, k2)[0][0]
            if types == 'A':
                if direct == 'P':
                    yx[j] = max(k1 - c, cr)
                else:
                    yx[j] = max(c - k1, cr)
            else:
                yx[j] = cr
        y[i] = yx
    return y


def delta_y(q, n1, t1, y, u, d, s, node=None):
    if node is None:
        node = [0, 0]
    h1 = t1 / n1
    delta = math.exp(-q * h1) * (y[node[0] + 1][node[1]] - y[node[0] + 1][node[1] + 1]) / (
            binary_s(s, u, d, [node[0] + 1, node[1]]) - binary_s(s, u, d, [node[0] + 1, node[1] + 1]))
    return delta


def gamma_y(q, n1, t1, y, u, d, s):
    gamma = (delta_y(q, n1, t1, y, u, d, s, [1, 0]) - delta_y(q, n1, t1, y, u, d, s, [1, 1])) / (
            binary_s(s, u, d, [1, 0]) - binary_s(s, u, d, [1, 1]))
    return gamma


def theta_y(n1, t1, q, y, u, d, s):
    eps = u * d * s - s
    h = t1 / n1
    theta = (y[2][1] - eps * delta_y(q, n1, t1, y, u, d, s) - 0.5 * eps ** 2 * gamma_y(q, n1, t1, y, u, d, s) -
             y[0][0]) / (2 * h)
    return theta


if __name__ == '__main__':
    r = math.log(1.01496)
    q = math.log(1.0058 * 1.0064 * 1.006 * 1.0058)
    weeks = pd.read_csv('AAPL.csv')['Adj Close']
    sigma = (weeks.apply(math.log) - weeks.shift(1).apply(math.log))[1:].std() * math.sqrt(52)
    s0 = 161.94
    n1=50
    n2=50
    for t1 in [0.25, 0.5, 0.75, 1]:
        for t2 in [0.25, 0.5, 0.75, 1]:
            if t1 + t2 == 0.5:  # March 2022
                for k1, k2 in zip([7, 11, 17], [170, 160, 150]):
                    y = binary_tree_Y(n1, n2, t1, t2, r, q, sigma, s0, k1, k2, types='A')
                    print(y[0][0])

                    h1 = t1 / n1
                    u = math.exp((r - q) * h1 + sigma * math.sqrt(h1))
                    d = math.exp((r - q) * h1 - sigma * math.sqrt(h1))
                    delta=delta_y(q, n1, t1, y, u, d, s0)
                    gamma=gamma_y(q, n1, t1, y, u, d, s0)
                    theta=theta_y(n1, t1, q, y, u, d, s0)

            elif t1 + t2 == 0.75:  # June 2022
                for k1, k2 in zip([11, 15, 21], [170, 160, 150]):
                    y = binary_tree_Y(n1, n2, t1, t2, r, q, sigma, s0, k1, k2, types='A')
                    print(y[0][0])
            elif t1 + t2 == 1:  # Sep 2022
                for k1, k2 in zip([14, 18, 23], [170, 160, 150]):
                    y = binary_tree_Y(n1, n2, t1, t2, r, q, sigma, s0, k1, k2, types='A')
                    print(y[0][0])

                    # h1 = t1 / n1
                    # u = math.exp((r - q) * h1 + sigma * math.sqrt(h1))
                    # d = math.exp((r - q) * h1 - sigma * math.sqrt(h1))
                    # delta=delta_y(q, n1, t1, y, u, d, s0)
                    # gamma=gamma_y(q, n1, t1, y, u, d, s0)
                    # theta=theta_y(n1, t1, q, y, u, d, s0)

            elif t1 + t2 == 1.25:  # Dec 2022
                for k1, k2 in zip([18, 22, 28], [170, 160, 150]):
                    y = binary_tree_Y(n1, n2, t1, t2, r, q, sigma, s0, k1, k2, types='A')
                    print(y[0][0])
            elif t1 + t2 == 1.5:  # March 2023
                for k1, k2 in zip([19, 23, 29], [170, 160, 150]):
                    y = binary_tree_Y(n1, n2, t1, t2, r, q, sigma, s0, k1, k2, types='A')
                    print(y[0][0])
            elif t1 + t2 == 1.75:  # June 2023
                for k1, k2 in zip([21, 26, 31], [170, 160, 150]):
                    y = binary_tree_Y(n1, n2, t1, t2, r, q, sigma, s0, k1, k2, types='A')
                    print(y[0][0])
            elif t1 + t2 == 2:  # Sep 2023
                for k1, k2 in zip([24, 28, 33], [170, 160, 150]):
                    y = binary_tree_Y(n1, n2, t1, t2, r, q, sigma, s0, k1, k2, types='A')
                    print(y[0][0])

                    # h1 = t1 / n1
                    # u = math.exp((r - q) * h1 + sigma * math.sqrt(h1))
                    # d = math.exp((r - q) * h1 - sigma * math.sqrt(h1))
                    # delta=delta_y(q, n1, t1, y, u, d, s0)
                    # gamma=gamma_y(q, n1, t1, y, u, d, s0)
                    # theta=theta_y(n1, t1, q, y, u, d, s0)
