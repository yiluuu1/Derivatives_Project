import math
from scipy.stats import norm

def eu_bs(s, k, t, sigma, r, q):
    """B-S formula"""
    d1 = (math.log((s * math.exp(-q * t)) / (k * math.exp(-r * t))) + 0.5 * (sigma ** 2) * t)/(sigma * math.sqrt(t))
    d2 = d1 - sigma * math.sqrt(t)
    pre = (s * math.exp(-q * t)) * norm.cdf(d1) - (k * math.exp(-r * t)) * norm.cdf(d2)
    return pre


def implied_sigma(s, k, t, r, q, c):
    """dichotomy to solve sigma reversely"""
    sigma_min = 0.00001
    sigma_max = 1000
    sigma_mid = (sigma_min + sigma_max) / 2
    call_min = eu_bs(s, k, t, sigma_min, r, q)
    call_max = eu_bs(s, k, t, sigma_max, r, q)
    call_mid = eu_bs(s, k, t, sigma_mid, r, q)
    diff = c - call_mid
    if c < call_min or c > call_max:
        print('error, the price of option is beyond the limit')
    else:
        while abs(diff) > 1e-6:
            if c > call_mid:
                sigma_min = sigma_mid
            else:
                sigma_max = sigma_mid
            sigma_mid = (sigma_min + sigma_max) / 2
            call_mid = eu_bs(s, k, t, sigma_mid, r,q)
            diff = c - eu_bs(s, k, t, sigma_mid, r,q)
    return sigma_mid

iv=implied_sigma(s=161.94, k=160, t=1, r=math.log(1.01496), q=math.log(1.0058 * 1.0064 * 1.006 * 1.0058), c=18.20)