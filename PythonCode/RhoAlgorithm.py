#Implementation of the Pollard's Rho algorithm
import random
from math import gcd

def TrialDivision(n):
    trial = 2
    while trial * trial < n:
        if n % trial == 0:
            return n // trial
        trial += 1

    return n

def polynomial(n, p):
    return (n*n + 1) % p

def PollardRhoFactorisation(n):
    factored = False
    used_seeds = set()
    factor = 1
    while not factored:
        seed = random.randint(1, n - 1)
        while seed in used_seeds:
            seed = random.randint(1, n - 1)

        u = seed
        v = seed
        u = polynomial(u, n)
        v = polynomial(polynomial(v, n), n)
        while gcd((u-v) % n, n) == 1 and u-v != 0:
            u = polynomial(u, n)
            v = polynomial(polynomial(v, n), n)

        if gcd((u-v) % n, n) > 1:
            factored = True
            factor = gcd((u-v) % n, n)

    return factor
