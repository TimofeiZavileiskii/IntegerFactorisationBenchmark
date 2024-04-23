import math
from bitarray import bitarray


def sieve_of_eratosthenes(B):
    bound = int(math.floor(math.sqrt(B)))
    sieve = bitarray('0') * B
    primes = []

    for i in range(2, bound+1):
        if not sieve[i]:
            primes.append(i)
            for ii in range(i, B, i):
                sieve[ii] = True

    for i in range(bound+1, B):
        if not sieve[i]:
            primes.append(i)

    return primes


def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)


def modinv(a, m):
    g, x, y = egcd(a, m)
    print(g)
    if g != 1:
        return None
    else:
        return x % m


def is_quadratic_residue(p, n):
    if n == 2:
        return True

    val = pow(n, (p-1)//2, p)
    return val == 1


def sieve_of_eratosthenes_check_residue(B, n):
    bound = int(math.floor(math.sqrt(B)))
    sieve = bitarray('0') * B
    primes = []

    for i in range(2, bound+1):
        if not sieve[i]:
            if is_quadratic_residue(i, n):
                primes.append(i)
            for ii in range(i, B, i):
                sieve[ii] = True

    for i in range(bound+1, B):
        if not sieve[i] and is_quadratic_residue(i, n):
            primes.append(i)

    return primes