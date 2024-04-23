import math
from bitarray import bitarray
from Utils import sieve_of_eratosthenes


def pollardsP1(n):
    a = 2
    B = int(math.pow(n, 1/3))
    primes = sieve_of_eratosthenes(B)

    for i, p in enumerate(primes):
        e = int(math.floor(math.log2(B)/math.log2(p)))
        f = pow(p, e)
        a = pow(a, f, n)
        if i % 10 == 0:
            q = math.gcd(a-1, n)
            if q > 1:
                return q
    return 1