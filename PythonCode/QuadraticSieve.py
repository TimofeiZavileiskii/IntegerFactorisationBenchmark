import math
from dataclasses import dataclass
from bitarray import bitarray
from math import gcd
from typing import List
from sympy import isprime
import galois
from Utils import sieve_of_eratosthenes, is_quadratic_residue, sieve_of_eratosthenes_check_residue


@dataclass
class Relation:
    x: int
    y: int
    powers: List[int]

    def relation(self, x, y, powers):
        self.x = x
        self.y = y
        self.powers = powers

    def multiply(self, other, n):
        self.x = (self.x * other.x)
        self.y = (self.y * other.y)

        for n, i in enumerate(other.powers):
            self.powers[n] += i

    def __str__(self):
        return f"Relation x:{self.x}, y:{self.y}, powers:{self.powers}"


# find roots in relation (start + x)^2 - n = 0 (mod p)
def get_roots(n, start, prime):
    n %= prime

    # Find a quadratic non-residue
    non_residue = 0
    for i in range(2, prime):
        if not is_quadratic_residue(prime, i):
            non_residue = i
            break
    #print(f"The non-residue is {non_residue}")
    f = prime - 1
    #print(f"F is {f}")
    e = 0
    while f % 2 == 0:
        f //= 2
        e += 1

    #print(f"The decomp is p = {f}*2^{e}")
    R = pow(n, f, prime)
    N = pow(non_residue, f, prime)
    j = 0

    for i in range(1, e):
        product = pow(R * pow(N, j, prime), pow(2, e-i-1), prime)
        if product == prime - 1:
            j += pow(2, i)
    root = pow(n, (f+1)//2, prime) * pow(N, j//2, prime) % prime

    #print(f"The root is {root}")

    roots = [(root - start) % prime, (prime - root - start) % prime]
    return roots


def poly(start, x, n):
    return (start + x)**2 - n


def polynomial_sieve(n, B, primes):
    print(f"Prime size is {len(primes)}")
    start = math.ceil(math.sqrt(n))
    get_numbers = 6*B

    test_numbers = [poly(start, i, n) for i in range(get_numbers)]
    powers = [[0 for ii in primes] for i in range(get_numbers)]

    #print(f"The sieved numbers: {test_numbers}")

    for prime_num, prime in enumerate(primes):
        roots = get_roots(n, start, prime)
        if roots[1] == roots[0]:
            roots.pop(1)

        for root in roots:
            for i in range(root, get_numbers, prime):
                power = 0
                while test_numbers[i] % prime == 0:
                    test_numbers[i] //= prime
                    power += 1
                if power == 0:
                    print("Poly is not factored =(")
                    print(f"Prime: {prime}, root: {root}, n:{n}, number:{test_numbers[i]}")
                powers[i][prime_num] = power
    #print(f"End power nums: {test_numbers}")
    output_matrix = []

    max_size = len(primes) + 6

    for i, test_number in enumerate(test_numbers):
        if test_numbers[i] == 1:
            output_matrix.append(Relation((start+i), (start+i)**2 - n, powers[i]))
            if len(output_matrix) > max_size:
                break
    return output_matrix


def verify_relation(n, primes, relation):
    factor_correct = True

    y = relation.y
    for i, prime in enumerate(primes):
        factor_power = 0
        while y % prime == 0:
            factor_power = (factor_power + 1) % 2
            y //= prime

        if relation.powers[i] != factor_power:
            factor_correct = False
            print(f"Wrong factorisation for prime{prime}")

    return (relation.x ** 2 - n == relation.y) and factor_correct


def construct_square(n, primes, power_matrix):
    GF = galois.GF(2)
    matrix = [[i & 1 for i in relation.powers] for relation in power_matrix]

    M_gf = GF(matrix)

    null_space = M_gf.left_null_space()[1]

    square_relation = Relation(1, 1, [0 for i in primes])

    for i, entry in enumerate(null_space):
        if entry == 1:
            #print(power_matrix[i])
            square_relation.multiply(power_matrix[i], n)

    y = 1
    for i, p in enumerate(square_relation.powers):
        y *= pow(primes[i], p//2)
    x = square_relation.x
    print(f"Nums: {x}, {y}")
    return x, y


def quadratic_sieve(n):
    B = int(math.exp(0.5 + math.sqrt(math.log(n) * math.log(math.log(n)))))
    print(f"The bound is {B}")
    primes = sieve_of_eratosthenes_check_residue(B, n)
    print(f"Primes obtained")

    power_matrix = polynomial_sieve(n, B, primes)
    print(f"Numbers sieved")

    sq1, sq2 = construct_square(n, primes, power_matrix)
    x = get_roots(sq1, 0, n)[0]
    y = get_roots(sq2, 0, n)[1]
    print(f"x:{x} x^2:{x ** 2} x^2 % n: {x ** 2 % n == y ** 2 % n}")

    print(f"X: {x} Y: {y}")
    print(f"Factor: {(x + y) % n}")
    factor = gcd((sq1 - sq2) % n, n)

    return factor