import math
from Utils import sieve_of_eratosthenes
from dataclasses import dataclass
from random import randint

@dataclass
class Point:
    x: any
    y: any
    curve: any

    def __add__(self, other):
        return self.curve.add_points(self, other)

    def __sub__(self, other):
        return self.curve.sub_points(self, other)


# https://stackoverflow.com/questions/4798654/modular-multiplicative-inverse-function-in-python

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)


def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise FactorisedException(g % m, 'modular inverse does not exist')
    else:
        return x % m


class FactorisedException(Exception):
    def __init__(self, factor, message):
        self.factor = factor
        super().__init__(message)


class Curve:
    def __init__(self, a, b, n):
        self.a = a
        self.b = b
        self.n = n

    def get_random_point(self):
        return Point(1, 1, self)

    def add_points(self, p1, p2):
        if p1.x == "inf":
            return p2
        if p2.x == "inf":
            return p1
        if p1.x == p2.x and p1.y == self.n - p2.y:
            p1.x = "inf"
            return p1

        s = 1
        if p1.x == p2.x:
            y2 = modinv((2 * p1.y) % self.n, self.n)
            s = ((3 * (p1.x ** 2) + self.a) * y2) % self.n
        else:
            x_inv = modinv((p2.x - p1.x) % self.n, self.n)
            s = (((p2.y - p1.y) % self.n) * x_inv) % self.n

        x3 = (s**2 - p1.x - p2.x) % self.n
        y3 = (s * (p1.x - x3) - p1.y) % self.n
        return Point(x3, y3, self)

    def sub_points(self, p1, p2):
        p2.y = self.n - p2.y
        return self.add_points(p1, p2)

    def mult_point(self, p1, mult):
        if mult == 1:
            return p1

        product = self.mult_point(p1, mult // 2)
        square = self.add_points(product, product)
        if mult & 1 == 1:
            return self.add_points(p1, square)
        else:
            return square


def ECM(n):
    B = int(math.sqrt(n))
    primes = sieve_of_eratosthenes(B)
    factor = 1

    while factor == 1:
        a = randint(1, 10000) % n
        x1 = randint(1, 10000) % n
        y1 = randint(1, 10000) % n
        b = (y1**2 - x1**3 - a*x1) % n

        g = math.gcd((4 * pow(a, 3, n) + 27 * (b**2)), n)
        if g == n:
            continue
        if g > 1:
            return g

        curve = Curve(a, b, n)
        point = Point(x1, y1, curve)
        try:
            for i, prime in enumerate(primes):
                e = int(math.floor(math.log(B)/math.log(prime)))
                f = pow(prime, e)
                point = curve.mult_point(point, f)

        except FactorisedException as e:
            factor = e.factor
            break
        print(curve)

    return factor
