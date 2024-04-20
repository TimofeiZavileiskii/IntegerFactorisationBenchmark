import sage.all as fn
import scipy.integrate as integrate
import scipy.optimize as optimize
from math import log, log2, exp, sqrt


def mu(a, b):
   # upper = min(a-1, 50)
    s, _ = integrate.quad(lambda x: fn.dickman_rho(x), a=(a-b), b=(a-1), limit=500)
  #  print(f"Integral is {s} for {a-b} {a-1}")
    return s/(a-b)


def curve_performance(p, k, b1, b2):
    a = log(p)/log(b1)
    b = log(b2)/log(b1)
    t = (b1 + (b2-b1)/k)

   # print(f"A: {a} B: {b} for p: {p} k: {k}")
    return mu(a, b)/t


def L(p):
    return exp((sqrt(2)/2)*sqrt(log(p)*log(log(p))))


def constraint(vars):
    x, y = vars
    return y - x


def optimised_params(x_range, y_range, p, k, factor):
    B1 = x_range[0]
    B2 = y_range[0]

    best_coords = [B1, B2]
    max_val = 0

    while(B1 < x_range[1]):
        B2 = B1
        while(B2 < y_range[1]):
            val = curve_performance(p, k, B1, B2)
            if val > max_val:
                max_val = val
                best_coords = [B1, B2]
                #print(f"Best coords B1: {B1} B2: {B2} at val {val}")
            B2 *= factor
        B1 *= factor
    return best_coords
            

def optimise_curve_performance(p, k):
    factor = 2**p
   # _ = curve_performance(p, k, L(p), L(p)*100)
#    return None
  #  to_optimize = lambda vars: -curve_performance(factor, k, vars[0], vars[1]) + 5
    inititial_guess = [L(factor), L(factor)*100]
    x_range = (2, factor/3)
    y_range = (3, factor/2)
    multiple = 1.1
    #constraint_defs = {"type": "ineq", "fun": constraint}
    #result = optimize.shgo(to_optimize, bounds=bounds, constraints=constraint_defs, n=1000, iters=5)
    result = optimised_params(x_range, y_range, factor, k, multiple)

    x_range = (result[0]/2, result[0]*2)
    y_range = (result[1]/2, result[1]*2)

    multiple = 1.05
    result = optimised_params(x_range, y_range, factor, k, multiple)

    print(f"Inital guess: {inititial_guess[0]} {inititial_guess[1]}")
    print(f"Optimised value for p {p}: B1={result[0]} B2={result[1]}")
    return result


def write_to_file(string, filename, mode):
    f = open(filename, mode)
    f.write(f"{string}\n")
    f.close()


def main():
    k = 27
    start = 10
    end = 315

    filename = "bounds.txt"

    write_to_file(f"{start}", filename, "w")

    b_values = []

    for i in range(start, end):
        b_values.append(optimise_curve_performance(i, k))
        write_to_file(f"{b_values[-1][0]},{b_values[-1][1]}", filename, "a")

    write_to_file(b_values, start, "bounds.txt")


if __name__ == "__main__":
    main()
