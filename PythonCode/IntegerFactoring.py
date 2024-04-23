from RhoAlgorithm import PollardRhoFactorisation, TrialDivision
import csv
import time
from QuadraticSieve import quadratic_sieve, sieve_of_eratosthenes, get_roots
from PollardsP1 import pollardsP1
from ECM import ECM
from sympy import Symbol, factor, factor_list


def run_benchmark():
    with open("../rsa_numbers.csv", "r") as f:
        reader = csv.reader(f)
        time_spent = []
        f_w = open("rsa_benchmark_python.csv", "w")
        f_w.close()
        for row in reader:
            f_w = open("rsa_benchmark_python.csv", "a")
            row = [int(i) for i in row]
            num_bits = row[0]
            f_w.write(f"{num_bits}")
            f_w.close()
            print(f"---------Factorising numbers of {num_bits} bits----------")
            for rsa_num in row[1:]:
                f_w = open("rsa_benchmark_python.csv", "a")
                f_w.write(",")
                f_w.close()
                print(f"Factorise {rsa_num}")
                start_time = time.time()
                factor = TrialDivision(rsa_num)
                end_time = time.time()
                time_for_factorisation = end_time - start_time
                time_spent.append(time_for_factorisation)
                other_factor = rsa_num // factor
                print(f"Number factorised {rsa_num} = {factor} * {other_factor}")
                print(f"Time spent {time_for_factorisation}")
                f_w = open("rsa_benchmark_python.csv", "a")
                f_w.write(f"{time_for_factorisation}")
                f_w.close()
            f_w = open("rsa_benchmark_python.csv", "a")
            f_w.write("\n")
            f_w.close()
            print(f"Average time for {num_bits} is {sum(time_spent)/len(time_spent)}")


def manual_input():
    n = int(input("Enter the number to factorise:\n"))
    p = PollardRhoFactorisation(n)

    if n % p != 0:
        print("Something went wrong")
        print(f"Factor {p}")
        raise Exception()

    q = n // p
    print(f"n = {p}*{q}")


def main():
    run_benchmark()


if __name__ == "__main__":
    main()
