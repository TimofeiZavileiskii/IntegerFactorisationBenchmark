import matplotlib.pyplot as plt
import numpy as np
import csv


def get_coordinates(benchmark_name, compute_average):
    bit_sizes = []
    time_averages = []

    with open(benchmark_name, 'r') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            if (len(row) < 3 and (not compute_average)) or (len(row) < 1):
                continue

            bit_size = int(row[0])
            if not compute_average:
                average_time = float(row[-2])
            else:
                average_time = float(row[-1])

            bit_sizes.append(bit_size)
            time_averages.append(average_time)
    return (bit_sizes, time_averages) 


def plot_benchmark(benchmark_names, benchmark_titles, compute_averages, plot_title, filename):
    bit_sizes = []
    times_averages = []

    for name, do_compute_average in zip(benchmark_names, compute_averages):
        (xs, ys) = get_coordinates(name, do_compute_average)
        bit_sizes.append(xs)
        times_averages.append(ys)

    for xs, ys, title in zip(bit_sizes, times_averages, benchmark_titles):
        plt.plot(xs, ys, label=title)

    plt.yscale('log', base=2)

    plt.xlabel('Bit Sizes')
    plt.ylabel('Average Time (s)')

    plt.title(plot_title)
    plt.legend()
    plt.savefig(filename)


def main():
    compute_average = [True, False, False]
    benchmark_names = ["rsa_benchmark_mathematica.csv", "benchmark_results_pari_1.csv", "rsa_benchmark_gmp_ecm.csv"]
    benchmark_titles = ["Wolfram Mathematica", "PARI", "GMP-ECM"]
    plot_title = "Mathematical Software"
    filename = "math_sf_graph.png"
    plot_benchmark(benchmark_names, benchmark_titles, compute_average, plot_title, filename)


if __name__ == "__main__":
    main()