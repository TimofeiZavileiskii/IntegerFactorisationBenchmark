import matplotlib.pyplot as plt
import numpy as np
import csv


def get_coordinates(benchmark_name, compute_average):
    bit_sizes = []
    time_averages = []

    with open(benchmark_name, 'r') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            if len(row) < 3:
                continue
            bit_size = int(row[0])
            if not compute_average:
                average_time = float(row[-2])
            else:
                s = 0
                for t in row[1:]:
                    if(t):
                        s += float(t)
                average_time = s/(len(row)-1)

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

    plt.xlabel('Bit Size')
    plt.ylabel('Average Time (s)')

    plt.title(plot_title)
    plt.legend()
    plt.savefig(filename)


def main():
    compute_average = [False, False, False]
    benchmark_names = ["benchmark_results_ecmm2_8.csv", "benchmark_results_emc_512.csv", "benchmark_results_ecmc_512_mont.csv"]
    benchmark_titles = ["8 threads stage 2", "GPU", "GPU Montgomery"]
    plot_title = "ECM GPU"
    filename = "ecm_gpu.png"
    plot_benchmark(benchmark_names, benchmark_titles, compute_average, plot_title, filename)


if __name__ == "__main__":
    main()