import matplotlib.pyplot as plt
import numpy as np
import csv


def get_coordinates(benchmark_name):
    bit_sizes = []
    time_averages = []

    with open(benchmark_name, 'r') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            row_b = [int(row[0])]
            row_f = [float(i) for i in row[1:]]
            row = row_b + row_f
            bit_sizes.append(row[0])
            bit_sizes_repeated = bit_sizes_repeated + ([row[0]] * len(row[1:]))
            times = times + row[1:]
            time_averages.append(sum(row[1:])/len(row[1:]))

    return (bit_sizes, time_averages) 


def plot_benchmark(benchmark_names, benchmark_titles, plot_title, filename):
    bit_sizes = []
    times_averages = []

    for name in benchmark_names:
        (xs, ys) = get_coordinates(name)
        bit_sizes.append(xs)
        times_averages.append(ys)

    for xs, ys, title in zip(bit_sizes, times_averages, benchmark_titles):
        plt.plot(xs, ys, label=title)

    plt.xscale('log', base=2)

    plt.xlabel('Bit Sizes')
    plt.ylabel('Average Time (s)')

    plt.title(plot_title)
    plt.legend()
    plt.save(filename)


def main():
    benchmark_names = ["benchmark_results_td_1.csv", "benchmark_results_td_8.csv", "benchmark_results_tdc.csv"]
    benchmark_titles = ["1 thread", "2 thread", "8 thread", "GPU"]
    plot_title = "Trial Division (C++)"
    filename = "td_cpp_graph.png"
    plot_benchmark(benchmark_names, benchmark_titles, plot_title, filename)


if __name__ == "__main__":
    main()