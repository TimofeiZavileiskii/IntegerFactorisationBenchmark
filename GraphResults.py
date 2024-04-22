import matplotlib.pyplot as plt
import numpy as np
import csv


def get_coordinates(benchmark_name):
    bit_sizes = []
    time_averages = []

    with open(benchmark_name, 'r') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            if len(row) < 3:
                continue

            bit_size = [int(row[0])]
            average_time = [float(row[-2])]
            bit_sizes.append(bit_size)
            time_averages.append(average_time)

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

    plt.yscale('log', base=2)

    plt.xlabel('Bit Sizes')
    plt.ylabel('Average Time (s)')

    plt.title(plot_title)
    plt.legend()
    plt.savefig(filename)


def main():
    benchmark_names = ["benchmark_results_td_1.csv", "benchmark_results_td_2.csv", "benchmark_results_td_8.csv", "benchmark_results_tdc_2048.csv"]
    benchmark_titles = ["1 thread", "2 thread", "8 thread", "GPU"]
    plot_title = "Trial Division (C++)"
    filename = "td_cpp_graph.png"
    plot_benchmark(benchmark_names, benchmark_titles, plot_title, filename)


if __name__ == "__main__":
    main()