# -*- coding: utf-8 -*-
import csv
import datetime
import os


def to_csv(
    output_file: str,
    host_name: str,
    backend: str,
    nx: int,
    num_runs: int,
    runtime_mean: float,
    runtime_stddev: float,
) -> None:
    if not os.path.exists(output_file):
        with open(output_file, "w") as csv_file:
            writer = csv.writer(csv_file, delimiter=",")
            writer.writerow(("date", "host", "backend", "nx", "num_runs", "mean", "stddev"))
    with open(output_file, "a") as csv_file:
        writer = csv.writer(csv_file, delimiter=",")
        writer.writerow(
            (
                datetime.date.today().strftime("%Y%m%d"),
                host_name,
                backend,
                nx,
                num_runs,
                runtime_mean,
                runtime_stddev,
            )
        )


def print_performance(runtimes: list[float]) -> tuple[float, float]:
    n = len(runtimes)
    mean = sum(runtimes) / n
    stddev = (sum((runtime - mean) ** 2 for runtime in runtimes) / (n - 1 if n > 1 else n)) ** 0.5
    print(f"Performance: Average runtime over {n} runs: {mean:.3f} \u00B1 {stddev:.3f} ms.")
    return mean, stddev
