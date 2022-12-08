# -*- coding: utf-8 -*-
from __future__ import annotations
import csv
import datetime
import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple


def to_csv(
    output_file: str,
    host_name: str,
    variant: str,
    num_cols: int,
    num_runs: int,
    runtime_mean: float,
    runtime_stddev: float,
) -> None:
    """Write mean and standard deviation of measured runtimes to a CSV file."""
    if not os.path.exists(output_file):
        with open(output_file, "w") as csv_file:
            writer = csv.writer(csv_file, delimiter=",")
            writer.writerow(("date", "host", "variant", "num_cols", "num_runs", "mean", "stddev"))
    with open(output_file, "a") as csv_file:
        writer = csv.writer(csv_file, delimiter=",")
        writer.writerow(
            (
                datetime.date.today().strftime("%Y%m%d"),
                host_name,
                variant,
                num_cols,
                num_runs,
                runtime_mean,
                runtime_stddev,
            )
        )


def print_performance(runtimes: list[float]) -> Tuple[float, float]:
    """Print means and standard deviation of measure runtimes to screen."""
    n = len(runtimes)
    mean = sum(runtimes) / n
    stddev = (sum((runtime - mean) ** 2 for runtime in runtimes) / (n - 1 if n > 1 else n)) ** 0.5
    print(f"Performance: Average runtime over {n} runs: {mean:.3f} \u00B1 {stddev:.3f} ms.")
    return mean, stddev
