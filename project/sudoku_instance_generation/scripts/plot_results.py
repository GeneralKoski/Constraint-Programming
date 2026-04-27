#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def plot_time_vs_clues(rows, out_path: Path):
    fig, ax = plt.subplots(figsize=(8, 5))
    by_strategy = defaultdict(list)
    for row in rows:
        by_strategy[row["strategy"]].append((int(row["clues"]), float(row["total_seconds"])))
    markers = {"random": "o", "symmetry": "s", "density": "^"}
    for strategy, points in by_strategy.items():
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        ax.scatter(xs, ys, label=strategy, marker=markers.get(strategy, "o"), alpha=0.7, s=40)
    ax.set_xlabel("Indizi rimanenti")
    ax.set_ylabel("Tempo totale di generazione (secondi)")
    ax.set_title("Tempo di generazione vs indizi rimanenti per strategia")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_strategy_comparison(rows, out_path: Path):
    by_strategy = defaultdict(list)
    for row in rows:
        by_strategy[row["strategy"]].append(int(row["clues"]))
    strategies = sorted(by_strategy.keys())
    means = [mean(by_strategy[s]) for s in strategies]
    mins = [min(by_strategy[s]) for s in strategies]
    maxs = [max(by_strategy[s]) for s in strategies]

    fig, ax = plt.subplots(figsize=(7, 5))
    x_pos = range(len(strategies))
    ax.bar(x_pos, means, yerr=[[m - lo for m, lo in zip(means, mins)], [hi - m for m, hi in zip(means, maxs)]],
           capsize=6, color=["#4C72B0", "#DD8452", "#55A467"], alpha=0.85)
    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(strategies)
    ax.set_ylabel("Indizi rimanenti (media, min-max)")
    ax.set_title("Strategia di rimozione vs numero finale di indizi")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_method_comparison(rows, out_path: Path):
    by_method = defaultdict(list)
    for row in rows:
        by_method[row["method"]].append(float(row["avg_check_seconds"]))
    methods = sorted(by_method.keys())
    means = [mean(by_method[m]) * 1000 for m in methods]

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar(methods, means, color=["#4C72B0", "#DD8452"], alpha=0.85)
    ax.set_ylabel("Tempo medio per check di unicità (ms)")
    ax.set_title("Confronto metodi: counting vs solve-and-block")
    ax.grid(True, axis="y", alpha=0.3)
    for i, v in enumerate(means):
        ax.text(i, v, f"{v:.2f} ms", ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Genera i grafici dal benchmark CSV")
    parser.add_argument("--input", default="results/full_benchmark.csv")
    parser.add_argument("--output-dir", default="report/assets")
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[1]
    csv_path = project_root / args.input
    out_dir = project_root / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = read_csv(csv_path)
    if not rows:
        print(f"Nessuna riga trovata in {csv_path}")
        return 1

    plot_time_vs_clues(rows, out_dir / "plot_time_vs_clues.png")
    plot_strategy_comparison(rows, out_dir / "plot_strategy_comparison.png")
    plot_method_comparison(rows, out_dir / "plot_method_comparison.png")

    print(f"Grafici salvati in {out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
