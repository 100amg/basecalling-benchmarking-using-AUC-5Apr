#!/usr/bin/env python3
"""
Visualisation of AUC results across models, Q-scores, and sample types.

Reads the CSVs produced by auc_model_qscore_analysis.py and generates:
  1. Line plots — mean AUC vs Q-score, one line per basecaller+model combination
  2. AHEAD vs EpiC comparison — paired bar chart for every condition
  3. Top 10 combinations — horizontal bar chart ranked by mean AUC

Usage:
    python plot_auc_results.py <auc_results_dir>

Example:
    python plot_auc_results.py /Volumes/Amishi_SSD/bio_data/5Apr/new_EpiC_runs/auc_results
"""

import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ── Aesthetics ────────────────────────────────────────────────────────────────
PALETTE = {
    "dorado_fast": "#2196F3",
    "dorado_hac":  "#4CAF50",
    "dorado_sup":  "#FF5722",
    "guppy_hac":   "#9C27B0",
    "deepmod2_sup": "#FF9800",
}

AHEAD_COLOR  = "#1565C0"
EPIC_COLOR   = "#AD1457"
GRID_COLOR   = "#E0E0E0"
BG_COLOR     = "#FAFAFA"

plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "axes.spines.top":  False,
    "axes.spines.right": False,
    "axes.facecolor":   BG_COLOR,
    "figure.facecolor": "white",
    "axes.grid":        True,
    "grid.color":       GRID_COLOR,
    "grid.linewidth":   0.8,
    "axes.labelsize":   11,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "legend.fontsize":  9,
    "legend.framealpha": 0.9,
})

# ── Helpers ───────────────────────────────────────────────────────────────────

def load_data(results_dir: "/Volumes/Amishi_SSD/bio_data/5Apr/new_EpiC_runs/auc_results"):
    summary_path    = results_dir / "model_qscore_summary.csv"
    per_sample_path = results_dir / "per_sample_auc.csv"

    if not summary_path.exists():
        print(f"ERROR: {summary_path} not found. Run auc_model_qscore_analysis.py first.")
        sys.exit(1)

    summary    = pd.read_csv(summary_path)
    per_sample = pd.read_csv(per_sample_path) if per_sample_path.exists() else None

    # Create a combined label for grouping
    summary["condition"] = summary["basecaller"] + "_" + summary["model"]
    summary = summary.sort_values(["condition", "qscore"]).reset_index(drop=True)

    return summary, per_sample


def condition_color(condition: str) -> str:
    return PALETTE.get(condition, "#607D8B")


def save(fig, path: Path, name: str):
    out = path / name
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Saved → {out}")
    plt.close(fig)


# ── Plot 1: Mean AUC vs Q-score per model ─────────────────────────────────────

def plot_auc_vs_qscore(summary: pd.DataFrame, out_dir: Path):
    conditions = sorted(summary["condition"].unique())
    fig, ax = plt.subplots(figsize=(10, 6))

    for cond in conditions:
        sub = summary[summary["condition"] == cond].sort_values("qscore")
        color = condition_color(cond)
        ax.plot(sub["qscore"], sub["mean_auc"],
                marker="o", linewidth=2, markersize=7,
                color=color, label=cond.replace("_", " ").upper())

        # Error band (std)
        if "std_auc" in sub.columns:
            ax.fill_between(
                sub["qscore"],
                sub["mean_auc"] - sub["std_auc"].fillna(0),
                sub["mean_auc"] + sub["std_auc"].fillna(0),
                alpha=0.12, color=color
            )

    ax.set_xlabel("Q-score Filter")
    ax.set_ylabel("Mean AUC")
    ax.set_title("Mean AUC vs Q-score Filter\n(shaded = ±1 std across samples)",
                 fontsize=13, pad=12)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.set_ylim(bottom=max(0, summary["mean_auc"].min() - 0.05), top=1.02)
    ax.legend(loc="lower left", title="Basecaller + Model")
    fig.tight_layout()
    save(fig, out_dir, "1_auc_vs_qscore.png")


# ── Plot 2: AHEAD vs EpiC mean AUC comparison ─────────────────────────────────

def plot_ahead_vs_epic(summary: pd.DataFrame, out_dir: Path):
    df = summary.dropna(subset=["ahead_mean_auc", "epic_mean_auc"]).copy()
    df["label"] = df["condition"].str.replace("_", " ").str.upper() + "\nQ" + df["qscore"].astype(str)
    df = df.sort_values("mean_auc", ascending=False).reset_index(drop=True)

    x      = np.arange(len(df))
    width  = 0.38
    fig, ax = plt.subplots(figsize=(max(12, len(df) * 0.7), 6))

    bars_ahead = ax.bar(x - width/2, df["ahead_mean_auc"],
                        width, color=AHEAD_COLOR, alpha=0.85, label="AHEAD")
    bars_epic  = ax.bar(x + width/2, df["epic_mean_auc"],
                        width, color=EPIC_COLOR,  alpha=0.85, label="EpiC")

    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean AUC")
    ax.set_title("AHEAD vs EpiC Mean AUC by Model + Q-score\n(sorted by overall mean AUC)",
                 fontsize=13, pad=12)
    ax.set_ylim(0, 1.08)
    ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.legend(title="Sample Type", loc="upper right")

    # Value labels on bars
    for bar in list(bars_ahead) + list(bars_epic):
        h = bar.get_height()
        if h > 0.05:
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.01,
                    f"{h:.2f}", ha="center", va="bottom", fontsize=6.5)

    fig.tight_layout()
    save(fig, out_dir, "2_ahead_vs_epic_auc.png")


# ── Plot 3: Top 10 combinations ───────────────────────────────────────────────

def plot_top10(summary: pd.DataFrame, out_dir: Path):
    top10 = summary.nlargest(10, "mean_auc").reset_index(drop=True)

    top10["label"] = (
        top10["condition"].str.replace("_", " ").str.upper()
        + "  Q" + top10["qscore"].astype(str)
    )

    top10 = top10.sort_values("mean_auc")

    colors = [condition_color(c) for c in top10["condition"]]

    fig, ax = plt.subplots(figsize=(12, 7))

    bars = ax.barh(top10["label"], top10["mean_auc"],
                   color=colors, alpha=0.88, height=0.65)

    # Error bars
    errs = top10["std_auc"].fillna(0)
    ax.errorbar(top10["mean_auc"], top10["label"],
                xerr=errs,
                fmt="none", color="black", capsize=4, linewidth=1.2)

    # Value labels (fixed)
    for bar, val, err in zip(bars, top10["mean_auc"], errs):
        ax.text(val + err + 0.01,
                bar.get_y() + bar.get_height()/2,
                f"{val:.3f}", va="center", fontsize=9)

    ax.set_xlabel("Mean AUC")
    ax.set_title("Top 10 Basecaller + Model + Q-score Combinations",
                 fontsize=13, pad=12)

    ax.set_xlim(0, min(1.15, top10["mean_auc"].max() + 0.12))
    ax.axvline(0.5, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)

    # Legend outside
    seen = {}
    for cond, col in zip(top10["condition"], colors):
        if cond not in seen:
            seen[cond] = col

    handles = [plt.Rectangle((0,0),1,1, color=c, alpha=0.88)
               for c in seen.values()]
    labels  = [k.replace("_", " ").upper() for k in seen.keys()]

    ax.legend(handles, labels,
              title="Basecaller + Model",
              loc="center left",
              bbox_to_anchor=(1.02, 0.5),
              fontsize=8)

    fig.subplots_adjust(right=0.75)
    fig.tight_layout()
    save(fig, out_dir, "3_top10_combinations.png")

# ── Plot 4: Per-sample AUC distribution (box plots) ──────────────────────────

def plot_sample_distribution(per_sample: pd.DataFrame, out_dir: Path):
    """Box plot showing AUC distribution across samples per condition."""
    per_sample = per_sample.copy()
    per_sample["condition"] = (
        per_sample["basecaller"] + "_" + per_sample["model"]
        + "_Q" + per_sample["qscore"].astype(str)
    )

    # Sort conditions by median AUC
    order = (per_sample.groupby("condition")["auc"]
             .median().sort_values(ascending=False).index.tolist())

    data = [per_sample[per_sample["condition"] == c]["auc"].values for c in order]
    colors = [condition_color("_".join(c.split("_")[:2])) for c in order]

    fig, ax = plt.subplots(figsize=(max(14, len(order) * 0.55), 8))

    bp = ax.boxplot(data, patch_artist=True,
                    medianprops=dict(color="black", linewidth=2))

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=60, ha="right", fontsize=7)

    ax.tick_params(axis="x", pad=8)

    ax.set_ylabel("AUC (per sample)")
    ax.set_title("Per-sample AUC Distribution Across All Conditions",
                 fontsize=13, pad=12)

    ax.set_ylim(0, 1.05)
    ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)

    # Layout fixes
    fig.subplots_adjust(bottom=0.3, right=0.8)

    fig.tight_layout()
    save(fig, out_dir, "4_per_sample_distribution.png")


# ── Main ──────────────────────────────────────────────────────────────────────

def main(results_dir: str):
    results_dir = Path(results_dir)
    if not results_dir.is_dir():
        print(f"ERROR: {results_dir} not found"); sys.exit(1)

    out_dir = results_dir / "plots"
    out_dir.mkdir(exist_ok=True)
    print(f"Saving plots to: {out_dir}\n")

    summary, per_sample = load_data(results_dir)

    print("Generating Plot 1: AUC vs Q-score per model...")
    plot_auc_vs_qscore(summary, out_dir)

    print("Generating Plot 2: AHEAD vs EpiC comparison...")
    plot_ahead_vs_epic(summary, out_dir)

    print("Generating Plot 3: Top 10 combinations...")
    plot_top10(summary, out_dir)

    if per_sample is not None:
        print("Generating Plot 4: Per-sample distribution...")
        plot_sample_distribution(per_sample, out_dir)

    print(f"\nAll plots saved to {out_dir}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_auc_results.py <auc_results_dir>")
        sys.exit(1)
    main(sys.argv[1])
