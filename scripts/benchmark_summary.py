"""Precision/recall/F1 per caller, per filter mode and per variant class.

Counts come from the vcfeval "split" output VCFs, so SNVs and indels are
scored from one haplotype-aware comparison instead of two separate runs.
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from vcfutils import CALLER_COLORS, CALLER_LABELS, metrics, read_vcf, var_class


def count_by_class(eval_dir):
    counts = {}
    for name in ["tp-baseline", "tp", "fp", "fn"]:
        for _, _, ref, alts in read_vcf(os.path.join(eval_dir, f"{name}.vcf.gz")):
            key = (var_class(ref, alts), name)
            counts[key] = counts.get(key, 0) + 1
    rows = []
    for vclass in ["SNV", "INDEL"]:
        c = {n: counts.get((vclass, n), 0) for n in ["tp-baseline", "tp", "fp", "fn"]}
        p, r, f1 = metrics(c["tp-baseline"], c["tp"], c["fp"], c["fn"])
        rows.append(dict(variant_class=vclass, truth_total=c["tp-baseline"] + c["fn"],
                         TP_baseline=c["tp-baseline"], TP_call=c["tp"], FP=c["fp"], FN=c["fn"],
                         precision=p, recall=r, F1=f1))
    return rows


def plot(df, out):
    pas = df[df["filter_mode"] == "pass"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    metrics_order = ["precision", "recall", "F1"]
    # Panels 1-2: dot plot of each metric, one dot per caller (PASS calls)
    for ax, vclass in zip(axes[:2], ["SNV", "INDEL"]):
        sub = pas[pas["variant_class"] == vclass]
        lo = min(sub[metrics_order].min().min(), 0.95)
        for j, caller in enumerate(["gatk", "bcftools"]):
            row = sub[sub["caller"] == caller].iloc[0]
            ys = [i + (j - 0.5) * 0.28 for i in range(3)]
            xs = [row[m] for m in metrics_order]
            ax.scatter(xs, ys, s=70, color=CALLER_COLORS[caller], label=CALLER_LABELS[caller],
                       zorder=3, edgecolor="white", linewidth=1.5)
            for x, y in zip(xs, ys):
                ax.annotate(f"{x:.4f}", (x, y), xytext=(7, 0), textcoords="offset points",
                            va="center", fontsize=8.5, color="#3d3d3a")
        ax.set_yticks(range(3), metrics_order)
        ax.invert_yaxis()
        ax.set_xlim(lo - (1 - lo) * 0.15, 1.0 + (1 - lo) * 0.25)
        n = int(sub["truth_total"].iloc[0])
        ax.set_title(f"{vclass}s  (truth n = {n:,})", loc="left", fontsize=11)
        ax.grid(axis="x", color="#e5e4df", linewidth=0.8)
        ax.set_axisbelow(True)
    # Panel 3: how hard filtering moves each caller in precision/recall space
    ax = axes[2]
    for caller in ["gatk", "bcftools"]:
        for vclass, marker in [("SNV", "o"), ("INDEL", "s")]:
            sub = df[(df["caller"] == caller) & (df["variant_class"] == vclass)].set_index("filter_mode")
            a, b = sub.loc["all"], sub.loc["pass"]
            ax.annotate("", xy=(b["recall"], b["precision"]), xytext=(a["recall"], a["precision"]),
                        arrowprops=dict(arrowstyle="->", color=CALLER_COLORS[caller], lw=1.5))
            ax.scatter([a["recall"]], [a["precision"]], marker=marker, s=50, facecolor="white",
                       edgecolor=CALLER_COLORS[caller], linewidth=1.5, zorder=3)
            ax.scatter([b["recall"]], [b["precision"]], marker=marker, s=60,
                       color=CALLER_COLORS[caller], edgecolor="white", linewidth=1.5, zorder=3,
                       label=f"{CALLER_LABELS[caller].split()[0]} {vclass}")
    ax.margins(0.12)
    ax.set_xlabel("recall")
    ax.set_ylabel("precision")
    ax.set_title("Effect of hard filters (hollow = unfiltered)", loc="left", fontsize=11)
    ax.grid(color="#e5e4df", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    # one shared caller legend under the dot plots
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower left", ncol=2, frameon=False, fontsize=9, bbox_to_anchor=(0.06, 0.0))
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle("HG002 vs GIAB v4.2.1 benchmark, high-confidence regions", x=0.01, ha="left", fontsize=12.5)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out, dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eval", nargs="+", required=True, help="caller:mode:vcfeval_dir")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-fig", required=True)
    args = ap.parse_args()
    rows = []
    for item in args.eval:
        caller, mode, path = item.split(":", 2)
        for r in count_by_class(path):
            rows.append(dict(caller=caller, filter_mode=mode, **r))
    df = pd.DataFrame(rows)
    df.to_csv(args.out_tsv, sep="\t", index=False, float_format="%.4f")
    plot(df, args.out_fig)


if __name__ == "__main__":
    main()
