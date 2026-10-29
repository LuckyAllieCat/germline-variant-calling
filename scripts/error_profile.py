"""Where does each caller go wrong? Stratify benchmark outcomes by read depth,
indel size and GIAB difficult-region context (PASS calls only)."""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from vcfutils import (CALLER_COLORS, CALLER_LABELS, SIZE_BINS, DepthTrack, IntervalSet,
                      indel_size, read_vcf, size_bin, var_class)

DEPTH_BINS = [0, 10, 20, 30, 40, 10_000]
DEPTH_LABELS = ["<10", "10-19", "20-29", "30-39", "40+"]


def site_table(caller, eval_dir, depth, repeats, lowmap):
    """One row per benchmark record with its outcome and genomic context."""
    rows = []
    for outcome, fname in [("TP", "tp-baseline"), ("FN", "fn"), ("FP", "fp")]:
        for chrom, pos, ref, alts in read_vcf(os.path.join(eval_dir, f"{fname}.vcf.gz")):
            vclass = var_class(ref, alts)
            rows.append(dict(caller=caller, outcome=outcome, variant_class=vclass, pos=pos,
                             depth=depth.at(chrom, pos),
                             size_bin=size_bin(indel_size(ref, alts)) if vclass == "INDEL" else "",
                             repeat=repeats.contains(chrom, pos),
                             lowmap=lowmap.contains(chrom, pos)))
    return pd.DataFrame(rows)


def recall(df):
    tp, fn = (df["outcome"] == "TP").sum(), (df["outcome"] == "FN").sum()
    return tp / (tp + fn) if tp + fn else float("nan"), int(tp + fn)


def summarise(sites):
    """Long-format table: caller, stratum, recall, error counts."""
    out = []
    for caller, d in sites.groupby("caller"):
        d = d.assign(depth_bin=pd.cut(d["depth"], DEPTH_BINS, right=False, labels=DEPTH_LABELS))
        strata = [("all", "all", d)]
        strata += [("depth", str(b), g) for b, g in d.groupby("depth_bin", observed=True)]
        strata += [("indel_size", b, g) for b, g in d[d["variant_class"] == "INDEL"].groupby("size_bin")]
        for name, col in [("tandem_repeat_or_homopolymer", "repeat"), ("low_mappability", "lowmap")]:
            strata += [(name, str(v), g) for v, g in d.groupby(col)]
        for kind, level, g in strata:
            for vclass in ["SNV", "INDEL"]:
                gv = g[g["variant_class"] == vclass]
                if gv.empty:
                    continue
                r, n = recall(gv)
                out.append(dict(caller=caller, stratum=kind, level=level, variant_class=vclass,
                                truth_n=n, recall=r, FN=int((gv["outcome"] == "FN").sum()),
                                FP=int((gv["outcome"] == "FP").sum())))
    return pd.DataFrame(out)


def plot(table, sites, out):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    callers = ["gatk", "bcftools"]
    width = 0.38
    # A: recall by depth (SNVs + indels pooled)
    ax = axes[0]
    for j, c in enumerate(callers):
        d = sites[(sites["caller"] == c) & (sites["outcome"] != "FP")]
        d = d.assign(b=pd.cut(d["depth"], DEPTH_BINS, right=False, labels=DEPTH_LABELS))
        rec = d.groupby("b", observed=False)["outcome"].apply(lambda s: (s == "TP").mean())
        ax.plot(range(len(DEPTH_LABELS)), rec.values, marker="o", lw=2, ms=7,
                color=CALLER_COLORS[c], label=CALLER_LABELS[c], markeredgecolor="white")
    n = d.groupby("b", observed=False).size()
    ax.set_xticks(range(len(DEPTH_LABELS)), [f"{l}\n(n={k:,})" for l, k in zip(DEPTH_LABELS, n.values)], fontsize=8)
    ax.set_xlabel("read depth at truth variant")
    ax.set_ylabel("recall")
    ax.set_title("A  Recall by read depth", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    # B: indel recall by indel size
    ax = axes[1]
    sub = table[(table["stratum"] == "indel_size")]
    bins = [b for b in SIZE_BINS if b in set(sub["level"])]
    for j, c in enumerate(callers):
        s = sub[sub["caller"] == c].set_index("level").reindex(bins)
        ax.bar([i + (j - 0.5) * width for i in range(len(bins))], s["recall"], width - 0.04,
               color=CALLER_COLORS[c], label=CALLER_LABELS[c])
    ns = sub[sub["caller"] == "gatk"].set_index("level").reindex(bins)["truth_n"]
    ax.set_xticks(range(len(bins)), [f"{b.replace(' ', chr(10))}\n(n={int(k)})" for b, k in zip(bins, ns)], fontsize=7.5)
    ax.set_ylabel("indel recall")
    ax.set_ylim(0, 1.05)
    ax.set_title("B  Indel recall by size (bp)", loc="left", fontsize=11)
    # C: share of errors that fall in difficult regions vs share of truth variants there
    ax = axes[2]
    truth = sites[(sites["caller"] == "gatk") & (sites["outcome"] != "FP")]
    groups = [("truth\nvariants", None, truth)]
    for c in callers:
        for o in ["FP", "FN"]:
            groups.append((f"{CALLER_LABELS[c].split()[0]}\n{o}", c, sites[(sites["caller"] == c) & (sites["outcome"] == o)]))
    xs = range(len(groups))
    vals = [g["repeat"].mean() * 100 if len(g) else 0 for _, _, g in groups]
    cols = ["#898781" if c is None else CALLER_COLORS[c] for _, c, _ in groups]
    ax.bar(xs, vals, 0.65, color=cols)
    for x, v, (_, _, g) in zip(xs, vals, groups):
        ax.annotate(f"{v:.0f}%\nn={len(g):,}", (x, v), xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=8, color="#3d3d3a")
    ax.set_xticks(list(xs), [g[0] for g in groups], fontsize=8.5)
    ax.set_ylabel("% in tandem repeats / homopolymers")
    ax.set_ylim(0, 118)
    ax.set_yticks(range(0, 101, 20))
    ax.set_title("C  Share in tandem repeats / homopolymers", loc="left", fontsize=11)
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
        a.grid(axis="y", color="#e5e4df", linewidth=0.8)
        a.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out, dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eval", nargs="+", required=True, help="caller:vcfeval_dir (PASS mode)")
    ap.add_argument("--depth", required=True, help="mosdepth per-base bed.gz")
    ap.add_argument("--repeats", required=True)
    ap.add_argument("--lowmap", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-sites", required=True)
    ap.add_argument("--out-fig", required=True)
    args = ap.parse_args()
    depth = DepthTrack(args.depth)
    repeats, lowmap = IntervalSet.from_bed(args.repeats), IntervalSet.from_bed(args.lowmap)
    sites = pd.concat([site_table(*e.split(":", 1), depth, repeats, lowmap) for e in args.eval])
    sites.to_csv(args.out_sites, sep="\t", index=False)
    table = summarise(sites)
    table.to_csv(args.out_tsv, sep="\t", index=False, float_format="%.4f")
    plot(table, sites, args.out_fig)


if __name__ == "__main__":
    main()
