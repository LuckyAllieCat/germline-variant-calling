"""Variant prioritization: consequence -> impact -> population frequency -> ClinVar.

This demonstrates a filtering workflow. It is NOT a clinical interpretation and
assigns no ACMG classification; ClinVar labels are shown as reported.
"""
import argparse
import gzip
import math
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from vcfutils import CALLER_LABELS, IntervalSet

# Sequence Ontology terms used by bcftools csq, mapped to VEP-style impact levels
IMPACT = {
    "frameshift": "HIGH", "stop_gained": "HIGH", "stop_lost": "HIGH", "start_lost": "HIGH",
    "splice_acceptor": "HIGH", "splice_donor": "HIGH", "transcript_ablation": "HIGH",
    "missense": "MODERATE", "inframe_insertion": "MODERATE", "inframe_deletion": "MODERATE",
    "inframe_altering": "MODERATE",
    "synonymous": "LOW", "splice_region": "LOW", "stop_retained": "LOW", "start_retained": "LOW",
}
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
CODING_TERMS = {t for t in IMPACT} | {"coding_sequence"}

# ClinVar review status -> gold stars (as displayed on the ClinVar website)
STARS = {
    "practice_guideline": 4, "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1, "criteria_provided,_single_submitter": 1,
}
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}


def parse_bcsq(bcsq):
    """Split a bcftools csq BCSQ string into one dict per transcript consequence.

    Format: consequence|gene|transcript|biotype|strand|aa_change|dna_change
    Entries starting with '@' point to a compound record and are skipped.
    """
    out = []
    if not bcsq or bcsq == ".":
        return out
    for entry in bcsq.split(","):
        if entry.startswith("@"):
            continue
        f = entry.split("|")
        terms = [t.lstrip("*") for t in f[0].split("&")]
        impact = max((IMPACT.get(t, "MODIFIER") for t in terms), key=IMPACT_RANK.get)
        out.append(dict(consequence="&".join(terms), impact=impact,
                        gene=f[1] if len(f) > 1 else "", transcript=f[2] if len(f) > 2 else "",
                        biotype=f[3] if len(f) > 3 else "", protein=f[5] if len(f) > 5 else "",
                        dna=f[6] if len(f) > 6 else "",
                        coding=any(t in CODING_TERMS for t in terms)))
    return out


def pick(csqs, mane):
    """Worst consequence over all transcripts, plus the MANE Select one (if any)."""
    if not csqs:
        return None, None
    worst = max(csqs, key=lambda c: (IMPACT_RANK[c["impact"]], c["biotype"] == "protein_coding",
                                     c["transcript"] in mane))
    mane_csq = next((c for c in csqs if c["transcript"] in mane), None)
    return worst, mane_csq


def mane_transcripts(gff):
    """Transcript IDs tagged MANE_Select in the Ensembl GFF3."""
    ids = set()
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            if "MANE_Select" in line:
                m = re.search(r"ID=transcript:([^;]+)", line)
                if m:
                    ids.add(m.group(1))
    return ids


def clinvar_class(clnsig):
    """First (primary) ClinVar significance term, e.g. 'Pathogenic|risk_factor' -> 'Pathogenic'."""
    if not clnsig or clnsig == ".":
        return "not_in_ClinVar"
    return re.split(r"[|,]", clnsig)[0]


def stars(revstat):
    return STARS.get(revstat, 0) if revstat and revstat != "." else ""


def load_gnomad(path, prefix):
    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in ["AF", "AF_grpmax", "nhomalt"]:
        df[col] = pd.to_numeric(df[col].replace(".", None), errors="coerce")
    df = df.drop_duplicates(["CHROM", "POS", "REF", "ALT"])
    df["POS"] = df["POS"].astype(int)
    return df.rename(columns={c: f"{prefix}_{c}" for c in ["FILTER", "AF", "AF_grpmax", "nhomalt"]})


def choose_caller(benchmark, requested):
    """'auto' = the caller with the higher mean of SNV and indel F1 (PASS calls)."""
    if requested != "auto":
        return requested
    b = pd.read_csv(benchmark, sep="\t")
    b = b[b["filter_mode"] == "pass"]
    return b.groupby("caller")["F1"].mean().idxmax()


def max_af(row):
    vals = [row.get(c) for c in ["gnomADe_AF", "gnomADe_AF_grpmax", "gnomADg_AF", "gnomADg_AF_grpmax"]]
    vals = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    return max(vals) if vals else 0.0  # absent from gnomAD -> treated as AF 0


def assign_tier(row, keep):
    cls = row["clinvar_significance"]
    if cls in PLP:
        return "1: ClinVar P/LP"
    if cls in keep:
        return "2: ClinVar VUS/conflicting"
    return "3: rare, not in ClinVar"


def annotate(table, mane, gnomad_e, gnomad_g):
    df = pd.read_csv(table, sep="\t", dtype=str)
    df["POS"] = df["POS"].astype(int)
    rows = []
    for rec in df.to_dict("records"):
        worst, mane_csq = pick(parse_bcsq(rec["BCSQ"]), mane)
        rep = mane_csq or worst
        rec.update(gene=worst["gene"] if worst else "",
                   worst_consequence=worst["consequence"] if worst else "intergenic",
                   impact=worst["impact"] if worst else "MODIFIER",
                   coding=bool(worst and worst["coding"]),
                   transcript=rep["transcript"] if rep else "",
                   mane_select=bool(mane_csq),
                   consequence=rep["consequence"] if rep else "",  # on MANE if available
                   protein_change=rep["protein"] if rep else "",
                   clinvar_significance=clinvar_class(rec["CLNSIG"]),
                   clinvar_stars=stars(rec["CLNREVSTAT"]))
        rows.append(rec)
    out = pd.DataFrame(rows)
    out = out.merge(gnomad_e, on=["CHROM", "POS", "REF", "ALT"], how="left")
    out = out.merge(gnomad_g, on=["CHROM", "POS", "REF", "ALT"], how="left")
    out["gnomAD_max_AF"] = out.apply(max_af, axis=1)
    return out


def truth_status(df, sites_path, caller):
    """Is the call a GIAB true positive? (position-level match on the benchmark output)"""
    s = pd.read_csv(sites_path, sep="\t")
    s = s[s["caller"] == caller]
    fp, tp = set(s.loc[s["outcome"] == "FP", "pos"]), set(s.loc[s["outcome"] == "TP", "pos"])
    return df["POS"].map(lambda p: "FP" if p in fp else "TP" if p in tp else "outside benchmark")


def plot(funnel, prio, caller, out):
    fig, axes = plt.subplots(1, 2, figsize=(14, max(4.2, 0.28 * len(prio) + 1.6)),
                             gridspec_kw={"width_ratios": [1, 1.15]})
    # A: the filtering funnel (one series, one colour)
    ax = axes[0]
    ys = range(len(funnel))
    ax.barh(list(ys), funnel["n"], color="#2a78d6", height=0.62)
    for y, (label, n) in enumerate(zip(funnel["step"], funnel["n"])):
        ax.annotate(f"{n:,}", (n, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=9)
    ax.set_yticks(list(ys), funnel["step"], fontsize=9)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.8, funnel["n"].max() * 4)
    ax.set_xlabel("variants (log scale)")
    ax.set_title(f"A  Prioritization funnel ({CALLER_LABELS[caller]} PASS calls)", loc="left", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    # B: prioritized variants by population frequency, coloured/shaped by ClinVar tier
    ax = axes[1]
    styles = {"1: ClinVar P/LP": ("#e87ba4", "D"), "2: ClinVar VUS/conflicting": ("#1baf7a", "s"),
              "3: rare, not in ClinVar": ("#898781", "o")}
    if prio.empty:
        ax.text(0.5, 0.5, "No variants passed all filters", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    else:
        p = prio.reset_index(drop=True)
        labels = [f"{g} {pc or c}" for g, pc, c in zip(p["gene"], p["protein_change"], p["consequence"])]
        x = p["gnomAD_max_AF"].clip(lower=1e-6)
        for tier, (col, mk) in styles.items():
            m = p["tier"] == tier
            if m.any():
                ax.scatter(x[m], p.index[m], color=col, marker=mk, s=60, edgecolor="white",
                           linewidth=1.2, zorder=3, label=tier)
        ax.set_yticks(range(len(p)), labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_xlim(5e-7, 0.05)
        ax.axvline(1e-6, color="#c3c2b7", lw=0.8, ls=":")
        ax.set_xlabel("max gnomAD v4.1 AF (absent plotted at 1e-6)")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        ax.grid(axis="x", color="#e5e4df", linewidth=0.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("B  Prioritized variants (research demo, not clinical)", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tables", nargs="+", required=True, help="caller:annotated.tsv")
    ap.add_argument("--caller", default="auto")
    ap.add_argument("--gnomad-exomes", required=True)
    ap.add_argument("--gnomad-genomes", required=True)
    ap.add_argument("--benchmark", required=True)
    ap.add_argument("--benchmark-sites", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--repeats", required=True, help="GIAB tandem repeat/homopolymer BED")
    ap.add_argument("--lowmap", required=True, help="GIAB low-mappability BED")
    ap.add_argument("--max-af", type=float, default=0.01)
    ap.add_argument("--impacts", nargs="+", default=["HIGH", "MODERATE"])
    ap.add_argument("--clinvar-keep", nargs="+", default=[])
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-coding", required=True)
    ap.add_argument("--out-funnel", required=True)
    ap.add_argument("--out-fig", required=True)
    a = ap.parse_args()

    caller = choose_caller(a.benchmark, a.caller)
    table = dict(t.split(":", 1) for t in a.tables)[caller]
    df = annotate(table, mane_transcripts(a.gff),
                  load_gnomad(a.gnomad_exomes, "gnomADe"), load_gnomad(a.gnomad_genomes, "gnomADg"))
    df["caller"] = caller
    df["GIAB_status"] = truth_status(df, a.benchmark_sites, caller)
    # Flag calls in regions where short reads are known to be unreliable
    rep, low = IntervalSet.from_bed(a.repeats), IntervalSet.from_bed(a.lowmap)
    df["difficult_region"] = [
        ";".join(n for n, s in [("repeat", rep), ("low_mappability", low)] if s.contains(c, p)) or "no"
        for c, p in zip(df["CHROM"], df["POS"])]

    # Step-by-step filters; each count feeds the funnel figure
    in_gene = df["gene"] != ""
    coding = in_gene & df["coding"]
    impact = coding & df["impact"].isin(a.impacts)
    rare = impact & (df["gnomAD_max_AF"] < a.max_af)
    not_benign = rare & ~df["clinvar_significance"].isin(BENIGN)
    plp_any = df["clinvar_significance"].isin(PLP)  # reported regardless of frequency
    keep = not_benign | plp_any
    funnel = pd.DataFrame({"step": ["PASS variants", "in a gene", "coding / splice site",
                                    f"{'/'.join(a.impacts)} impact", f"rare (gnomAD AF < {a.max_af:g})",
                                    "not ClinVar benign", "prioritized (incl. any ClinVar P/LP)"],
                           "n": [len(df), in_gene.sum(), coding.sum(), impact.sum(), rare.sum(),
                                 not_benign.sum(), keep.sum()]})

    prio = df[keep].copy()
    prio["tier"] = prio.apply(assign_tier, axis=1, keep=set(a.clinvar_keep))
    prio = prio.sort_values(["tier", "gnomAD_max_AF", "POS"])
    cols = ["tier", "gene", "CHROM", "POS", "REF", "ALT", "GT", "DP", "QUAL", "worst_consequence", "impact", "consequence",
            "protein_change", "transcript", "mane_select", "clinvar_significance", "clinvar_stars",
            "CLNDN", "CLINVAR_ID", "gnomADe_AF", "gnomADg_AF", "gnomADe_AF_grpmax", "gnomADg_AF_grpmax",
            "gnomAD_max_AF", "GIAB_status", "difficult_region", "caller"]
    prio[cols].rename(columns={"CLNDN": "clinvar_condition"}).to_csv(
        a.out_table, sep="\t", index=False, float_format="%.3g")
    coding_cols = [c for c in cols if c != "tier"]
    df[coding][coding_cols].to_csv(a.out_coding, sep="\t", index=False, float_format="%.3g")
    funnel.to_csv(a.out_funnel, sep="\t", index=False)
    plot(funnel, prio, caller, a.out_fig)


if __name__ == "__main__":
    main()
