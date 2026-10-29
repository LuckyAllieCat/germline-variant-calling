"""Unit tests for the small pieces of logic the results depend on."""
import gzip
import math

import pandas as pd

from prioritize import assign_tier, clinvar_class, max_af, parse_bcsq, pick, stars
from vcfutils import IntervalSet, indel_size, metrics, read_vcf, size_bin, var_class


def test_var_class():
    assert var_class("A", ["G"]) == "SNV"
    assert var_class("A", ["G", "T"]) == "SNV"
    assert var_class("A", ["AT"]) == "INDEL"
    assert var_class("AT", ["A", "G"]) == "INDEL"  # multi-allelic with one indel allele
    assert var_class("A", ["*", "G"]) == "SNV"     # spanning-deletion allele ignored


def test_indel_size_and_bins():
    assert indel_size("A", ["ATTT"]) == 3
    assert indel_size("ATTTTT", ["A"]) == -5
    assert size_bin(1) == "ins 1"
    assert size_bin(-20) == "del >15"
    assert size_bin(-6) == "del 6-15"


def test_metrics():
    p, r, f1 = metrics(tp_base=90, tp_call=90, fp=10, fn=10)
    assert math.isclose(p, 0.9) and math.isclose(r, 0.9) and math.isclose(f1, 0.9)


def test_interval_set_merges_and_uses_vcf_coordinates():
    s = IntervalSet([("chr20", 10, 20), ("chr20", 15, 30)])
    assert s.contains("chr20", 11)       # 1-based 11 = 0-based 10 (inside)
    assert s.contains("chr20", 30)       # 0-based 29, last base
    assert not s.contains("chr20", 10)   # 0-based 9, just before
    assert not s.contains("chr20", 31)
    assert not s.contains("chr1", 15)


def test_read_vcf(tmp_path):
    p = tmp_path / "x.vcf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        fh.write("chr20\t100\t.\tA\tG,T\t50\tPASS\t.\n")
    assert list(read_vcf(str(p))) == [("chr20", 100, "A", ["G", "T"])]


def test_parse_bcsq_and_pick_prefers_mane():
    bcsq = ("missense|HNF4A|ENST1|protein_coding|+|117V>117M|44401260G>A,"
            "stop_gained|HNF4A|ENST2|protein_coding|+|90Q>90*|44401260G>A,@44401300")
    csqs = parse_bcsq(bcsq)
    assert [c["impact"] for c in csqs] == ["MODERATE", "HIGH"]
    worst, mane = pick(csqs, mane={"ENST1"})
    assert worst["transcript"] == "ENST2"           # filtering uses the worst consequence
    assert mane["protein"] == "117V>117M"           # reporting uses the MANE transcript
    assert parse_bcsq(".") == []


def test_compound_and_starred_terms():
    c = parse_bcsq("*missense&inframe_altering|G|T|protein_coding|+|1A>1B|1A>B")[0]
    assert c["consequence"] == "missense&inframe_altering" and c["impact"] == "MODERATE"
    assert parse_bcsq("intron|G|T|protein_coding")[0]["coding"] is False


def test_clinvar_helpers():
    assert clinvar_class("Pathogenic|risk_factor") == "Pathogenic"
    assert clinvar_class(".") == "not_in_ClinVar"
    assert stars("reviewed_by_expert_panel") == 3
    assert stars(".") == ""


def test_max_af_and_tiers():
    assert max_af({"gnomADe_AF": 0.001, "gnomADg_AF_grpmax": 0.02}) == 0.02
    assert max_af({"gnomADe_AF": float("nan")}) == 0.0  # absent = never seen
    keep = {"Uncertain_significance"}
    assert assign_tier(pd.Series({"clinvar_significance": "Likely_pathogenic"}), keep).startswith("1")
    assert assign_tier(pd.Series({"clinvar_significance": "Uncertain_significance"}), keep).startswith("2")
    assert assign_tier(pd.Series({"clinvar_significance": "not_in_ClinVar"}), keep).startswith("3")
