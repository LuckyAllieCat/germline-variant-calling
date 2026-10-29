"""Small, dependency-free helpers shared by the analysis scripts."""
import bisect
import gzip

# Fixed colours per caller (categorical slots 1 and 2), used by every figure
CALLER_COLORS = {"gatk": "#2a78d6", "bcftools": "#eb6834"}
CALLER_LABELS = {"gatk": "GATK HaplotypeCaller", "bcftools": "bcftools call"}


def read_vcf(path):
    """Yield (chrom, pos, ref, alts) for every record of a (b)gzipped VCF."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            yield f[0], int(f[1]), f[3], f[4].split(",")


def var_class(ref, alts):
    """SNV if every allele is one base long, otherwise INDEL (incl. complex)."""
    alts = [a for a in alts if a not in ("*", ".")]
    if len(ref) == 1 and all(len(a) == 1 for a in alts):
        return "SNV"
    return "INDEL"


def indel_size(ref, alts):
    """Signed length change of the largest-effect allele (+ insertion, - deletion)."""
    sizes = [len(a) - len(ref) for a in alts if a not in ("*", ".")]
    return max(sizes, key=abs) if sizes else 0


def size_bin(size):
    """Group indel sizes into readable bins, deletions first."""
    if size == 0:
        return "complex"
    kind = "ins" if size > 0 else "del"
    n = abs(size)
    label = "1" if n == 1 else "2-5" if n <= 5 else "6-15" if n <= 15 else ">15"
    return f"{kind} {label}"


SIZE_BINS = ["del >15", "del 6-15", "del 2-5", "del 1", "ins 1", "ins 2-5", "ins 6-15", "ins >15"]


def metrics(tp_base, tp_call, fp, fn):
    """Precision, recall and F1 the way vcfeval/hap.py define them."""
    precision = tp_call / (tp_call + fp) if tp_call + fp else float("nan")
    recall = tp_base / (tp_base + fn) if tp_base + fn else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else float("nan")
    return precision, recall, f1


class IntervalSet:
    """Merged BED intervals with fast point lookup (0-based half-open)."""

    def __init__(self, intervals):
        self.starts, self.ends = {}, {}
        by_chrom = {}
        for chrom, s, e in intervals:
            by_chrom.setdefault(chrom, []).append((s, e))
        for chrom, ivs in by_chrom.items():
            merged = []
            for s, e in sorted(ivs):
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.starts[chrom] = [s for s, _ in merged]
            self.ends[chrom] = [e for _, e in merged]

    @classmethod
    def from_bed(cls, path):
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            rows = (line.split("\t")[:3] for line in fh if line.strip() and not line.startswith(("#", "track")))
            return cls((c, int(s), int(e)) for c, s, e in rows)

    def contains(self, chrom, pos):
        """True if 1-based VCF position `pos` falls inside an interval."""
        starts = self.starts.get(chrom)
        if not starts:
            return False
        i = bisect.bisect_right(starts, pos - 1) - 1
        return i >= 0 and pos - 1 < self.ends[chrom][i]


class DepthTrack:
    """Per-base depth from a mosdepth *.per-base.bed.gz file."""

    def __init__(self, path):
        self.starts, self.depths = {}, {}
        with gzip.open(path, "rt") as fh:
            for line in fh:
                chrom, s, _, d = line.split("\t")
                self.starts.setdefault(chrom, []).append(int(s))
                self.depths.setdefault(chrom, []).append(int(d))

    def at(self, chrom, pos):
        starts = self.starts.get(chrom)
        if not starts:
            return 0
        return self.depths[chrom][bisect.bisect_right(starts, pos - 1) - 1]
