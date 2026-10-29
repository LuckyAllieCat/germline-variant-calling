# Benchmark both call sets against GIAB with RTG vcfeval (haplotype-aware
# comparison, so different representations of the same variant still match).

rule rtg_sdf:
    input:
        REF,
    output:
        directory(f"data/ref/{CHROM}_sdf"),
    log:
        "logs/rtg_format.log",
    shell:
        "rtg format -o {output} {input} > {log} 2>&1"


rule vcfeval:
    input:
        truth=rules.truth_set.output.vcf,
        bed=rules.truth_set.output.bed,
        query="results/vcf/{caller}.filtered.vcf.gz",
        sdf=rules.rtg_sdf.output,
    output:
        directory("results/benchmark/{caller}_{mode}"),
    params:
        # "all" mode ignores the FILTER column, i.e. scores the unfiltered calls
        extra=lambda w: "--all-records" if w.mode == "all" else "",
    log:
        "logs/vcfeval/{caller}_{mode}.log",
    threads: 2
    shell:
        """
        rm -rf {output}
        RTG_MEM=2g rtg vcfeval -b {input.truth} -c {input.query} -e {input.bed} -t {input.sdf} \
          -o {output} --ref-overlap --threads {threads} {params.extra} > {log} 2>&1
        """


rule benchmark_summary:
    input:
        expand("results/benchmark/{caller}_{mode}", caller=CALLERS, mode=MODES),
    output:
        tsv="results/tables/benchmark_summary.tsv",
        fig="figures/benchmark_comparison.png",
    params:
        evals=" ".join(f"{c}:{m}:results/benchmark/{c}_{m}" for c in CALLERS for m in MODES),
    shell:
        "python scripts/benchmark_summary.py --eval {params.evals} --out-tsv {output.tsv} --out-fig {output.fig}"


rule error_profile:
    input:
        evals=expand("results/benchmark/{caller}_pass", caller=CALLERS),
        depth=rules.alignment_qc.output.perbase,
        repeats="data/strat/repeats.bed",
        lowmap="data/strat/lowmap.bed",
    output:
        tsv="results/tables/error_profile.tsv",
        sites="results/benchmark/benchmark_sites.tsv",
        fig="figures/error_profile.png",
    params:
        evals=" ".join(f"{c}:results/benchmark/{c}_pass" for c in CALLERS),
    shell:
        """
        python scripts/error_profile.py --eval {params.evals} --depth {input.depth} \
          --repeats {input.repeats} --lowmap {input.lowmap} \
          --out-tsv {output.tsv} --out-sites {output.sites} --out-fig {output.fig}
        """
