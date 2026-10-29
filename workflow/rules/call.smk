# Variant calling with GATK HaplotypeCaller and bcftools, scattered over
# region chunks, then hard-filtered and left-normalised.

rule gatk_haplotypecaller:
    input:
        bam=rules.align.output.bam,
        fa=REF,
        dict=rules.reference.output.dict,
    output:
        temp("results/vcf/chunks/gatk.{chunk}.vcf.gz"),
    params:
        reg=lambda w: CHUNKS[w.chunk],
    log:
        "logs/gatk/hc.{chunk}.log",
    threads: 2
    shell:
        """
        gatk --java-options "-Xmx2g" HaplotypeCaller -R {input.fa} -I {input.bam} \
          -L {params.reg} --native-pair-hmm-threads {threads} -O {output} > {log} 2>&1
        """


rule bcftools_call:
    input:
        bam=rules.align.output.bam,
        fa=REF,
    output:
        temp("results/vcf/chunks/bcftools.{chunk}.vcf.gz"),
    params:
        reg=lambda w: CHUNKS[w.chunk],
        c=config["bcftools_call"],
    log:
        "logs/bcftools/call.{chunk}.log",
    shell:
        """
        (bcftools mpileup -f {input.fa} -r {params.reg} -a FORMAT/AD,FORMAT/DP \
            -q {params.c[min_mapq]} -Q {params.c[min_baseq]} -d {params.c[max_depth]} -Ou {input.bam} \
          | bcftools call -m -v -Oz -o {output}) 2> {log}
        """


rule concat_chunks:
    input:
        lambda w: expand("results/vcf/chunks/{caller}.{chunk}.vcf.gz", caller=w.caller, chunk=CHUNKS),
    output:
        "results/vcf/{caller}.raw.vcf.gz",
    shell:
        "bcftools concat -Oz -o {output} {input} && tabix -f -p vcf {output}"


def gatk_filter_args():
    # "!vc.isSNP()" sends indels and mixed sites to the indel thresholds
    args = []
    for vtype, test in [("snp", "vc.isSNP()"), ("indel", "!vc.isSNP()")]:
        for name, expr in config["gatk_filters"][vtype].items():
            args.append(f"--filter-name {vtype.upper()}_{name} --filter-expression '{test} && {expr}'")
    return " ".join(args)


rule gatk_filter:
    input:
        vcf="results/vcf/gatk.raw.vcf.gz",
        fa=REF,
    output:
        "results/vcf/gatk.filtered.vcf.gz",
    params:
        args=gatk_filter_args(),
    log:
        "logs/gatk/filter.log",
    shell:
        """
        gatk VariantFiltration -R {input.fa} -V {input.vcf} {params.args} \
          -O results/vcf/gatk.tmp.vcf.gz > {log} 2>&1
        bcftools norm -f {input.fa} -Oz -o {output} results/vcf/gatk.tmp.vcf.gz 2>> {log}
        tabix -f -p vcf {output}; rm -f results/vcf/gatk.tmp.vcf.gz*
        """


rule bcftools_filter:
    # Soft filters: low QUAL, low/excess depth (excess depth hints at collapsed
    # repeats or paralogs) and SNVs right next to indels.
    input:
        vcf="results/vcf/bcftools.raw.vcf.gz",
        fa=REF,
        depth=rules.alignment_qc.output.summary,
    output:
        "results/vcf/bcftools.filtered.vcf.gz",
    params:
        f=config["bcftools_filters"],
    log:
        "logs/bcftools/filter.log",
    shell:
        """
        mean=$(awk '$1=="total_region" {{print $4}}' {input.depth})
        maxdp=$(awk -v m=$mean -v k={params.f[max_depth_factor]} 'BEGIN{{printf "%d", m*k}}')
        echo "mean depth $mean -> max DP $maxdp" > {log}
        bcftools filter -m+ -s LowQual -e 'QUAL<{params.f[min_qual]}' {input.vcf} \
          | bcftools filter -m+ -s LowDP -e 'INFO/DP<{params.f[min_depth]}' \
          | bcftools filter -m+ -s HighDP -e "INFO/DP>$maxdp" \
          | bcftools filter -m+ -s SnpGap -g {params.f[snp_gap]} \
          | bcftools norm -f {input.fa} -Oz -o {output} 2>> {log}
        tabix -f -p vcf {output}
        """
