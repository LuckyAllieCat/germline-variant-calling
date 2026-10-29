# Annotation of the PASS calls: consequence (bcftools csq + Ensembl), ClinVar,
# and gnomAD frequencies (fetched remotely only for candidate sites).
# Both call sets are annotated (cheap); prioritize.py picks the better caller.

# Consequences that need a population frequency. bcftools' "~" operator does not
# support "|" alternation, so the terms are OR-ed as separate expressions.
CODING_TERMS = ["missense", "stop_gained", "stop_lost", "start_lost", "frameshift",
                "inframe", "splice_donor", "splice_acceptor"]
CANDIDATE_EXPR = " || ".join([f'INFO/BCSQ~"{t}"' for t in CODING_TERMS] + ['ID!="."'])


rule pass_split:
    # Keep PASS calls and split multi-allelic sites so ClinVar/gnomAD match per allele
    input:
        vcf="results/vcf/{caller}.filtered.vcf.gz",
        fa=REF,
    output:
        "results/annotation/{caller}.pass.vcf.gz",
    shell:
        """
        bcftools view -f PASS {input.vcf} | bcftools norm -m -any -f {input.fa} -Oz -o {output} 2> /dev/null
        tabix -f -p vcf {output}
        """


rule consequence:
    # --local-csq scores each variant on its own (like VEP), since calls are unphased
    input:
        vcf=rules.pass_split.output,
        fa=REF,
        gff=rules.gene_models.output,
    output:
        "results/annotation/{caller}.csq.vcf.gz",
    log:
        "logs/annotate/{caller}.csq.log",
    shell:
        """
        bcftools csq -f {input.fa} -g {input.gff} --local-csq --phase a -Oz -o {output} {input.vcf} 2> {log}
        tabix -f -p vcf {output}
        """


rule clinvar_annotate:
    input:
        vcf=rules.consequence.output,
        clinvar=rules.clinvar.output,
    output:
        "results/annotation/{caller}.annotated.vcf.gz",
    shell:
        """
        bcftools annotate -a {input.clinvar} -c ID,INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN \
          -Oz -o {output} {input.vcf}
        tabix -f -p vcf {output}
        """


rule candidate_sites:
    # Sites worth a population-frequency lookup: protein/splice-altering or in ClinVar
    input:
        rules.clinvar_annotate.output,
    output:
        "results/annotation/{caller}.candidate_sites.tsv",
    shell:
        """
        bcftools query -i '{CANDIDATE_EXPR}' -f '%CHROM\\t%POS\\n' {input} > {output}
        """


rule gnomad:
    # Remote queries against the gnomAD v4.1 chr20 sites VCFs (5-12 GB files, of
    # which only the blocks around the candidate sites are downloaded)
    input:
        expand("results/annotation/{caller}.candidate_sites.tsv", caller=CALLERS),
    output:
        "results/annotation/gnomad_{cohort}.tsv",
    params:
        url=lambda w: config["urls"][f"gnomad_{w.cohort}"],
    log:
        "logs/annotate/gnomad_{cohort}.log",
    threads: 4
    retries: 2  # remote object stores occasionally drop a connection
    shell:
        """
        cat {input} > {output}.sites
        bash scripts/query_gnomad.sh {output}.sites {params.url} {output} {threads} 2> {log}
        rm -f {output}.sites
        """


rule annotated_table:
    input:
        rules.clinvar_annotate.output,
    output:
        "results/annotation/{caller}.annotated.tsv",
    shell:
        """
        echo -e "CHROM\\tPOS\\tREF\\tALT\\tQUAL\\tCLINVAR_ID\\tBCSQ\\tCLNSIG\\tCLNREVSTAT\\tCLNDN\\tGT\\tDP" > {output}
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%ID\\t%INFO/BCSQ\\t%INFO/CLNSIG\\t%INFO/CLNREVSTAT\\t%INFO/CLNDN[\\t%GT\\t%DP]\\n' \
          {input} >> {output}
        """


rule prioritize:
    input:
        tables=expand("results/annotation/{caller}.annotated.tsv", caller=CALLERS),
        exomes="results/annotation/gnomad_exomes.tsv",
        genomes="results/annotation/gnomad_genomes.tsv",
        benchmark=rules.benchmark_summary.output.tsv,
        sites=rules.error_profile.output.sites,
        gff=rules.gene_models.output,
        repeats="data/strat/repeats.bed",
        lowmap="data/strat/lowmap.bed",
    output:
        table="results/tables/prioritized_variants.tsv",
        coding="results/tables/coding_variants_annotated.tsv",
        funnel="results/tables/prioritization_funnel.tsv",
        fig="figures/prioritization_summary.png",
    params:
        p=config["prioritization"],
        tables=" ".join(f"{c}:results/annotation/{c}.annotated.tsv" for c in CALLERS),
        impacts=" ".join(config["prioritization"]["impacts"]),
        clinvar=" ".join(config["prioritization"]["clinvar_keep"]),
    shell:
        """
        python scripts/prioritize.py --tables {params.tables} --caller {config[annotate_caller]} \
          --gnomad-exomes {input.exomes} --gnomad-genomes {input.genomes} \
          --benchmark {input.benchmark} --benchmark-sites {input.sites} --gff {input.gff} \
          --repeats {input.repeats} --lowmap {input.lowmap} \
          --max-af {params.p[max_af]} --impacts {params.impacts} --clinvar-keep {params.clinvar} \
          --out-table {output.table} --out-coding {output.coding} --out-funnel {output.funnel} --out-fig {output.fig}
        """


rule report:
    input:
        benchmark=rules.benchmark_summary.output.tsv,
        errors=rules.error_profile.output.tsv,
        prioritized=rules.prioritize.output.table,
        funnel=rules.prioritize.output.funnel,
        fastp=rules.fastp.output.json,
        flagstat=rules.alignment_qc.output.flagstat,
        depth=rules.alignment_qc.output.summary,
        thresholds=rules.alignment_qc.output.thresholds,
        dup=rules.align.output.dup,
    output:
        "results/tables/run_summary.md",
    shell:
        """
        python scripts/make_report.py --region {REGION} --benchmark {input.benchmark} --errors {input.errors} \
          --prioritized {input.prioritized} --funnel {input.funnel} --fastp {input.fastp} \
          --flagstat {input.flagstat} --depth {input.depth} --thresholds {input.thresholds} \
          --markdup {input.dup} --out {output}
        """
