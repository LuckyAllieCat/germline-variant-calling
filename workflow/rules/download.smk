# Download rules. Remote indexed files are queried for the region only.
# curl aborts a transfer that stalls below 1 kB/s for 60 s, then retries.

rule reference:
    retries: 2  # public servers occasionally stall or drop connections
    output:
        fa=REF,
        fai=REF + ".fai",
        dict=REF.replace(".fa", ".dict"),
    log:
        "logs/download/reference.log",
    shell:
        """
        curl -sSfL --retry 3 --retry-delay 10 --speed-limit 1000 --speed-time 60 {config[urls][reference]} | gunzip -c > {output.fa} 2> {log}
        samtools faidx {output.fa}
        samtools dict -o {output.dict} {output.fa} 2>> {log}
        """


rule region_bed:
    output:
        "results/region.bed",
    shell:
        "printf '%s\\t%d\\t%d\\n' {CHROM} $(({START} - 1)) {END} > {output}"


rule extract_reads:
    # Pull the region's reads from the remote GIAB BAM and turn them back into
    # paired FASTQ, so the pipeline re-does alignment from scratch.
    retries: 2  # public servers occasionally stall or drop connections
    output:
        r1=f"data/reads/{SAMPLE}_R1.fastq.gz",
        r2=f"data/reads/{SAMPLE}_R2.fastq.gz",
    params:
        sub=lambda w: (f"--subsample {config['downsample_fraction']} --subsample-seed 42"
                       if float(config["downsample_fraction"]) < 1 else ""),
    log:
        "logs/download/extract_reads.log",
    threads: 4
    shell:
        """
        out=$(pwd)
        # htslib saves the remote .bai in the working dir, so run inside data/reads
        cd data/reads
        # -F 0xB00: drop secondary, supplementary and QC-fail records
        samtools view -u -F 0xB00 {params.sub} {config[urls][bam]} {REGION} 2> $out/{log} \
          | samtools collate -u -O - tmp_collate \
          | samtools fastq -n -0 /dev/null -s /dev/null -c 6 \
              -1 $out/{output.r1} -2 $out/{output.r2} - 2>> $out/{log}
        """


rule truth_set:
    # GIAB benchmark VCF (region query) and high-confidence BED clipped to the region
    retries: 2  # public servers occasionally stall or drop connections
    output:
        vcf="data/truth/truth.vcf.gz",
        bed="data/truth/hiconf_region.bed",
    log:
        "logs/download/truth.log",
    shell:
        """
        tabix -h {config[urls][truth_vcf]} {REGION} 2> {log} | bgzip -c > {output.vcf}
        tabix -f -p vcf {output.vcf}
        rm -f $(basename {config[urls][truth_vcf]}).tbi
        curl -sSfL --retry 3 --retry-delay 10 --speed-limit 1000 --speed-time 60 {config[urls][truth_bed]} \
          | awk -v c={CHROM} -v s=$(({START} - 1)) -v e={END} 'BEGIN{{OFS="\\t"}}
              {{ a=$2+0; b=$3+0 }}  # force numeric comparison
              $1==c && b>s+0 && a<e+0 {{ if(a<s+0) a=s; if(b>e+0) b=e; print $1,a,b }}' \
          > {output.bed} 2>> {log}
        """


rule stratification:
    # GIAB difficult-region BEDs (no index available, so stream and keep the region)
    retries: 2  # public servers occasionally stall or drop connections
    output:
        "data/strat/{name}.bed",
    params:
        url=lambda w: config["urls"][f"strat_{w.name}"],
    shell:
        """
        curl -sSfL --retry 3 --retry-delay 10 --speed-limit 1000 --speed-time 60 {params.url} | gunzip -c \
          | awk -v c={CHROM} -v s=$(({START} - 1)) -v e={END} 'BEGIN{{OFS="\\t"}}
              {{ a=$2+0; b=$3+0 }}  # force numeric comparison
              $1==c && b>s+0 && a<e+0 {{ if(a<s+0) a=s; if(b>e+0) b=e; print $1,a,b }}' \
          > {output}
        """


rule clinvar:
    # ClinVar uses "20" not "chr20"; query the region, rename, split and normalise
    retries: 2  # public servers occasionally stall or drop connections
    input:
        fa=REF,
    output:
        "data/annotation/clinvar.vcf.gz",
    params:
        num=CHROM.removeprefix("chr"),
        reg=f"{CHROM.removeprefix('chr')}:{START}-{END}",
    log:
        "logs/download/clinvar.log",
    shell:
        """
        cd data/annotation
        url={config[urls][clinvar]}
        curl -sfI $url > /dev/null || url={config[urls][clinvar_archive]}
        echo "ClinVar source: $url" > ../../{log}
        printf '%s\\t%s\\n' {params.num} {CHROM} > chr_map.txt
        bcftools view -r {params.reg} $url 2>> ../../{log} \
          | bcftools annotate --rename-chrs chr_map.txt \
          | bcftools norm -m -any -f ../../{input.fa} -c w -Oz -o clinvar.vcf.gz 2>> ../../{log}
        tabix -f -p vcf clinvar.vcf.gz
        rm -f clinvar_*.tbi
        """


rule gene_models:
    # Ensembl GFF3 for one chromosome, with the "chr" prefix added
    retries: 2  # public servers occasionally stall or drop connections
    output:
        "data/annotation/genes.gff3.gz",
    shell:
        """
        curl -sSfL --retry 3 --retry-delay 10 --speed-limit 1000 --speed-time 60 {config[urls][gff3]} | gunzip -c \
          | awk 'BEGIN{{OFS="\\t"}} /^#/ {{print; next}} {{$1="chr"$1; print}}' \
          | gzip -c > {output}
        """
