# Read QC, alignment, duplicate marking and alignment QC.

rule fastp:
    input:
        r1=rules.extract_reads.output.r1,
        r2=rules.extract_reads.output.r2,
    output:
        r1=f"results/trimmed/{SAMPLE}_R1.fastq.gz",
        r2=f"results/trimmed/{SAMPLE}_R2.fastq.gz",
        json=f"results/qc/fastp/{SAMPLE}.fastp.json",
        html=f"results/qc/fastp/{SAMPLE}.fastp.html",
    log:
        "logs/fastp.log",
    threads: 4
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
          {config[fastp][extra]} -w {threads} -j {output.json} -h {output.html} 2> {log}
        """


rule bwa_index:
    input:
        REF,
    output:
        multiext(REF, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        "logs/bwa_index.log",
    shell:
        "bwa-mem2 index {input} > {log} 2>&1"


rule align:
    # bwa-mem2 -> fixmate (adds mate tags markdup needs) -> sort -> markdup
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        idx=rules.bwa_index.output,
        fa=REF,
    output:
        bam=f"results/bam/{SAMPLE}.markdup.bam",
        bai=f"results/bam/{SAMPLE}.markdup.bam.bai",
        dup=f"results/qc/{SAMPLE}.markdup.txt",
    params:
        rg=f"@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}\\tPL:ILLUMINA\\tLB:{SAMPLE}_lib1",
    log:
        "logs/align.log",
    threads: THREADS
    shell:
        """
        (bwa-mem2 mem -t {threads} -R '{params.rg}' {input.fa} {input.r1} {input.r2} \
          | samtools fixmate -m -u - - \
          | samtools sort -u -@ 2 -m 1G -T results/bam/tmp_sort - \
          | samtools markdup -@ 2 -f {output.dup} - {output.bam}) 2> {log}
        samtools index {output.bam}
        """


rule alignment_qc:
    input:
        bam=rules.align.output.bam,
        bed=rules.region_bed.output,
    output:
        flagstat=f"results/qc/{SAMPLE}.flagstat.txt",
        stats=f"results/qc/{SAMPLE}.stats.txt",
        summary=f"results/qc/mosdepth/{SAMPLE}.mosdepth.summary.txt",
        perbase=f"results/qc/mosdepth/{SAMPLE}.per-base.bed.gz",
        thresholds=f"results/qc/mosdepth/{SAMPLE}.thresholds.bed.gz",
    threads: 4
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        samtools stats -t {input.bed} {input.bam} > {output.stats}
        # per-base depth is kept: it is reused to explain FN/FP calls later
        mosdepth -t {threads} --by {input.bed} --thresholds 10,20,30 \
          results/qc/mosdepth/{SAMPLE} {input.bam}
        """


rule multiqc:
    input:
        rules.fastp.output.json,
        rules.alignment_qc.output,
        rules.align.output.dup,
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    shell:
        "multiqc -f -o results/qc -n multiqc_report.html results/qc > {log} 2>&1"
