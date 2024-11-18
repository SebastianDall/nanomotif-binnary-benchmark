
rule nanomotif_motif_discovery:
    conda:
        "../envs/nanomotif.yaml"
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        bins = lambda wildcards: config["samples"][wildcards.sample]["benchmarks"]["original_contig_bin"],
    output:
        o = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "bin-motifs.tsv"),
        motifs = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin","motifs.tsv"),
        m = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin","motifs-scored.tsv"),
    threads:
        40
    resources:
        mem = "40G",
        walltime = "1-00:00:00",
        nodetype = config["partition"],
    shell:
        """
        nanomotif motif_discovery \
            -t {threads} \
            --out {BASELINE}/motif_discovery/{wildcards.sample}/original_contig_bin/ \
            {input.a} {input.p} {input.bins}
        """
