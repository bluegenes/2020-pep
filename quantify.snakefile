



rule paladin_index:
    input: join(assembly_dir, "{assembly}.fasta")
    output: join(paladin_dir, "{assembly}.fasta.bwt"),
    params:
        reference_type= index_params.get('reference_type', '3'),
        gff = index_params.get('gff_file', '')
    log: join(logs_dir, 'paladin', "{assembly}_index.log"),
    benchmark: join(logs_dir, 'paladin', "{assembly}_index.benchmark"),
    conda: "environment.yml"
    script: 'paladin-index.py'

rule paladin_align:
    input:
        unpack(get_paladin_input),
        index = join(paladin_dir, "{assembly}.fasta.bwt"),
    output:
        join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.bam"),
    threads: 20
    params:
        f = alignment_params.get('f','125'),
        extra = alignment_params.get('extra', '')
    log: join(logs_dir, 'paladin', "{sample}_{unit}_x_{assembly}.paladin.log"),
    benchmark: join(logs_dir, 'paladin', "{sample}_{unit}_x_{assembly}.paladin.benchmark"),
    conda: "environment.yml"
    script: 'paladin-align.py'

rule samtools_sort_paladin:
    input: join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.bam")
    output: join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam")
    conda: "environment.yml"
    log: join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.log")
    benchmark: join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.benchmark")
    threads: 5
    shell:"""
    samtools sort -@ {threads} {input} -o {output}
    """

rule samtools_flagstat_paladin:
    input: join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam")
    output: join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam.flagstat")
    log: join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.flagstat.log")
    benchmark: join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.flagstat.benchmark")
    conda: "environment.yml"
    shell:"""
    samtools flagstat {input} > {output}
    """

rule samtools_index_paladin:
    input: join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam")
    output: join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam.bai")
    conda: "environment.yml"
    log: join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.index.log")
    benchmark: join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.index.benchmark")
    shell:"""
    samtools index {input}
    """
