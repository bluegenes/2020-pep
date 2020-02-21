import os


## nucleotide-level quantification: salmon
rule salmon_index:
    """
    Index the transcriptome for salmon quantification
    """
    input:
        os.path.join(assembly_dir, "{assembly}.fasta") 
    output:
        directory(os.path.join(quant_dir, "{assembly}.salmonindex"))
    message:
        """--- Indexing the transcriptome with Salmon ---"""
    threads: 10
    params:
        extra = salmon_params['index_params'].get('extra', '')
    log:os.path.join(logs_dir, 'salmon','{assembly}_index.log')
    benchmark:os.path.join(logs_dir, 'salmon','{assembly}_index.benchmark')
    conda: os.path.join(wrappers_dir, "salmon-env.yml")
	script: os.path.join(wrappers_dir, "salmon-index.wrapper.py")

rule salmon_quant:
    """
    Quantify transcripts with Salmon
    """
    input: ## add reads
        r1=
        r2=
        #unpack(get_sample_no_combine),
        index = os.path.join(quant_dir, "{assembly}.salmonindex") 
    output:
        quant = os.path.join(quant_dir,"{sample}_{unit}_x_{assembly}.salmon", "quant.sf"),
        lib = os.path.join(quant_dir, "{sample}_{unit}_x_{assembly}.salmon", "lib_format_counts.json")
    message:
        """--- Quantifying transcripts with Salmon ---"""
    params:
        libtype = salmon_params['quant_params'].get('libtype', 'A'),
        extra = salmon_params['quant_params'].get('extra', '')
    threads: 20
    log:os.path.join(logs_dir, 'salmon/{sample}_{unit}_x_{assembly}.log')
    benchmark:os.path.join(logs_dir, 'salmon/{sample}_{unit}_x_{assembly}.benchmark')
    conda: os.path.join(wrappers_dir, "salmon-env.yml")
    script: os.path.join(wrappers_dir, "salmon-quant.wrapper.py")


## protein-level quantification
rule paladin_index:
    input: os.path.join(assembly_dir, "{assembly}.fasta")
    output: os.path.join(paladin_dir, "{assembly}.fasta.bwt"),
    params:
        reference_type= index_params.get('reference_type', '3'),
        gff = index_params.get('gff_file', '')
    log: os.path.join(logs_dir, 'paladin', "{assembly}_index.log"),
    benchmark: os.path.join(logs_dir, 'paladin', "{assembly}_index.benchmark"),
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    script: os.path.join(wrappers_dir, 'paladin-index.py')

rule paladin_align:
    input:
        unpack(get_paladin_input),
        index = os.path.join(paladin_dir, "{assembly}.fasta.bwt"),
    output:
        os.path.join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.bam"),
    threads: 20
    params:
        f = alignment_params.get('f','125'),
        extra = alignment_params.get('extra', '')
    log: os.path.join(logs_dir, 'paladin', "{sample}_{unit}_x_{assembly}.paladin.log"),
    benchmark: os.path.join(logs_dir, 'paladin', "{sample}_{unit}_x_{assembly}.paladin.benchmark"),
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    script: os.path.join(wrappers_dir, 'paladin-align.py')

# make a samtools env instead of putting it in paladin (will use for mmseqs2 output as well)
rule samtools_sort_paladin:
    input: os.path.join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.bam")
    output: os.path.join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam")
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    log: os.path.join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.log")
    benchmark: os.path.join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.benchmark")
    threads: 5
    shell:"""
    samtools sort -@ {threads} {input} -o {output}
    """

rule samtools_flagstat_paladin:
    input: os.path.join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam")
    output: os.path.join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam.flagstat")
    log: os.path.join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.flagstat.log")
    benchmark: os.path.join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.flagstat.benchmark")
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    shell:"""
    samtools flagstat {input} > {output}
    """

rule samtools_index_paladin:
    input: os.path.join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam")
    output: os.path.join(paladin_dir,"{sample}_{unit}_x_{assembly}.paladin.sort.bam.bai")
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    log: os.path.join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.index.log")
    benchmark: os.path.join(logs_dir, 'paladin',"{sample}_{unit}_x_{assembly}.paladin.sort.bam.index.benchmark")
    shell:"""
    samtools index {input}
    """

rule salmon_quant_paladin:
    input:
        cdhit="outputs/cd-hit95/{nbhd}.cdhit95.faa",
        bam=rules.paladin_align.output
    output: "outputs/salmon/{nbhd}_quant/quant.sf"
    params:
        out="outputs/salmon/{nbhd}_quant"
    conda: ENV
    shell:'''
    salmon quant -t {input.cdhit} -l A -a {input.bam} -o {params.out}
    '''


rule mmseqs2:
# https://github.com/soedinglab/MMseqs2/wiki#mapping-very-similar-sequences-using-mmseqs-map
# also cluster w/ linclust: mmseqs linclust inDB outDB tmp OR rbh
   input: os.path.join(out_dir, "{assembly}_plass.fa")
   output: os.path.join(out_dir, "mmseqs2","{sample}_x_{assembly}.plass.bam"),
   conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
   log: os.path.join(logs_dir, "mmseqs2", "{sample}_x_{assembly}.plass.log")
   benchmark: os.path.join(logs_dir, "mmseqs2", "{sample}_x_{assembly}.plass.benchmark")
   shell:
       """
       mmseqs map queryDB targetDB resultDB tmp 
       """

