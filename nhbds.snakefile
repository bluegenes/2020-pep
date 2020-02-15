import os
import sys
import yaml
import pandas as pd
from pep_utils import read_samples, write_yaml

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

sample_namelist = config.get("samplelist", "ten_haptophytes.txt")
SAMPLES = [x.strip().split('\t')[0] for x in open(sample_namelist, "r")]
#SAMPLES.remove("MMETSP0009")
#SAMPLES = ["MMETSP0143"]

info_csv = config.get("info_csv", "all_mmetsp_elvers.csv")
samplesDF = read_samples(info_csv)
out_dir = config.get("out_dir", "orthopep_out")
data_dir = os.path.join(out_dir, "input_data")
logs_dir = os.path.join(out_dir, "logs")
wrappers_dir = "wrappers"
sgc_configdir = os.path.join(out_dir, "spacegraphcats", "sgc_config")

pep_dir="/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"

ksizes= config.get("ksizes", ["31"])
radiuses= config.get("radiuses", ["1"])


rule all:
    input:
        expand(os.path.join(out_dir, "mmseqs2", "{sample}_sgc_k{k}_r{r}_x_plass.mmap.sam"), sample=SAMPLES, k=ksizes, r=radiuses)
        #expand(os.path.join(out_dir, "plass", "{sample}_plass.cdhit100.fa"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}/first_doms.txt"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "blast", "{sample}_plass.cdhit100.psq"),sample=SAMPLES)
#        expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}", "{sample}_sgc_contigs_x_plass.paladin.sort.bam{end}"), sample=SAMPLES, k=ksizes, r=radiuses, end=[".bai", ".flagstat"]),
        #expand(os.path.join(out_dir, "plass", "{sample}_split", "split_fasta.done"), sample=SAMPLES)

rule split_plass_fasta:
    input: os.path.join(out_dir, "plass", "{sample}_plass.fa")
    #input: rules.plass.output
    output: os.path.join(out_dir, "plass", "{sample}_split", "split_fasta.done"),
    log: os.path.join(logs_dir, "split_fasta", "{sample}_plass_split.log")
    benchmark: os.path.join(logs_dir, "split_fasta", "{sample}_plass_split.benchmark")
    conda: os.path.join(wrappers_dir, "khmer-env.yml") # uses screed
    script: os.path.join(wrappers_dir, "split-fasta.py")

# split interleaved sgc file
rule split_sgc_pairs:
    input:
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0/{sample}.cutpolyA.trim.abundtrim.fq.gz.cdbg_ids.reads.fa.gz")
    output:
        r1=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}.cdbg_ids.reads_1.fq.gz"),
        r2=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}.cdbg_ids.reads_2.fq.gz"),
        orphans=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}.cdbg_ids.reads_orphans.fq.gz"),
    conda: os.path.join(wrappers_dir, "khmer-env.yml")
    log: os.path.join(logs_dir, "khmer", "{sample}_k{ksize}_r{radius}.cdbg_ids.reads_orphans.log")
    benchmark: os.path.join(logs_dir, "khmer", "{sample}_k{ksize}_r{radius}.cdbg_ids.reads_orphans.benchmark")
    threads: 3
    shell:
        """
        split-paired-reads.py {input} --gzip -1 {output.r1} -2 {output.r2} -0 {output.orphans} > {log} 2>&1
        """

## pear merge r1, r2 files
rule pear_merge_sgc:
    """
    Merge PE reads with PEAR, for input into PLASS
    """
    input:
        r1=rules.split_sgc_pairs.output.r1,
        r2=rules.split_sgc_pairs.output.r2,
    output:
        assembled = os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", '{sample}.cdbg_ids.pear_assembled.fq.gz'),
        discarded = os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", '{sample}.cdbg_ids.pear_discarded.fq.gz'),
        unassembled_r1 = os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", '{sample}.cdbg_ids.pear_unassembled_r1.fq.gz'),
        unassembled_r2 = os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", '{sample}.cdbg_ids.pear_unassembled_r2.fq.gz'),
    message:
        """--- Merging paired reads using PEAR  ---"""
    resources:
        mem_mb=4000 # 4G
    params:
        pval = "0.01",
        extra = "",
    threads: 6
    log: os.path.join(logs_dir, "pear", "{sample}_k{ksize}_r{radius}.sgc_cdbg_ids.log")
    benchmark: os.path.join(logs_dir, "pear", "{sample}_k{ksize}_r{radius}.sgc_cdbg_ids.benchmark")
    conda: os.path.join(wrappers_dir, 'pear-env.yml')
    script: os.path.join(wrappers_dir, 'pear-wrapper.py')


#paladin map reads --> sample plass assembly
rule plass_paladin_index:
    input: os.path.join(out_dir, "plass", "{sample}_plass.fa") 
    #input: rules.plass.output 
    output: os.path.join(out_dir, "plass", "{sample}_plass.fasta.bwt"),
    params:
        reference_type="3",
        gff = "", 
    log: os.path.join(logs_dir, 'paladin', "{sample}_plass_index.log"),
    benchmark: os.path.join(logs_dir, 'paladin', "{sample}_plass_index.benchmark"),
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    script: os.path.join(wrappers_dir, 'paladin-index.py')

rule paladin_map_sgc_contigs_x_plass:
    input:
        index=rules.plass_paladin_index.output,
        r=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "contigs.fa.gz"),
    output:
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "{sample}_sgc_contigs_x_plass.paladin.bam")
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    threads: 20
    log: os.path.join(logs_dir, "paladin", "{sample}_k{ksize}_r{radius}_sgc_contigs_x_plass.log")
    benchmark:os.path.join(logs_dir, "paladin", "{sample}_k{ksize}_r{radius}_sgc_contigs_x_plass.benchmark")
    params:
        f = "125", 
        extra = "",
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    script: os.path.join(wrappers_dir, 'paladin-align.py')

rule paladin_map_search_reads_x_plass_assembly:
    input:
        index=rules.plass_paladin_index.output,
        r=rules.pear_merge_sgc.output.assembled,
    output:
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_x_plass.paladin.bam")
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    threads: 20
    log: os.path.join(logs_dir, "paladin", "{sample}_k{ksize}_r{radius}_nbhd_x_assembly.log")
    benchmark:os.path.join(logs_dir, "paladin", "{sample}_k{ksize}_r{radius}_nbhd_x_assembly.benchmark")
    params:
        f = "125", 
        extra = "",
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    script: os.path.join(wrappers_dir, 'paladin-align.py')


rule mmseqs_createdb_plass:
   input: os.path.join(out_dir, "plass", "{sample}_plass.fa")
   output: os.path.join(out_dir, "plass", "{sample}_plass.mmseqsDB") 
   log: os.path.join(logs_dir, "mmseqs2", "{sample}_plass.createDB.log")
   benchmark: os.path.join(logs_dir, "mmseqs2", "{sample}_plass.createDB.benchmark")
   conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
   shell:
       """
       mmseqs createdb {input} {output}
       """

rule mmseqs_createdb_sgc_contigs:
   input: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "contigs.fa.gz") 
   output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "contigs.mmseqsDB") 
   log: os.path.join(logs_dir, "mmseqs", "{sample}_k{ksize}_r{radius}_contigs.createDB.log")
   benchmark: os.path.join(logs_dir, "mmseqs", "{sample}_k{ksize}_r{radius}_contigs.createDB.benchmark")
   conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
   shell:
       """
       mmseqs createdb {input} {output} tmp
       """

rule mmseqs_map:
# https://github.com/soedinglab/MMseqs2/wiki#mapping-very-similar-sequences-using-mmseqs-map
# also cluster w/ linclust: mmseqs linclust inDB outDB tmp OR rbh
    input: 
        plass=os.path.join(out_dir, "plass", "{sample}_plass.mmseqsDB"),
        sgc=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "contigs.mmseqsDB")
    output: os.path.join(out_dir, "mmseqs2","{sample}_sgc_k{ksize}_r{radius}_x_plass.mmap.resultsDB"),
    log: os.path.join(logs_dir, "mmseqs2", "{sample}_sgc_k{ksize}_r{radius}_x_plass.mmap.log")
    benchmark: os.path.join(logs_dir, "mmseqs2", "{sample}_sgc_k{ksize}_r{radius}_x_plass.mmap.benchmark")
    conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
    shell:
        """
        mmseqs map {input.sgc} {input.plass} resultDB tmp
        """

rule mmseqs_results_to_sam:
    input: 
        target_db=os.path.join(out_dir, "plass", "{sample}_plass.mmseqsDB"),
        query_db=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "contigs.mmseqsDB"),
        results_db=rules.mmseqs_map.output
    output: os.path.join(out_dir, "mmseqs2", "{sample}_sgc_k{ksize}_r{radius}_x_plass.mmap.sam"),
    log: os.path.join(logs_dir, "mmseqs2", "{sample}_sgc_k{ksize}_r{radius}_x_plass.convert.log")
    benchmark: os.path.join(logs_dir, "mmseqs2", "{sample}_sgc_k{ksize}_r{radius}_x_plass.convert.benchmark")
    conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
    threads: 2
    shell:
        """
        mmseqs convertalis --threads {threads} --format-mode 1 {input.query_db} {input.target_db} {input.results_db} {output}
        """
# mmseqs extractorfs YourCrassphageDB crassphageOrfs --min-length
# http://wwwuser.gwdg.de/~compbiol/molbio_course/2018/worksheet-day2.pdf

# make a samtools env instead of putting it in paladin (will use for mmseqs2 output as well)
rule samtools_sort_paladin:
    input: rules.paladin_map_sgc_contigs_x_plass.output
    #input: rules.paladin_map_nbhdreads_x_plass_assembly.output 
    output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "{sample}_sgc_contigs_x_plass.paladin.sort.bam")
    #output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_x_plass.paladin.sort.bam")
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    log: os.path.join(logs_dir, 'paladin', "{sample}_k{ksize}_r{radius}_sgc_contigs_x_plass.samtools_sort.log")
    benchmark: os.path.join(logs_dir, 'paladin', "{sample}_k{ksize}_r{radius}_sgc_contigs_x_assembly.samtools_sort.benchmark") 
    threads: 5
    shell:"""
    samtools sort -@ {threads} {input} -o {output}
    """

rule samtools_flagstat_paladin:
    input: rules.samtools_sort_paladin.output 
    output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "{sample}_sgc_contigs_x_plass.paladin.sort.bam.flagstat") 
    #output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_x_plass.paladin.sort.bam.flagstat")
    log: os.path.join(logs_dir, 'paladin', "{sample}_k{ksize}_r{radius}_nbhd_x_assembly.samtools_flagstat.log") 
    benchmark: os.path.join(logs_dir, 'paladin', "{sample}_k{ksize}_r{radius}_nbhd_x_assembly.samtools_flagstat.benchmark") 
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    shell:"""
    samtools flagstat {input} > {output}
    """

rule samtools_index_paladin:
    input: rules.samtools_sort_paladin.output
    output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "{sample}_sgc_contigs_x_plass.paladin.sort.bam.bai") 
    #output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_x_plass.paladin.sort.bam.bai") 
    conda: os.path.join(wrappers_dir, "paladin-env.yml")
    log: os.path.join(logs_dir, 'paladin', "{sample}_k{ksize}_r{radius}_nbhd_x_assembly.samtools_index.log") 
    benchmark: os.path.join(logs_dir, 'paladin', "{sample}_k{ksize}_r{radius}_nbhd_x_assembly.samtools_index.benchmark") 
    shell:"""
    samtools index {input}
    """

#rule paladin_quant:
#    input:
#        assembly="outputs/cd-hit95/{nbhd}.cdhit95.faa",
#        bam=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}", "{sample}_sgc_contigs_x_plass.paladin.sort.bam")
#    output: os.path.join(out_dir, "salmon", "{nhbd}_quant/quant.sf")
#    params:
#        out="outputs/salmon/{nbhd}_quant"
#    conda: "salmon-env.yml"
#    shell:'''
#    salmon quant -t {input.cdhit} -l A -a {input.bam} -o {params.out}
#    '''

## do something like this?
# https://github.com/taylorreiter/hu-snake/blob/6dadc5f4f230f89a93f4855523cb213213462b12/PLASS/variant_snakemake/Snakefile
#rule extract_plass_readnames_from_bam:
# extract names of reads for each aa seq bam section that matched 
# PFAM domain of interest
#    output: dynamic("outputs/bam_subsets/{plass_match}-NAMES.txt") 
#    input:
#        plass_names = expand("outputs/hmmscan/{nbhd}.{trim}trim.plass.{pfam_base}-hmmscanT100-NAMES.txt", nbhd = NBHD, trim = TRIM, pfam_base=PFAM_BASE), 
#        bai = expand("outputs/paladin/{nbhd}.{trim}trim.sort.bam.bai", nbhd = NBHD, trim = TRIM), 
#        bam = expand("outputs/paladin/{nbhd}.{trim}trim.sort.bam", nbhd = NBHD, trim = TRIM) 
#    run:
#        import pysam
#        import re
#        
#        with open(str(input.plass_names)) as f:
#            plass_hmmscan_matches = f.readlines()
#
#        plass_hmmscan_matches = [x.strip() for x in plass_hmmscan_matches] 
#
#        samfile = pysam.AlignmentFile(str(input.bam), 'rb')
##
#        for match in plass_hmmscan_matches:
#            reads = []
#            for read in samfile.fetch(match):
#                name = re.sub('.*:', '', read.qname)
#                reads.append(name)
#                with open(f"outputs/bam_subsets/{match}-NAMES.txt", 'w') as outfile:
#                    for s in reads:
#                        outfile.write("%s\n" % s)
#        samfile.close()

#blastp -query {input.assembly} -db {params.prefix} -num_threads {threads} -out {output} -outfmt 6 -evalue 1e-5

