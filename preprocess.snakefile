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
SAMPLES.remove("MMETSP0090")
#SAMPLES = ["MMETSP0090"]

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
        #expand(os.path.join(out_dir, "preprocess", "cutadapt", "{sample}.polyAabundtrim.fq.gz"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "preprocess", "cutadapt", "{sample}.polyAabundtrim_2.fq.gz"), sample=SAMPLES)
        #expand(os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_2.fq.gz"), sample=SAMPLES)
        #expand(os.path.join(out_dir, "preprocess", "{sample}_1.polyAtrim.fq.gz"), sample=SAMPLES)
        #expand(os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}.codingpep.fa"), sample=SAMPLES, molecule= "protein", ksize=[5,7]),
        #expand(os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}.codingpep.fa"), sample=SAMPLES, molecule= "dayhoff", ksize=[9,11,13,15]),
       # expand(os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}.codingpep.fa"), sample=SAMPLES, molecule= "hp", ksize=[13,15,17,19,21]),
        #expand(os.path.join(out_dir, "sourmash_compare", "codingpep_k{k}_{encoding}_compare.csv"), k=["5","7","11","13","15","17","19","21"], encoding= ["protein", "dayhoff", "hp"])
        expand(os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}.codingpep.fa"), sample=SAMPLES, molecule=["hp"], ksize=["30"])
        #expand(os.path.join(out_dir, "sourmash", "codingpep_{molecule}_k{ksize}", "codingpep_k{k}_{encoding}_compare.csv"), molecule=["hp"], ksize=["30"], k=["5","7","11","13","15","17","19","21"], encoding= ["protein", "dayhoff", "hp"])

# grab the data
rule ftp_get_fq:
    output: 
        r1=os.path.join(data_dir,"{sample}_1.fq.gz"),
        r2=os.path.join(data_dir,"{sample}_2.fq.gz")
    log: os.path.join(logs_dir,"get_data", "{sample}.log")
    conda: os.path.join(wrappers_dir, "sratools-env.yml")
    params: 
        srr=lambda wildcards: samplesDF.loc[wildcards.sample]['unit'],
        out_dir=data_dir,
        compression_level="-9"
    shadow: "shallow"
    threads: 9
    shell:
        """
        fasterq-dump {params.srr} -O {params.out_dir} -e {threads} -p > {log} 2>&1
        pigz -c -p {threads} {params.compression_level} {params.out_dir}/{params.srr}_1.fastq > {output.r1}
        pigz -c -p {threads} {params.compression_level} {params.out_dir}/{params.srr}_2.fastq > {output.r2}
        rm -rf {params.out_dir}/{params.srr}_*.fastq
        """

rule polyA_trim:
    input:
        r1 = os.path.join(data_dir, "{sample}_1.fq.gz"),
        r2 = os.path.join(data_dir, "{sample}_2.fq.gz"),
    output:
        r1=os.path.join(out_dir, "preprocess", "cutadapt", "{sample}_1.cutpolyA.fq.gz"),
        r2=os.path.join(out_dir, "preprocess", "cutadapt", "{sample}_2.cutpolyA.fq.gz"),
    log: os.path.join(logs_dir, "cutadapt", "{sample}_pe.log")
    benchmark: os.path.join(logs_dir, "cutadapt", "{sample}_pe.benchmark")
    threads: 20
    conda: os.path.join(wrappers_dir, "cutadapt-env.yml")
    shell:
        """
        cutadapt -a 'A{{10}}' -A 'A{{10}}' --cores {threads} --minimum-length 31 -o {output.r1} -p {output.r2} {input.r1} {input.r2}
        """

rule adapter_trim:
    input:
        r1=os.path.join(out_dir, "preprocess", "cutadapt", "{sample}_1.cutpolyA.fq.gz"),
        r2=os.path.join(out_dir, "preprocess", "cutadapt", "{sample}_2.cutpolyA.fq.gz"),
        adapters = os.path.join(wrappers_dir, "TruSeq3-PE.fa")
    output:
        r1 = os.path.join(out_dir, "preprocess", "trimmomatic", "{sample}_1.cutpolyA.trim.fq.gz"),
        r2 = os.path.join(out_dir, "preprocess", "trimmomatic", "{sample}_2.cutpolyA.trim.fq.gz"),
        r1_unpaired = os.path.join(out_dir, "preprocess", "trimmomatic", "{sample}_1.cutpolyA.trimse.fq.gz"),
        r2_unpaired = os.path.join(out_dir, "preprocess", "trimmomatic","{sample}_2.cutpolyA.trimse.fq.gz"),
    log: os.path.join(logs_dir, "trimmomatic", "{sample}_pe.log")
    benchmark: os.path.join(logs_dir, "trimmomatic", "{sample}_pe.benchmark")
    threads: 32
    params:
        trimmer=["ILLUMINACLIP:wrappers/TruSeq3-PE.fa:2:0:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:31"],
        compression_level="-9"
    conda: os.path.join(wrappers_dir, "trimmomatic-env.yml")
    script: os.path.join(wrappers_dir, "trimmomatic-pe.py")

rule kmer_trim:
    input:
        os.path.join(out_dir, "preprocess", "trimmomatic", "{sample}_1.cutpolyA.trim.fq.gz"),
        os.path.join(out_dir, "preprocess", "trimmomatic", "{sample}_2.cutpolyA.trim.fq.gz"),
    output:
        paired=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim.fq.gz"),
        single=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim.se.fq.gz")
    log: os.path.join(logs_dir, "khmer", "{sample}_abundtrim.log")
    benchmark: os.path.join(logs_dir, "khmer", "{sample}_abundtrim.benchmark")
    params:
        k = "20",
        Z = "18",
        C = "3",
    resources:
        mem=30000000000 #30GB
        #mem=20000000000 #20GB
    conda: os.path.join(wrappers_dir, "khmer-env.yml")
    shell:
        """
        interleave-reads.py {input} | trim-low-abund.py -V -k {params.k} -Z {params.Z} -C {params.C} -M {resources.mem} - -o - | \
        (extract-paired-reads.py --gzip -p {output.paired} -s {output.single}) > {log} 2>&1
        """

rule split_pairs:
    input:
        os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim.fq.gz")
    output:
        r1=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_1.fq.gz"),
        r2=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_2.fq.gz"),
        orphans=os.path.join(out_dir, "preprocess", "cutadapt", "{sample}.cutpolyA.trim.abundtrim_orphans.fq.gz")
    conda: os.path.join(wrappers_dir, "khmer-env.yml")
    log: os.path.join(logs_dir, "khmer", "{sample}_polyAabundtrim_split_pairs.log")
    benchmark: os.path.join(logs_dir, "khmer", "{sample}_polyAabundtrim_split_pairs.benchmark")
    threads: 3
    shell:
        """
        split-paired-reads.py {input} --gzip -1 {output.r1} -2 {output.r2} -0 {output.orphans} > {log} 2>&1
        """

rule pear_read_merging:
    """
    Merge PE reads with PEAR, for input into extract_coding, trinity, etc
    """
    input:
        r1=rules.split_pairs.output.r1,
        r2=rules.split_pairs.output.r2
        #r1=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_1.fq.gz"),
        #r2=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_2.fq.gz"),
    output:
        assembled = os.path.join(out_dir, "preprocess", "pear", '{sample}.pear_assembled.fq.gz'),
        discarded = os.path.join(out_dir, "preprocess", "pear", '{sample}.pear_discarded.fq.gz'),
        unassembled_r1 = os.path.join(out_dir, "preprocess", "pear", '{sample}.pear_unassembled_r1.fq.gz'),
        unassembled_r2 = os.path.join(out_dir, "preprocess", "pear", '{sample}.pear_unassembled_r2.fq.gz'),
    message:
        """--- Merging paired reads using PEAR  ---"""
    resources:
        mem_mb=4000 # 4G
    params:
        pval = "0.01",
        extra = "",
    threads: 6
    log: os.path.join(logs_dir, "pear", "{sample}.log")
    benchmark: os.path.join(logs_dir, "pear", "{sample}.benchmark")
    conda: os.path.join(wrappers_dir, 'pear-env.yml')
    script: os.path.join(wrappers_dir, 'pear-wrapper.py')


moltypeD = {"hp": "hydrophobic-polar", "protein": "protein", "dayhoff": "dayhoff"}
rule extract_coding:
    input: rules.pear_read_merging.output.assembled
    output: 
        coding_prot=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}.codingpep.fa"),
        noncoding_nucl=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}.noncoding.fa"),
        low_complexity_nucl=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}.lowcomplexnucl.fa"),
        csv=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}.csv"),
    log: os.path.join(logs_dir, "khtools", "{sample}_{molecule}_k{ksize}.extract_coding.log")
    benchmark: os.path.join(logs_dir, "khtools", "{sample}_{molecule}_k{ksize}.extract_coding.benchmark")
    params:
        molecule=lambda w: moltypeD[w.molecule],
        ksize=lambda w: int(w.ksize),
    conda: os.path.join(wrappers_dir, "khtools-env.yml")
    shell:
        """
        khtools extract-coding --verbose --molecule {params.molecule} --peptide-ksize {params.ksize} --noncoding-nucleotide-fasta {output.noncoding_nucl} --low-complexity-nucleotide-fasta {output.low_complexity_nucl} --csv {output.csv} {input} > {output.coding_prot} 2> {log}
        """

rule sourmash_compute_reads:
    input: rules.extract_coding.output.coding_prot
    output: os.path.join(out_dir, "preprocess", "sourmash", "{sample}_{molecule}_k{ksize}.codingpep.sig")
    params:
        k=[5,7,11,13,15,17,19,21],
        scaled=2000,
        compute_moltypes=["protein", "dayhoff", "hp"],
        input_is_protein=True,
        track_abundance=True,
    log: os.path.join(logs_dir, "sourmash", "{sample}_{molecule}_k{ksize}_codingpep_compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_{molecule}_k{ksize}_codingpep_compute.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compute.wrapper.py")

# build all compare matrices: np and csv output
rule sourmash_compare_cosine: # use abundances
    input: sigs=expand(os.path.join(out_dir, "preprocess", "sourmash", "{sample}_{{molecule}}_k{{ksize}}.codingpep.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "sourmash", "codingpep_{molecule}_k{ksize}", "coding_pep_k{k}_{encoding}_cosine_compare.np"),
        csv=os.path.join(out_dir, "sourmash", "codingpep_{molecule}_k{ksize}", "codingpep_k{k}_{encoding}_cosine_compare.csv")
    params:
        include_encodings = lambda w: f"{w.encoding}",
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: f"{w.k}",
    log: os.path.join(logs_dir, "sourmash", "codingpep_{molecule}_k{ksize}_coding_pep_k{k}_{encoding}_cosine_compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "codingpep_{molecule}_k{ksize}_coding_pep_k{k}_{encoding}_cosine_compare.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule sourmash_compare_jaccard: #ignore abundances
    input: sigs=expand(os.path.join(out_dir, "preprocess", "sourmash", "{sample}_{{molecule}}_k{{ksize}}.codingpep.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "sourmash", "codingpep_{molecule}_k{ksize}", "coding_pep_k{k}_{encoding}_jaccard_compare.np"),
        csv=os.path.join(out_dir, "sourmash", "codingpep_{molecule}_k{ksize}", "codingpep_k{k}_{encoding}_jaccard_compare.csv")
    params:
        ignore_abundance = True,
        include_encodings = lambda w: f"{w.encoding}",
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: f"{w.k}",
    log: os.path.join(logs_dir, "sourmash", "codingpep_{molecule}_k{ksize}_coding_pep_k{k}_{encoding}_jaccard_compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "codingpep_{molecule}_k{ksize}_coding_pep_k{k}_{encoding}_jaccard_compare.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

#rule compare_tsne:
    #input: rules.sourmash_compare_jaccard.output.csv
#    input: rules.plass_sourmash_compare_jaccard.output.csv
    #output: os.path.join(out_dir, "tsne", "codingpep_k{k}_{encoding}_jaccard_tsne.pdf")
#    output: os.path.join(out_dir, "tsne", "codingpep_k{k}_{encoding}_jaccard_tsne.pdf")
        #report("plots/celltype-tsne.seed={seed}.pdf", caption="../report/celltype-tsne.rst", category="Dimension Reduction")
#    log: os.path.join(logs_dir, "codingpep_k{k}_{encoding}_jaccard_tsne.log")
#    benchmark: os.path.join(logs_dir, "codingpep_k{k}_{encoding}_jaccard_tsne.benchmark")
#    conda: os.path.join(wrappers_dir, "sklearn.yml")
#    shell:
#        """
#
 #       """
