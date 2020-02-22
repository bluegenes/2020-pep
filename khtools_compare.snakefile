import os
import sys
import yaml
import pandas as pd
from pep_utils import read_samples, write_yaml

sample_namelist = config.get("samplelist", "ten_haptophytes.txt")
SAMPLES = [x.strip().split('\t')[0] for x in open(sample_namelist, "r")]
#SAMPLES.remove("MMETSP0090")
#SAMPLES = ["MMETSP0090"]

info_csv = config.get("info_csv", "all_mmetsp_elvers.csv")
samplesDF = read_samples(info_csv)
out_dir = config.get("out_dir", "kcompare_out")
#data_dir = os.path.join(out_dir, "input_data")
logs_dir = os.path.join(out_dir, "logs")
wrappers_dir = "wrappers"

#pep_dir="/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"
pep_dir="/home/ntpierce/2020-pep/pep_fasta"
#pep_dir="/pylon5/mc5phkp/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"

ASSEMBLIES=["transdecoder"]

rule all:
    #input: os.path.join(out_dir, "all_by_all_kmer_comparison.parquet")
    input: os.path.join(out_dir, "busco_kmer_comparison.parquet")

def get_all_pep(w):
    #most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
    #odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    fastalist = []
    for sample in SAMPLES:
        if sample == "MMETSP0251":
            fastalist+=[os.path.join(pep_dir, f"{sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep.fasta")]
        else:
            fastalist+=[os.path.join(pep_dir, f"{sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep.fasta")]
    return fastalist

# olga guesses:
#protein, ksize 7
#hsdm17, ksize 8
#dayhoff, ksize 12
#gbmr4, 16
#hp, ksize 31

rule khtools_comparekmers:
    input: get_all_pep
    output: os.path.join(out_dir, "all_by_all_kmer_comparison.parquet")
    conda: os.path.join(wrappers_dir, "khtools-kcompare.yml")
    log: os.path.join(logs_dir, "all_by_all_kmer_comparison.log")
    benchmark: os.path.join(logs_dir, "all_by_all_kmer_comparison.benchmark")
    threads: 20
    params:
        min_k=5,
        max_k=33,
    shell:
        """
        khtools compare-kmers --processes {threads} --ksize-min {params.min_k} --ksize-max {params.max_k} --parquet {output} --intermediate-parquet {input} 2> {log}
        """

def get_buscofastas(w):
    buscofastas=[]
    for sample in SAMPLES:
        for asb in ASSEMBLIES:
            fastapath=os.path.join(out_dir, "busco_fastas", f"{sample}_{asb}", f"{sample}_{asb}" +" _{busco_id}.fasta")
            sample_fastas=expand(fastapath, busco_id=glob_wildcards(fastapath).busco_id)
            if len(sample_fastas) >= 1:
                buscofastas+=sample_fastas
    return buscofastas

rule khtools_comparebuscokmers:
    #input: expand(os.path.join(out_dir, "busco_fastas", "{sample}_{assemb}", "{sample}_{assemb}_{buscoid}.fasta"), sample=SAMPLES, assemb=["transdecoder"], busco_id=glob_wildcards(os.path.join(out_dir, "busco_fastas", "{sample}_{assemb}", "{sample}_{assemb}_{buscoid}.fasta"))
    input: get_buscofastas
    output: os.path.join(out_dir, "busco_kmer_comparison.parquet")
    conda: os.path.join(wrappers_dir, "khtools-kcompare.yml")
    log: os.path.join(logs_dir, "busco_kmer_comparison.log")
    benchmark: os.path.join(logs_dir, "busco_kmer_comparison.benchmark")
    threads: 20
    params:
        min_k=5,
        max_k=33,
    shell:
        """
        khtools compare-kmers --processes {threads} --ksize-min {params.min_k} --ksize-max {params.max_k} --parquet {output} --intermediate-parquet {input} 2> {log}
        """

