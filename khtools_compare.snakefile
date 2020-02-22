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

pep_dir="/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"
#pep_dir="/pylon5/mc5phkp/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"

rule all:
    input: os.path.join(out_dir, "all_by_all_kmer_comparison.parquet")

def get_all_pep(w):
    #most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
    #odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    fastalist = []
    for sample in SAMPLES:
        if sample == "MMETSP0251":
            fastalist+=[os.path.join(pep_dir, f"{sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")]
        else:
            fastalist+=[os.path.join(pep_dir, f"{sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")]
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
    threads: 20 #20
    params:
        min_k=5,
        max_k=33,
    shell:
        """
        khtools compare-kmers --processes {threads} --ksize-min {params.min_k} --ksize-max {params.max_k} --parquet {output} --intermediate-parquet {input}
        """
