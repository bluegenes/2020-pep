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
#SAMPLES.remove("MMETSP0090")
#SAMPLES = ["MMETSP0090"]

info_csv = config.get("info_csv", "all_mmetsp_elvers.csv")
samplesDF = read_samples(info_csv)
out_dir = config.get("out_dir", "orthopep_out")
data_dir = os.path.join(out_dir, "input_data")
logs_dir = os.path.join(out_dir, "logs")
wrappers_dir = "wrappers"
sgc_configdir = os.path.join(out_dir, "spacegraphcats", "sgc_config")

pep_dir="/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"


rule all:
    input: 
        #expand(os.path.join(out_dir, "busco", "{sample}_{assemb}_busco_euk"), sample=SAMPLES, assemb=["trinity", "plass", "jpep", "transdecoder"]),
        os.path.join(out_dir, "busco", "busco_summaries", "busco_figure.png")

def get_jpep(w):
    #most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
    #odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    if w.sample == "MMETSP0251":
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")
    else:
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")

rule busco_trans:
    input: os.path.join(out_dir, "trinity", "{sample}_trinity.fa")
    output: directory(os.path.join(out_dir, "busco", "{sample}_trinity_busco_euk"))
    log: os.path.join(logs_dir, "busco", "{sample}_trinity_busco.log")
    benchmark: os.path.join(logs_dir, "busco","{sample}_trinity_busco.benchmark")
    threads: 8
    params:
        mode="transcriptome",
        lineage="eukaryota_odb10",
        #auto_lineage="euk",
    conda: os.path.join(wrappers_dir, "busco-env.yml")
    script: os.path.join(wrappers_dir, "busco-wrapper.py")

rule busco_pep:
    input: os.path.join(out_dir, "transdecoder", "{sample}_trinity.transdecoder.pep")
    output: directory(os.path.join(out_dir, "busco", "{sample}_transdecoder_busco_euk"))
    log: os.path.join(logs_dir, "busco", "{sample}_transdecoder_busco.log")
    benchmark: os.path.join(logs_dir, "busco","{sample}_transdecoder_busco.benchmark")
    threads: 8
    params:
        mode="protein",
        lineage="eukaryota_odb10",
        #auto_lineage="euk",
    conda: os.path.join(wrappers_dir, "busco-env.yml")
    script: os.path.join(wrappers_dir, "busco-wrapper.py")

rule busco_plass:
    input: os.path.join(out_dir, "plass", "{sample}_plass.fa")
    output: directory(os.path.join(out_dir, "busco", "{sample}_plass_busco_euk"))
    log: os.path.join(logs_dir, "busco","{sample}_plass_busco.log")
    benchmark: os.path.join(logs_dir, "busco","{sample}_plass_busco.benchmark")
    threads: 8
    params:
        mode="protein",
        lineage="eukaryota_odb10",
        #auto_lineage="euk",
        extra=""
    conda: os.path.join(wrappers_dir, "busco-env.yml")
    script: os.path.join(wrappers_dir, "busco-wrapper.py")

rule busco_jpep:
    input: get_jpep 
    output: directory(os.path.join(out_dir, "busco", "{sample}_jpep_busco_euk"))
    log: os.path.join(logs_dir, "busco","{sample}_jpep_busco.log")
    benchmark: os.path.join(logs_dir, "busco","{sample}_jpep_busco.benchmark")
    threads: 8
    params:
        mode="protein",
        lineage="eukaryota_odb10",
        #auto_lineage="euk",
        extra=""
    conda: os.path.join(wrappers_dir, "busco-env.yml")
    script: os.path.join(wrappers_dir, "busco-wrapper.py")


localrules: plot_busco_summaries

def get_short_summary(w):
    sums = []
    for sample in SAMPLES:
        for assemb in ["trinity", "plass", "jpep"]:
            short_sum = os.path.join(out_dir, "busco", "{sample}_{assemb}_busco", "short_summary.specific.eukaryota_odb10.{sample}_{assemb}_busco.txt")
            if not os.path.exists(short_sum):
                short_sum = os.path.join(out_dir, "busco", "{sample}_{assemb}_busco", "short_summary.generic.eukaryota_odb10.{sample}_{assemb}_busco.txt")
               # short_sum = os.path.join(out_dir, "busco", "{sample}_{assemb}_busco", "auto_lineage", "run_eukaryota_odb10", "short_summary.txt")
            if os.path.exists(short_sum):
                sums+=[short_sum]
    return sums

rule plot_busco_summaries:
    input: expand(os.path.join(out_dir, "busco", "{sample}_{assemb}_busco_euk", "short_summary.specific.eukaryota_odb10.{sample}_{assemb}_busco_euk.txt"), sample=SAMPLES, assemb=["trinity", "plass", "jpep", "transdecoder"])
    #input: 
        #summary_files=get_short_summary
    output: os.path.join(out_dir, "busco", "busco_summaries", "busco_figure.png")
    conda: os.path.join(wrappers_dir, "busco-env.yml")
    log: os.path.join(logs_dir, "busco", "plot_summaries.log")
    params:
        outdir= os.path.join(out_dir, "busco", "busco_summaries")
    shell:
        """
        mkdir -p {params.outdir}
        cp {input} {params.outdir}
        generate_plot.py --working_directory {params.outdir} 2> {log}
        """
