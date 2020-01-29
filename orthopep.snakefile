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

info_csv = config.get("info_csv", "all_mmetsp_elvers.csv")
samplesDF = read_samples(info_csv)
# subset by samples --> not necessary, bc just generating samples from rule all
#samplesDF[samplesDF['sample'].isin(samples)]
out_dir = config.get("out_dir", "orthopep_out")
data_dir = os.path.join(out_dir, "input_data")
logs_dir = os.path.join(out_dir, "logs")
wrappers_dir = "wrappers"
sgc_configdir = "sgc_config"

ksizes= config.get("ksizes", ["31"])
radiuses= config.get("radiuses", ["1"])

rule all:
    input: 
        expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}/first_doms.txt"), sample=SAMPLES, k=ksizes, r=radiuses),
        expand(os.path.join(out_dir, "plass", "{sample}_plass.fa"), sample=SAMPLES),

# grab the data
rule ftp_get_fq1:
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

adapter_file= "https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa"
rule grab_adapters:
    input: HTTP.remote(adapter_file, static=True, keep_local=True, allow_redirects=True)
    output: os.path.join(wrappers_dir, "TruSeq3-PE.fa")
    log: os.path.join(logs_dir, "get_data", "download_adapters.log")
    shell: "mv {input} {output} 2> {log}"

rule adapter_trim:
    input:
        r1 = os.path.join(data_dir, "{sample}_1.fq.gz"),
        r2 = os.path.join(data_dir, "{sample}_2.fq.gz"),
        adapters = os.path.join(wrappers_dir, "TruSeq3-PE.fa")
    output:
        r1 = os.path.join(out_dir, "trimmed", "{sample}_1.trim.fq.gz"),
        r2 = os.path.join(out_dir, "trimmed", "{sample}_2.trim.fq.gz"),
        r1_unpaired = os.path.join(out_dir, "trimmed", "{sample}_1.trimse.fq.gz"),
        r2_unpaired = os.path.join(out_dir, "trimmed", "{sample}_2.trimse.fq.gz"),
    log: os.path.join(logs_dir, "trimmomatic", "{sample}_pe.log")
    benchmark: os.path.join(logs_dir, "trimmomatic", "{sample}_pe.benchmark")
    threads: 32
    params:
        trimmer=["ILLUMINACLIP:{snakemake.input.adapters}:2:0:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:31"],
        compression_level="-9"
    conda: os.path.join(wrappers_dir, "trimmomatic-env.yml")
    script: os.path.join(wrappers_dir, "trimmomatic-pe.py")

rule kmer_trim:
    input: 
        os.path.join(out_dir, "trimmed", "{sample}_1.trim.fq.gz"),
        os.path.join(out_dir, "trimmed", "{sample}_2.trim.fq.gz"),
    output: 
        paired=os.path.join(out_dir, "khmer", "{sample}.abundtrim.fq.gz")
    log: os.path.join(logs_dir, "khmer", "{sample}_abundtrim.log")
    benchmark: os.path.join(logs_dir, "khmer", "{sample}_abundtrim.benchmark")
    resources:
        mem_mb=6000000000 #6GB
    conda: os.path.join(wrappers_dir, "khmer-env.yml")
    shell:
        """
        interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M {resources.mem_mb} -V - -o {output} 2> {log}
        """

rule kmer_split_pairs:
    input:
        os.path.join(out_dir, "khmer", "{sample}.abundtrim.fq.gz")
    output:
        r1=os.path.join(out_dir, "khmer", "{sample}.abundtrim_1.fq.gz"),
        r2=os.path.join(out_dir, "khmer", "{sample}.abundtrim_2.fq.gz")
    conda: os.path.join(wrappers_dir, "khmer-env.yml")
    threads: 3
    shell:
        """
        split-paired-reads.py {input} --gzip -1 {output.r1} -2 {output.r2} > {log} 2>&1
        """

rule plass:
    input:
        left=os.path.join(out_dir, "khmer", "{sample}.abundtrim_1.fq.gz"),
        right=os.path.join(out_dir, "khmer", "{sample}.abundtrim_2.fq.gz")
    output: 
        os.path.join(out_dir, "plass", "{sample}_plass.fa")
    log: 
        os.path.join(logs_dir, "plass", "{sample}_plass.log")
    benchmark: 
        os.path.join(logs_dir, "plass", "{sample}_plass.benchmark")
    threads: 16
    conda: os.path.join(wrappers_dir, "plass-env.yml")
    script: os.path.join(wrappers_dir, "plass-wrapper.py")


rule write_sgc_config:
    output: os.path.join(sgc_configdir, "{sample}_k{ksize}_r{radius}.yml"),
    run:
        configD = {"catlas_base": str(wildcards.sample), "radius": str(wildcards.radius), "ksize": str(wildcards.ksize), "search": [os.path.join(out_dir, "plass", str(wildcards.sample) + "_plass.fa")]}
        with open(str(output), "w") as out:
            yaml.dump(configD, stream=out, indent=2, default_flow_style=False)

rule spacegraphcats_build:
    input:
        reads=os.path.join(out_dir, "khmer", "{sample}.abundtrim.fq.gz"),
        config=os.path.join(sgc_configdir, "{sample}_k{ksize}_r{radius}.yml"),
    output:
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/bcalm.{sample}_k{ksize}_r{radius}.k31.unitigs.fa"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/catlas.csv"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/cdbg.gxt"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/contigs.fa.gz"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/contigs.fa.gz.indices"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/contigs.fa.gz.info.csv"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/contigs.fa.gz.mphf"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/first_doms.txt"),
    log: 
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0.log")
    benchmark: 
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0.benchmark")
    params:
        outdir = os.path.join(out_dir, "spacegraphcats")
    conda: os.path.join(wrappers_dir, "spacegraphcats-env.yml")
    shell:
        """
        python -m spacegraphcats {input.config} build --nolock --outdir={params.outdir}  
        """

rule spacegraphcats_extract_reads_contigs:
    input:
        reads=os.path.join(out_dir, "khmer", "{sample}.abundtrim.fq.gz"),
        config=os.path.join(sgc_configdir, "{sample}_k{ksize}_r{radius}.yml"),
    output:
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0/{sample}.cdbg_ids.reads.fa.gz"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0/{sample}.cdbg_ids.contigs.fa.gz")
    log: 
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0.log")
    benchmark: 
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0.benchmark")
    params:
        outdir = os.path.join(out_dir, "spacegraphcats")
    conda: os.path.join(wrappers_dir, "spacegraphcats-env.yml")
    shell:
        """
        python -m spacegraphcats {input.config} build extract_contigs extract_reads --nolock --outdir={params.outdir}  
        """
