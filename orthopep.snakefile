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
out_dir = config.get("out_dir", "orthopep_out")
data_dir = os.path.join(out_dir, "input_data")
logs_dir = os.path.join(out_dir, "logs")
wrappers_dir = "wrappers"
sgc_configdir = os.path.join(out_dir, "spacegraphcats", "sgc_config")

ksizes= config.get("ksizes", ["31"])
radiuses= config.get("radiuses", ["1"])

rule all:
    input: 
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}/first_doms.txt"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "plass", "{sample}_plass.fa"), sample=SAMPLES),
        expand(os.path.join(out_dir, "transdecoder", "{sample}.fasta.transdecoder.pep"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}_search_oh0/{sample}.cdbg_ids.contigs.fa.gz"), sample=SAMPLES, k=ksizes, r=radiuses)
        #expand(os.path.join(out_dir, "trimmed", "{sample}_1.polyAtrim.fq.gz"), sample=SAMPLES)

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
        trimmer=["ILLUMINACLIP:wrappers/TruSeq3-PE.fa:2:0:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:31"],
        compression_level="-9"
    conda: os.path.join(wrappers_dir, "trimmomatic-env.yml")
    script: os.path.join(wrappers_dir, "trimmomatic-pe.py")

rule polyA_trim:
    input: 
        r1=os.path.join(out_dir, "trimmed", "{sample}_1.trim.fq.gz"),
        r2=os.path.join(out_dir, "trimmed", "{sample}_2.trim.fq.gz"),
    output:
        r1=os.path.join(out_dir, "trimmed", "{sample}_1.polyAtrim.fq.gz"),
        r2=os.path.join(out_dir, "trimmed", "{sample}_2.polyAtrim.fq.gz"),
    log: os.path.join(logs_dir, "prinseq", "{sample}_pe.log")
    benchmark: os.path.join(logs_dir, "prinseq", "{sample}_pe.benchmark")
    threads: 10
    params:
        out_good=lambda w: os.path.join(out_dir, "trimmed","{sample}_polyAtrim"),
        out_bad="null",
        polyAlen=5,
        extra=""
    conda: os.path.join(wrappers_dir, "prinseq-env.yml")
    shell:
        """
        prinseq -fastq {input.r1} -fastq2 {input.r2} -out_good {params.out_good} -out_bad {params.out_bad} -log {log} -trim_tail_left {params.polyAlen} -trim_tail_right {params.polyAlen} -stats_all >>{log}
        mv {params.out_good}_1.f*q.gz {ouput.r1}
        mv {params.out_good}_2.f*q.gz {ouput.r2}
        """

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
        interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M {resources.mem_mb} -V - -o {output} > {log} 2>&1
        """

rule kmer_split_pairs:
    input:
        os.path.join(out_dir, "khmer", "{sample}.abundtrim.fq.gz")
    output:
        r1=os.path.join(out_dir, "khmer", "{sample}.abundtrim_1.fq.gz"),
        r2=os.path.join(out_dir, "khmer", "{sample}.abundtrim_2.fq.gz"),
        orphans=os.path.join(out_dir, "khmer", "{sample}.abundtrim_orphans.fq.gz")
    conda: os.path.join(wrappers_dir, "khmer-env.yml")
    log: os.path.join(logs_dir, "khmer", "{sample}_split_pairs.log")
    benchmark: os.path.join(logs_dir, "khmer", "{sample}_split_pairs.benchmark")
    threads: 3
    shell:
        """
        split-paired-reads.py {input} --gzip -1 {output.r1} -2 {output.r2} -0 {output.orphans} > {log} 2>&1
        """

rule trinity:
    input:
        left=os.path.join(out_dir, "khmer", "{sample}.abundtrim_1.fq.gz"),
        right=os.path.join(out_dir, "khmer", "{sample}.abundtrim_2.fq.gz")
    output: 
        fasta=os.path.join(out_dir, "trinity", "{sample}_trinity.fa"),
        gene_trans_map=os.path.join(out_dir, "trinity", "{sample}_trinity.gene_trans_map")
    message:
        """ --- Assembling read data with Trinity 2.9.1 --- """
    log: 
        os.path.join(logs_dir, "trinity", "{sample}_trinity.log")
    benchmark: 
        os.path.join(logs_dir, "trinity", "{sample}_trinity.benchmark")
    resources:
        mem_mb=100000 # 100GB
    params:
        extra=""
    shadow: "shallow"
    threads: 32
    conda: os.path.join(wrappers_dir, "trinity-env.yml")
    script: os.path.join(wrappers_dir, "trinity-wrapper.py")

rule transdecoder_longorfs:
    input:
        fasta=os.path.join(out_dir, "trinity", "{sample}_trinity.fa") 
    output:
        os.path.join(out_dir, "transdecoder", "{sample}.transdecoder_dir/longest_orfs.pep")
    log:
        os.path.join(logs_dir, "transdecoder", "{sample}.transdecoder-longorfs.log")
    benchmark: 
        os.path.join(logs_dir, "transdecoder", "{sample}.transdecoder-longorfs.benchmark")
    params:
        extra= " -m 80 "
    threads: 8
    conda: os.path.join(wrappers_dir, "transdecoder-env.yml") 
    wrapper: os.path.join(wrappers_dir, "transdecoder-longorfs.wrapper.py")
      

rule transdecoder_predict:
    input:
        fasta=os.path.join(out_dir, "trinity", "{sample}_trinity.fa"),
        longorfs=os.path.join(out_dir, "transdecoder", "{sample}.transdecoder_dir/longest_orfs.pep")
    output:
        os.path.join(out_dir, "transdecoder", "{sample}.fasta.transdecoder.bed"),
        os.path.join(out_dir, "transdecoder", "{sample}.fasta.transdecoder.cds"),
        os.path.join(out_dir, "transdecoder", "{sample}.fasta.transdecoder.pep"),
        os.path.join(out_dir, "transdecoder","{sample}.fasta.transdecoder.gff3")
    log:
        os.path.join(logs_dir, "transdecoder", "{sample}.transdecoder-predict.log")
    benchmark:
        os.path.join(logs_dir, "transdecoder", "{sample}.transdecoder-predict.benchmark")
    params:
        extra="" 
    threads: 8
    conda: os.path.join(wrappers_dir, "transdecoder-env.yml") 
    wrapper: os.path.join(wrappers_dir, "transdecoder-predict.wrapper.py")


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
    params:
        #tmpdir="tmp",
        extra=" -v 3"
    shadow: "shallow"
    threads: 28
    conda: os.path.join(wrappers_dir, "plass-env.yml")
    script: os.path.join(wrappers_dir, "plass-wrapper.py")


rule cluster_plass:
    input: 
        os.path.join(out_dir, "plass", "{sample}_plass.fa")
    output: 
        os.path.join(out_dir, "cdhit", "{sample}_plass_p{perc_id}.clstr")
    log: 
        os.path.join(logs_dir, "cdhit", "{sample}_plass_p{perc_id}.log")
    benchmark: 
        os.path.join(logs_dir, "cdhit", "{sample}_plass_p{perc_id}.benchmark")
    resources:
        mem_mb=16000 # 16GB
    threads: 10
    params:
        coverage=60,
        perc_id=lambda w: {wildcards.perc_id},
        word_size=5, # 5 for thresholds 0.7->1.0
        prefix=lambda wildcards: os.path.dirname(output),
        extra=""
    conda: os.path.join(wrappers_dir, "cdhit-env.yml")
    shell:
        """
         cd-hit -i {input} -o {params.prefix} -c {params.perc_id} -n {params.word_size}  -d 0 -aS {params.coverage}
         -aL {params.coverage} -T {threads} -M {resources.mem_mb} -o {params.prefix}  {params.extra} &> {log}
         mv {params.prefix} {output[0]} 2>> {log}
        """
         #other options
         #cd-hit -i {input} -o essMarkerGenes/marker-{wildcards.marker}_cd-hit -c {CDHIT_PERC} -d 0 -T {THREADS}
         #cd-hit-est -o cdhit -c 0.98 -i Trinity.fasta -p 1 -d 0 -b 3 -T 10

rule write_sgc_config:
    output: os.path.join(sgc_configdir, "{sample}_k{ksize}_r{radius}.yml"),
    run:
        configD = {"catlas_base": str(wildcards.sample), "radius": str(wildcards.radius), "ksize": str(wildcards.ksize), "search": [os.path.join(out_dir, "plass", str(wildcards.sample) + "_plass.fa")], "input_sequences": [os.path.join(out_dir, "khmer", str(wildcards.sample) + ".abundtrim.fq.gz")]}
        with open(str(output), "w") as out:
            yaml.dump(configD, stream=out, indent=2, default_flow_style=False)

rule spacegraphcats_build:
    input:
        reads=os.path.join(out_dir, "khmer", "{sample}.abundtrim.fq.gz"),
        config=os.path.join(sgc_configdir, "{sample}_k{ksize}_r{radius}.yml"),
    output:
 #       os.path.join(out_dir, "spacegraphcats", "{sample}/bcalm.{sample}.k{size}.unitigs.fa"),
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
