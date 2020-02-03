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

ksizes= config.get("ksizes", ["31"])
radiuses= config.get("radiuses", ["1"])

rule all:
    input: 
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}/first_doms.txt"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "plass", "{sample}_plass.fa"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "preprocess", "cutadapt", "{sample}.polyAabundtrim.fq.gz"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "preprocess", "cutadapt", "{sample}.polyAabundtrim_2.fq.gz"), sample=SAMPLES)
        #expand(os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_2.fq.gz"), sample=SAMPLES)
        expand(os.path.join(out_dir, "transdecoder", "{sample}.fasta.transdecoder.pep"), sample=SAMPLES),
        expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}_search_oh0/results.csv"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}_search_oh0/{sample}.cdbg_ids.contigs.fa.gz"), sample=["MMETSP0286"], k=ksizes, r=radiuses)
        #expand(os.path.join(sgc_configdir, "{sample}_k{k}_r{r}.yml"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "preprocess", "{sample}_1.polyAtrim.fq.gz"), sample=SAMPLES)
        expand(os.path.join(out_dir, "plass", "{sample}_plass.cdhit100.fa"), sample=SAMPLES),
        expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}/first_doms.txt"), sample=SAMPLES, k=ksizes, r=radiuses)

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

#adapter_file= "https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa"
#rule grab_adapters:
#    input: HTTP.remote(adapter_file, static=True, keep_local=True, allow_redirects=True)
#    output: os.path.join(wrappers_dir, "TruSeq3-PE.fa")
#    log: os.path.join(logs_dir, "get_data", "download_adapters.log")
#    shell: "mv {input} {output} 2> {log}"
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

rule trinity:
    input:
        left=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_1.fq.gz"),
        right=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_2.fq.gz"),
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
    threads: 10
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
    threads: 10
    conda: os.path.join(wrappers_dir, "transdecoder-env.yml") 
    wrapper: os.path.join(wrappers_dir, "transdecoder-predict.wrapper.py")

rule pear_read_merging:
    """
    Merge PE reads with PEAR, for input into PLASS
    """
    input:
        r1=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_1.fq.gz"),
        r2=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_2.fq.gz"),
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

### plass workflow = the next three rules. Going to be used both for main plass assemblies
### AND for assembling neighborhoods. Pull out --> separate subworkflow
rule plass:
    input:
        single=os.path.join(out_dir, "preprocess", "pear", '{sample}.pear_assembled.fq.gz')
        #left=os.path.join(out_dir, "khmer", "{sample}.abundtrim_1.fq.gz"),
        #right=os.path.join(out_dir, "khmer", "{sample}.abundtrim_2.fq.gz")
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
    threads: 32
    conda: os.path.join(wrappers_dir, "plass-env.yml")
    script: os.path.join(wrappers_dir, "plass-wrapper.py")

localrules: plass_remove_stop, plass_eliminate_identical_contigs

rule plass_remove_stop:
    input: os.path.join(out_dir, "plass", "{sample}_plass.fa")
    output: os.path.join(out_dir, "plass", "{sample}_plass.nostop.fa")
    log: os.path.join(logs_dir, "plass", "{sample}_rm_plass_stop.log")
    benchmark: os.path.join(logs_dir, "plass", "{sample}_rm_plass_stop.benchmark")
    conda: os.path.join(wrappers_dir, "khmer-env.yml") # should have screed in it, which is what we need    
    script: os.path.join(wrappers_dir, "remove-stop-plass.py")

rule plass_eliminate_identical_contigs:
    input: os.path.join(out_dir, "plass", "{sample}_plass.nostop.fa")
    output: os.path.join(out_dir, "plass", "{sample}_plass.cdhit100.fa")
    log: os.path.join(logs_dir, "plass", "{sample}_cdhit100.log")
    benchmark: os.path.join(logs_dir, "plass", "{sample}_cdhit100.benchmark")
    conda: os.path.join(wrappers_dir, "cdhit-env.yml")
    shell:
        """
         cd-hit -c 1 -i {input} -o  {output}
        """

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

def get_search(w):
    #most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
    #odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    if w.sample == "MMETSP0251":
        johnson_pep=os.path.join(pep_dir, f"{w.sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")
    else:
        johnson_pep=os.path.join(pep_dir, f"{w.sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")
    plass_pep=os.path.join(out_dir, "plass", f"{w.sample}_plass.fa")
    reads=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim.fq.gz")
    return [reads, johnson_pep]

localrules: write_sgc_config

rule write_sgc_config:
    input: get_search
    output: os.path.join(sgc_configdir, "{sample}_k{ksize}_r{radius}.yml"),
    run:
        configD = {"catlas_base": str(wildcards.sample), "radius": str(wildcards.radius), "ksize": str(wildcards.ksize), "search": list(map(str, input)), "input_sequences": [os.path.join(out_dir, "preprocess", "khmer", str(wildcards.sample) + ".cutpolyA.trim.abundtrim.fq.gz")]}
        with open(str(output), "w") as out:
            yaml.dump(configD, stream=out, indent=2, default_flow_style=False)

rule spacegraphcats_build:
    input:
        reads=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim.fq.gz"),
        config=os.path.join(sgc_configdir, "{sample}_k{ksize}_r{radius}.yml"),
    output:
        # this one confuses snakemake bc it doesn't have the radius wildcard
        #os.path.join(out_dir, "spacegraphcats", "{sample}/bcalm.{sample}.k{size}.unitigs.fa"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/catlas.csv"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/cdbg.gxt"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/contigs.fa.gz"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/contigs.fa.gz.indices"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/contigs.fa.gz.info.csv"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/contigs.fa.gz.mphf"),
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/first_doms.txt"),
    log: 
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_build.log")
    benchmark: 
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_build.benchmark")
    params:
        outdir = os.path.join(out_dir, "spacegraphcats")
    conda: os.path.join(wrappers_dir, "spacegraphcats-env.yml")
    shell:
        """
        python -m spacegraphcats {input.config} build --nolock --outdir={params.outdir}  
        """

rule spacegraphcats_extract_reads_contigs:
    input:
        search=get_search,
        config=os.path.join(sgc_configdir, "{sample}_k{ksize}_r{radius}.yml"),
        catlas=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}/catlas.csv")
    output:
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0/{search_file}.cdbg_ids.reads.fa.gz"), search_file = input.search),
        #os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0/{pep}.cdbg_ids.contigs.fa.gz")
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0/results.csv")
    log: 
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0.log")
    benchmark: 
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0.benchmark")
    params:
        outdir = os.path.join(out_dir, "spacegraphcats")
    conda: os.path.join(wrappers_dir, "spacegraphcats-env.yml")
    shell:
        """
        python -m spacegraphcats {input.config} extract_contigs extract_reads --nolock --outdir={params.outdir}
        """
