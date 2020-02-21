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
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}_search_oh0/results.csv"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}_search_oh0/{sample}.cdbg_ids.contigs.fa.gz"), sample=["MMETSP0286"], k=ksizes, r=radiuses)
        #expand(os.path.join(sgc_configdir, "{sample}_k{k}_r{r}.yml"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}/first_doms.txt"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}", "{sample}_sgc_contigs_x_plass.paladin.sort.bam{end}"), sample=SAMPLES, k=ksizes, r=radiuses, end=[".bai", ".flagstat"]),


subworkflow preprocess:
    snakefile: "preprocess.snakefile"
    configfile: "dummy_config.yml"

def get_search(w):
    #most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
    #odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    if w.sample == "MMETSP0251":
        johnson_pep=os.path.join(pep_dir, f"{w.sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")
    else:
        johnson_pep=os.path.join(pep_dir, f"{w.sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")
    plass_pep=os.path.join(out_dir, "plass", f"{w.sample}_plass.fa")
    reads=os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim.fq.gz")
    return [reads, johnson_pep, plass_pep]

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
        plass_pep=rules.plass_eliminate_identical_contigs.output
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

rule assemble_nbhd:
    input:
        single=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}.cutpolyA.trim.abundtrim.fq.gz.cdbg_ids.reads.fa.gz")
    output:
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_plass.fa")
    log:
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_nbhd_plass.log")
    benchmark:
        os.path.join(logs_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_nbhd_plass.benchmark")
    params:
        #tmpdir="tmp",
        min_len=20,
        extra=" -v 3"
    #shadow: "shallow"
    threads: 10
    conda: os.path.join(wrappers_dir, "plass-env.yml")
    script: os.path.join(wrappers_dir, "plass-wrapper.py")

rule plass_nbhd_remove_stop:
    input: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_plass.fa")
    output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_plass.nostop.fa")
    log: os.path.join(logs_dir, "plass", "{sample}_k{ksize}_r{radius}_nbhd_rm_plass_stop.log")
    benchmark: os.path.join(logs_dir, "plass", "{sample}_k{ksize}_r{radius}_nbhd_rm_plass_stop.benchmark")
    conda: os.path.join(wrappers_dir, "khmer-env.yml") # should have screed in it, which is what we need
    script: os.path.join(wrappers_dir, "remove-stop-plass.py")

rule plass_nbhd_eliminate_identical_contigs:
    input: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_plass.nostop.fa")
    output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_plass.cdhit100.fa")
    log: os.path.join(logs_dir, "plass", "{sample}_k{ksize}_r{radius}_nbhd_cdhit100.log")
    benchmark: os.path.join(logs_dir, "plass", "{sample}_k{ksize}_r{radius}_nbhd_cdhit100.benchmark")
    conda: os.path.join(wrappers_dir, "cdhit-env.yml")
    shell:
        """
         cd-hit -c 1 -i {input} -o  {output}
        """

rule makeblastdb_plass_assembly:
    input: 
        os.path.join(out_dir, "plass", "{sample}_plass.cdhit100.fa")
    output:
        os.path.join(out_dir, "blast", "{sample}_plass.cdhit100.psq")
    log:
        os.path.join(logs_dir, "blast", "{sample}_plass_makeblastdb.log")
    benchmark:
        os.path.join(logs_dir, "blast", "{sample}_plass_makeblastdb.benchmark")
    params:
        #prefix = lambda wildcards: input.rsplit(".psq", 1)[0],
        prefix = lambda w: os.path.join( out_dir, "blast", w.sample +"_plass.cdhit100"),
    threads: 4
    conda: os.path.join(wrappers_dir, "blast-env.yml")
    shell:
        """
        makeblastdb -in {input} -dbtype prot -out {params.prefix} -logfile {log}
        """

rule blast_nbhd_x_plass_assembly:
    input:
        assembly=rules.plass_nbhd_eliminate_identical_contigs.output,
        #assembly=os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_plass.cdhit100.fa"),
        #db=os.path.join(out_dir, "blast", "{sample}_plass.cdhit100.psq")
        db=rules.makeblastdb_plass_assembly.output
    output:
        os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_x_assembly.out6")
    log: 
        os.path.join(logs_dir, "blast", "{sample}_k{ksize}_r{radius}_nbhd_x_assembly.log")
    benchmark:
        os.path.join(logs_dir, "blast", "{sample}_k{ksize}_r{radius}_nbhd_x_assembly.benchmark")
    conda: os.path.join(wrappers_dir, "blast-env.yml")
    params:
        prefix = lambda w: os.path.join(out_dir, "blast", w.sample + "_plass.cdhit100")
    shell:
        """
        blastp -query {input.assembly} -db {params.prefix} -num_threads {threads} -out {output} -outfmt 6 -evalue 1e-5
        """

rule compute_nbhd:
    input: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_plass.cdhit100.fa")
    output: os.path.join(out_dir, "spacegraphcats", "{sample}_k{ksize}_r{radius}_search_oh0", "{sample}_nbhd_plass.cdhit100.sig")
    params:
        k=[31],
        scaled=2000,
        compute_moltypes=["protein", "dayhoff", "hp"],
        input_is_protein=True,
        track_abundance=True,
    log: os.path.join(logs_dir, "sourmash", "{sample}_k{ksize}_r{radius}_nbhd_plass_compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_k{ksize}_r{radius}_nbhd_plass_compute.benchmark")
    conda: "sourmash-3.1.0.yml"
    script: "sourmash-compute.wrapper.py"

#include: "nhbds.snakefile"

# build all compare matrices: np and csv output
#rule sourmash_compare:
#    input: sigs=expand(os.path.join(compute_dir, "{sample}.sig"), sample= SAMPLES)
#    output:
#        np=os.path.join(compare_dir, "mmetsp_k{k}_{encoding}_compare.np"),
#        csv=os.path.join(compare_dir, "mmetsp_k{k}_{encoding}_compare.csv")
#    params:
#        include_encodings = lambda w: f"{w.encoding}",
#        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
#        k = lambda w: f"{w.k}",
#    log: os.path.join(logs_dir, "mmetsp_k{k}_{encoding}_compare.log")
#    benchmark: os.path.join(logs_dir, "mmetsp_k{k}_{encoding}_compare.benchmark")
#    conda: "sourmash-3.1.0.yml"
#    script: "sourmash-compare.wrapper.py"

