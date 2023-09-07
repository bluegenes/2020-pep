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

# for all mmetsp samples
info_csv = config.get("info_csv", "all_mmetsp_elvers.csv")
samplesDF = read_samples(info_csv)
MMETSP_SAMPLES = samplesDF["sample"].tolist()

out_dir = config.get("out_dir", "orthopep_out")
data_dir = os.path.join(out_dir, "input_data")
logs_dir = os.path.join(out_dir, "logs")
wrappers_dir = "wrappers"
sgc_configdir = os.path.join(out_dir, "spacegraphcats", "sgc_config")

pep_dir="/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"

ksizes= config.get("ksizes", ["31"])
radiuses= config.get("radiuses", ["1"])

moltypeD = {"hp": "hydrophobic-polar", "protein": "protein", "dayhoff": "dayhoff"}
pepRefD = {"sprot": "/home/ntpierce/2020-pep/khtools_testing/uniprot_sprot.fasta.gz", "merc": "/home/ntpierce/2020-pep/khtools_testing/MERC.fasta.gz"}



rule all:
    input:
        expand(os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pepref}.codingpep.diamond_nr.out"), sample=SAMPLES, molecule= "dayhoff", ksize=[9,11,13,15],                        pepref=["sprot"] ,merc"]),
        expand(os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pepref}.codingpep.diamond_nr.out"), sample=SAMPLES, molecule= "protein", ksize=[5,7], pepref=["sprot", "merc"]),   #        expand(os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pepref}.codingpep.diamond_nr.out"), sample=SAMPLES, molecule= "hp", ksize=[13,15,17,19,21], pepref=["sprot",       "merc"])


rule sencha_index:
    input:
        pep_ref= lambda w: pepRefD[w.pep_ref]
    output:
        index=os.path.join(out_dir, "preprocess", "khtools", "reference", "ref{pep_ref}.codingpep.fa"),
    log: os.path.join(logs_dir, "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.extract_coding.log")
    benchmark: os.path.join(logs_dir, "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.extract_coding.benchmark")
    params:
        molecule=lambda w: moltypeD[w.molecule],
    wildcard_constraints:
        ksize=["\d+"],
    conda: os.path.join(wrappers_dir, "khtools-env.yml")
    shell:
        #khtools extract-coding --verbose --molecule {params.molecule} --peptide-ksize {wildcards.ksize} --noncoding-nucleotide-fasta {output.noncoding_nucl} --low-complexity-nucleotide-fasta                     {output.             low_complexity_nucl} --csv {output.csv} {input.pep_ref} {input.reads} > {output.coding_prot} 2> {log}
        """
        khtools index --alphabet {wildcards.alphabet} --peptide-ksize {wildcards.ksize} > {output} 2> {log}
        """

rule sencha_translate:
    input:
        reads=rules.pear_read_merging.output.assembled,
        pep_ref= lambda w: pepRefD[w.pep_ref]
    output:
        coding_prot=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.codingpep.fa"),
        noncoding_nucl=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.noncoding.fa"),
        low_complexity_nucl=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.lowcomplexnucl.fa"),
        csv=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.csv"),
    log: os.path.join(logs_dir, "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.extract_coding.log")
    benchmark: os.path.join(logs_dir, "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.extract_coding.benchmark")
    params:
        molecule=lambda w: moltypeD[w.molecule],
    #wildcard_constraints:
    #    ksize=["\d+"],
    conda: os.path.join(wrappers_dir, "khtools-env.yml")
    shell:
        """
        khtools extract-coding --verbose --molecule {params.molecule} --peptide-ksize {wildcards.ksize} --noncoding-nucleotide-fasta {output.noncoding_nucl} --low-complexity-nucleotide-fasta {output.             low_complexity_nucl} --csv {output.csv} {input.pep_ref} {input.reads} > {output.coding_prot} 2> {log}
        """

rule diamond_makedb_nr:
    input: "/group/ctbrowngrp/pierce/databases/nr.gz"
    output: "/group/ctbrowngrp/pierce/databases/nr.dmnd"
    conda: os.path.join(wrappers_dir, "diamond-env.yml")
    log: os.path.join(logs_dir, "diamond", "nr_makedb.log")
    benchmark: os.path.join(logs_dir, "diamond", "nr_makedb.benchmark")
    shell:
        """
        diamond makedb --in {input} --db {output} 2> {log}
        """
# for orthofinder, output needs to be of form: Blast{species}_{species}.txt.gz
rule diamond_blastp_extract_coding:
    input:
        pep = rules.extract_coding.output.coding_prot,
        db = rules.diamond_makedb_nr.output
    output: os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.codingpep.diamond_nr.out")
    log: os.path.join(logs_dir, "diamond", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.codingpep.diamond_nr.log")
    benchmark: os.path.join(logs_dir, "diamond", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.codingpep.diamond_nr.benchmark")
    shadow: "shallow"
    conda: os.path.join(wrappers_dir, "diamond-env.yml")
    shell:
        """
        diamond blastp -d {input.db} -q {input.pep} -o {output} 2> {log}
        """
