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

#ASSEMBLIES=["trinity", "jpep", "transdecoder"] #plass
ASSEMBLIES=["transdecoder"] #plass
rule all:
    input: 
        #expand(os.path.join(out_dir, "busco", "{sample}_{assemb}_busco_euk"), sample=SAMPLES, assemb=["trinity", "plass", "jpep", "transdecoder"]),
        #os.path.join(out_dir, "busco", "busco_summaries", "busco_figure.png"),
        #directory(os.path.join(out_dir, "busco_fastas"))
        expand(os.path.join(out_dir, "busco_sigs", "buscohashes_above_{min_count}_scaled{scaled}_{encoding}_k{k}.sig"), scaled=10, encoding="dayhoff", k=31, min_count=2),
        expand(os.path.join(out_dir, "busco_filtsigs_plots", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np.matrix.pdf"), scaled=10, encoding="dayhoff", k=31, min_count=2)

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

# start with single copy
# orthopep_out/busco/MMETSP0091_transdecoder_busco_euk/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences
checkpoint copy_busco_fastas:
    input: os.path.join(out_dir, "busco", "{sample}_{assemb}_busco_euk", "short_summary.specific.eukaryota_odb10.{sample}_{assemb}_busco_euk.txt")
    output: directory(os.path.join(out_dir, "busco_fastas", "{sample}_{assemb}")) #, "{sample}_{assemb}_{buscoid}.fasta")
    #conda: os.path.join(wrappers_dir, "khmer-env.yml")
    log: os.path.join(logs_dir, "copy_busco_fastas", "{sample}_{assemb}.log")
    params:
        basename = lambda w: w.sample + "_" + w.assemb,
        source_dir= lambda w: os.path.join(out_dir, "busco", f"{w.sample}_{w.assemb}_busco_euk", "run_eukaryota_odb10", "busco_sequences", "single_copy_busco_sequences"),
        dest_dir = os.path.join(out_dir, "busco_fastas", "{sample}_{assemb}")
    script: os.path.join(wrappers_dir, "copy_busco_fastas.py") 


rule sourmash_compute:
    input: fasta=os.path.join(out_dir, "busco_fastas", "{sample}_{assemb}", "{sample}_{assemb}_{buscoid}.fasta")
    output: os.path.join(out_dir, "busco_sigs", "{sample}_{assemb}", "{sample}_{assemb}_{buscoid}_scaled{scaled}_{encoding}_k{k}.sig")
    params:
        k= lambda w: w.k,
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: w.encoding,
        input_is_protein=True,
        track_abundance=True,
        #singleton=True,
    log: os.path.join(logs_dir, "sourmash", "{sample}_{assemb}_{buscoid}_scaled{scaled}_{encoding}_k{k}_compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_{assemb}_{buscoid}_scaled{scaled}_{encoding}_k{k}_compute.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compute.wrapper.py")  

def aggregate_sigs(w):
    """
    aggregate the file names of the unknown number of files generated at the copy_busco_fastas step
    """
    sigfiles=[]
    for sp in SAMPLES:
        for asb in ASSEMBLIES:
            checkpoint_output=checkpoints.copy_busco_fastas.get(sample=sp, assemb=asb).output[0]
            sigpath= os.path.join(out_dir, "busco_sigs", f"{sp}_{asb}", f"{sp}_{asb}" +"_{buscoid}_scaled{scaled}_{encoding}_k{k}.sig")
            sigfiles+=expand(sigpath, sample=sp, assemb=asb, buscoid=glob_wildcards(os.path.join(checkpoint_output, f"{sp}_{asb}" + "_{busco_id}.fasta")).busco_id, scaled=w.scaled, encoding=w.encoding, k=w.k)
    return sigfiles


# maybe do this by assembly? Likely same unique/err hashes across my four diff "assemblies." Or maybe also try mincount=5?
rule drop_unique_hashes:
    input: aggregate_sigs
    output:
        hashes=os.path.join(out_dir, "busco_sigs", "buscohashes_above_{min_count}_scaled{scaled}_{encoding}_k{k}.hashes.txt"),
        sig=os.path.join(out_dir, "busco_sigs", "buscohashes_above_{min_count}_scaled{scaled}_{encoding}_k{k}.sig"),
    log: os.path.join(logs_dir, "sourmash", "buscohashes_above_{min_count}_scaled{scaled}_{encoding}_k{k}.log")
    benchmark: os.path.join(logs_dir, "sourmash", "buscohashes_above_{min_count}_scaled{scaled}_{encoding}_k{k}.benchmark")
    params:
        scaled= lambda w: w.scaled,
        ksize= lambda w: w.k,
        molecule= lambda w: w.encoding,
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "drop_unique_hashes.py")

#checkpoint intersect_to_drop_unique:
rule intersect_to_drop_unique:
    input:
        keep_hashes=os.path.join(out_dir, "busco_sigs", "buscohashes_above_{min_count}_scaled{scaled}_{encoding}_k{k}.sig"),
        sig=os.path.join(out_dir, "busco_sigs", "{sample}_{assemb}", "{sample}_{assemb}_{buscoid}_scaled{scaled}_{encoding}_k{k}.sig")
    output:
        filt=os.path.join(out_dir, "busco_filtsigs", "{sample}_{assemb}_{buscoid}_scaled{scaled}_{encoding}_k{k}_above{min_count}.sig"),
        filt_renamed=os.path.join(out_dir, "busco_filtsigs", "{sample}_{assemb}_{buscoid}_scaled{scaled}_{encoding}_k{k}_above{min_count}_renamed.sig")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash signature intersect -o {output.filt} --abundances-from {input.sig} -k {wildcards.k} {input.sig} {input.keep_hashes}
        sourmash signature rename -o {output.filt_renamed} -k {wildcards.k} {input.sig} {wildcards.sample}_above{wildcards.min_count}_renamed
        """

def aggregate_filt_sigs(w):
    """
    aggregate the file names of the unknown number of files generated at the split_fasta step
    """
    sigfiles=[]
    for sp in SAMPLES:
        for asb in ASSEMBLIES:
            checkpoint_output=checkpoints.copy_busco_fastas.get(sample=sp, assemb=asb).output[0]
            sigpath= os.path.join(out_dir, "busco_filtsigs", f"{sp}_{asb}" +"_{buscoid}_scaled{scaled}_{encoding}_k{k}_above{min_count}_renamed.sig")
            sigfiles+=expand(sigpath, sample=sp, assemb=asb, buscoid=glob_wildcards(os.path.join(checkpoint_output, f"{sp}_{asb}" + "_{buscoid}.fasta")).buscoid, scaled=w.scaled, encoding=w.encoding, k=w.k, min_count=w.min_count)
    return sigfiles

rule sourmash_compare_cosine_nounique:
    input: 
        sigs=aggregate_filt_sigs, 
    output:
        np=os.path.join(out_dir, "busco_filtsigs_compare", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.np"),
        csv=os.path.join(out_dir, "busco_filtsigs_compare", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.csv"),
    params:
        include_encodings = lambda w: w.encoding,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.k,
    log: os.path.join(logs_dir, "sourmash", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.log")
    benchmark: os.path.join(logs_dir, "sourmash", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule sourmash_compare_jaccard_nounique:
    input: 
        sigs=aggregate_filt_sigs,
    output:
        np=os.path.join(out_dir, "busco_filtsigs_compare", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np"),
        csv=os.path.join(out_dir, "busco_filtsigs_compare", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.csv"),
    params:
        ignore_abundance=True,
        include_encodings = lambda w: w.encoding,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.k,
    log: os.path.join(logs_dir, "sourmash", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.log")
    benchmark: os.path.join(logs_dir, "sourmash", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

# sourmash plot each compare matrix numpy output
rule sourmash_plot_cosine_nounique:
    input: os.path.join(out_dir, "busco_filtsigs_compare", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.np"),
    output: os.path.join(out_dir, "busco_filtsigs_plots", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.np.matrix.pdf"),
    params:
        plot_dir=os.path.join(out_dir, "busco_filtsigs_plots")
    log: os.path.join(logs_dir, "sourmash", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_cosine_plot.log")
    benchmark: os.path.join(logs_dir, "sourmash", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_cosine_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """

rule sourmash_plot_jaccard_nounique:
    input: os.path.join(out_dir, "busco_filtsigs_compare", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np"),
    output: os.path.join(out_dir, "busco_filtsigs_plots", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np.matrix.pdf"),
    params:
        plot_dir=os.path.join(out_dir, "busco_filtsigs_plots")
    log: os.path.join(logs_dir, "sourmash", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_jaccard_plot.log")
    benchmark: os.path.join(logs_dir, "sourmash", "buscohashes_scaled{scaled}_{encoding}_k{k}_above{min_count}_jaccard_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """
