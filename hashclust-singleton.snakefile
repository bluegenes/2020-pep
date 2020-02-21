"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  hashclust-singleton.snakefile --use-conda -n
"""

import os
import re
from pep_utils import read_samples


configfile: "config_mmetsp_pepsmash_singleton.yml"

# read in elvers csv, grab sample names
samples_csv = config.get("samples", None)
if not samples_csv:
    sys.stderr.write(f"\n\tError: Please provide samples file via yml config file\n\n")
    sys.exit(-1)

samplesDF = read_samples(samples_csv)
SAMPLES = samplesDF["sample"].tolist()
SAMPLES.remove("MMETSP0754") # pep file is empty!!!
SAMPLES = SAMPLES[0:5]
#print(SAMPLES)
pep_dir = config["pep_dir"]
scaled_vals=config.get("scaled", [1])
ksizes=config.get("ksizes", ["7", "11", "17"])
encodings=config.get("encodings", ["protein", "hp", "dayhoff"])

#build some output dirs
out_dir = config.get("out_dir", "mmetsp_pepsmash")
logs_dir = os.path.join(out_dir, "logs")
#compute_dir = os.path.join(out_dir, "johnson_pep", "singleton_sigs")
compute_dir = os.path.join(out_dir, "johnson_pep", "splitpep_sigs")
compare_dir = os.path.join(out_dir, "johnson_pep", "singleton_compare")
plots_dir = os.path.join(out_dir, "johnson_pep", "singleton_plots")
wrappers_dir = config.get("wrappers_dir", "wrappers") 

compare_exts= [".np", ".csv"]

#sigfiles = expand(os.path.join(compute_dir, "{sample}.sig"), sample= SAMPLES)
#comparefiles = expand(os.path.join(compare_dir, "mmetsp_k{k}_{encoding}_compare{end}"), k=ksizes, encoding=encodings, end=compare_exts)
#plotfiles = expand(os.path.join(plots_dir, "mmetsp_k{k}_{encoding}_compare_{type}.np.matrix.pdf"), k=ksizes, encoding=encodings, type=["cosine", "jaccard"], end=compare_exts)

rule all:
    input: 
#        expand(os.path.join(compute_dir, "all_hashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.hashes.txt"), sample=SAMPLES, k=ksizes, scaled=scaled_vals,encoding=encodings, min_count=2)
#        expand(os.path.join(compute_dir, "{sample}_jpep_scaled{scaled}_{encoding}_k{k}.sig"), sample=SAMPLES, k=ksizes, scaled=scaled_vals,encoding=encodings)
        expand(os.path.join(compute_dir, "all_hashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.sig"), k=ksizes, scaled=scaled_vals, encoding=encodings, min_count=2),
        #expand(os.path.join(plots_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np.matrix.pdf"), k=ksizes, scaled=scaled_vals, encoding=encodings, min_count=2)

def get_pep(w):
    #most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
    #odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    if w.sample == "MMETSP0251":
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")
    else:
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")

## can either split the fasta or split the singleton sig. let's do the fasta, then compute each contig separately

# this checkpoint creates an {tnum} variable (splits fasta to produce {sample}_{tnum}.fasta)
checkpoint split_pep:
    input: get_pep 
    output: 
        directory(os.path.join(out_dir, "johnson_pep", "split-fastas", "{sample}_split")),
        #os.path.join(out_dir, "johnson_pep", "split-fastas", "{sample}_split", "split_fasta.done"),
    params:
        prefix=lambda w: w.sample,
    log: os.path.join(logs_dir, "split_fasta", "{sample}_jpep_split.log")
    benchmark: os.path.join(logs_dir, "split_fasta", "{sample}_jpep_split.benchmark")
    conda: os.path.join(wrappers_dir, "khmer-env.yml") # requires screed
    script: os.path.join(wrappers_dir, "split-fasta.py")


#def get_fasta(w):
#    checkpoint_output=checkpoints.split_pep.get(**w).output[0]
#    fastapath=os.path.join(out_dir, "johnson_pep", "split-fastas", "{sample}_split", "{sample}_{tnum}.fasta")
#    sp = w.sample
#    return(expand(fastapath, tnum=glob_wildcards(os.path.join(checkpoint_output, f"{sp}" + "_{tnum}.fasta")).tnum, sample=w.sample))
        #import pdb;pdb.set_trace()


rule sourmash_compute:
    input: fasta=os.path.join(out_dir, "johnson_pep", "split-fastas", "{sample}_split", "{sample}_{tnum}.fa")
         #split_done=os.path.join(out_dir, "johnson_pep", "split-fastas", "{sample}_split", "split_fasta.done"),
    output: os.path.join(compute_dir, "{sample}_split", "{sample}_{tnum}_jpep_scaled{scaled}_{encoding}_k{k}.sig")
    params:
        k= lambda w: w.k,
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: w.encoding,
        input_is_protein=True,
        track_abundance=True,
        #singleton=True,
    log: os.path.join(logs_dir, "sourmash", "{sample}_{tnum}_jpep_scaled{scaled}_{encoding}_k{k}_compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_{tnum}_jpep_scaled{scaled}_{encoding}_k{k}_compute.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compute.wrapper.py")

def aggregate_sigs(w):
    """
    aggregate the file names of the unknown number of files generated at the split_fasta step
    """
    sigfiles=[]
    for sp in SAMPLES:
        checkpoint_output=checkpoints.split_pep.get(sample=sp).output[0] #(**wildcards) 
        sigpath= os.path.join(compute_dir, f"{sp}_split", f"{sp}_{{tnum}}_jpep_scaled{{scaled}}_{{encoding}}_k{{k}}.sig")
        #sigfiles+=expand(os.path.join(compute_dir, "{sample}_split", "{sample}_{tnum}_jpep_scaled{scaled}_{encoding}_k{k}.sig"), sample=sp, tnum=glob_wildcards(os.path.join(checkpoint_output, f"{sp}_{{tnum}}.fasta")).tnum, scaled=w.scaled, encoding=w.encoding, k=w.k)
        #tnum=glob_wildcards(os.path.join(checkpoint_output, f"{sp}" + "_{tnum}.fasta")).tnum
        #print(tnum)
        #import pdb;pdb.set_trace()
        sigfiles+=expand(sigpath, sample=sp, tnum=glob_wildcards(os.path.join(checkpoint_output, f"{sp}" + "_{tnum}.fa")).tnum, scaled=w.scaled, encoding=w.encoding, k=w.k)
        #print(sigfiles)
        #import pdb;pdb.set_trace()
    return sigfiles 

rule drop_unique_hashes_jpep:
# rule modified from @taylorreiter get_greater_than_1_filt_sigs at https://github.com/taylorreiter/ibd/blob/master/Snakefile
    input: aggregate_sigs
    #input: expand(os.path.join(compute_dir, "{sample}_split", "{sample}_{{tnum}}_jpep_scaled{{scaled}}_{{encoding}}_k{{k}}.sig"), sample=SAMPLES)
    output:
        hashes=os.path.join(compute_dir, "all_hashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.hashes.txt"),
        sig=os.path.join(compute_dir, "all_hashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.sig"),
    log: os.path.join(logs_dir, "sourmash", "allhashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.log")
    benchmark: os.path.join(logs_dir, "sourmash", "allhashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.benchmark")
    params:
        scaled= lambda w: w.scaled,
        ksize= lambda w: w.k,
        molecule= lambda w: w.encoding,
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "drop_unique_hashes.py")

# how does intersect work with --singleton?? It doesn't.
rule intersect_to_drop_unique:
    input:
        keep_hashes=os.path.join(out_dir, "hashclust", "all_hashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.sig"),
        sig=os.path.join(compute_dir, "{sample}_split", "{sample}_{tnum}_jpep_scaled{scaled}_{encoding}_k{k}.sig")
    output:
        filt=os.path.join(compute_dir, "{sample}_split_filt", "{sample}_{tnum}_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}.sig"),
        filt_renamed=os.path.join(compute_dir, "{sample}_split_filt", "{sample}_{tnum}_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_renamed.sig")
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
    #sigpath= os.path.join(compute_dir, f"{sp}_split", f"{sp}_{{tnum}}_jpep_scaled{{scaled}}_{{encoding}}_k{{k}}_above{{min_count}}_renamed.sig")
    for sp in SAMPLES:
        checkpoint_output=checkpoints.split_pep.get(sample=sp).output[0] #(**wildcards)
        sigfiles+=expand(os.path.join(compute_dir, f"{sp}_split_filt", "{sample}_{tnum}_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_renamed.sig"), sample=sp, tnum=glob_wildcards(os.path.join(checkpoint_output,  f"{sp}" +"_{tnum}.fa")).tnum, scaled=w.scaled, encoding=w.encoding, k=w.k, min_count=w.min_count)
        #import pdb;pdb.set_trace()
    return sigfiles

rule sourmash_compare_cosine_nounique:
    input: sigs=aggregate_filt_sigs
    #input: sigs= expand(os.path.join(compute_dir, "{sample}_split", "{sample}_{tnum}_jpep_scaled{{scaled}}_{{encoding}}_k{{k}}_above{{min_count}}_renamed.sig"), sample=SAMPLES)
    output: 
        np=os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.np"),
        csv=os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.csv"),
    params:
        include_encodings = lambda w: w.encoding,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.k,
    log: os.path.join(logs_dir, "compare", "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.log")
    benchmark: os.path.join(logs_dir, "compare", "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule sourmash_compare_jaccard_nounique:
    input: sigs=aggregate_filt_sigs
    #input: sigs= expand(os.path.join(compute_dir, "{sample}_jpep_scaled{{scaled}}_{{encoding}}_k{{k}}_above{{min_count}}_renamed.sig"), sample=SAMPLES)
    output: 
        np=os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np"),
        csv=os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.csv"),
    params:
        ignore_abundance=True,
        include_encodings = lambda w: w.encoding,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.k,
    log: os.path.join(logs_dir, "compare", "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.log")
    benchmark: os.path.join(logs_dir, "compare", "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

# sourmash plot each compare matrix numpy output
rule sourmash_plot_cosine_nounique:
    input: os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.np"), 
    output: os.path.join(plots_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.np.matrix.pdf"),
    params:
        plot_dir=plots_dir 
    log: os.path.join(logs_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_cosine_plot.log")
    benchmark: os.path.join(logs_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_cosine_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """

rule sourmash_plot_jaccard_nounique:
    input: os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np"), 
    output: os.path.join(plots_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np.matrix.pdf"),
    params:
        plot_dir=plots_dir 
    log: os.path.join(logs_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_jaccard_plot.log")
    benchmark: os.path.join(logs_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_jaccard_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """
