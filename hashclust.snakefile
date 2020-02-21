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

ksizes= config.get("ksizes", [5,7,11,13,15,17,19,21,25,31,35,41])
radiuses= config.get("radiuses", ["1"])

rule all:
    input:
        #expand(os.path.join(out_dir, "preprocess", "cutadapt", "{sample}.polyAabundtrim.fq.gz"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "hashclust", "all_hashes_above_{min_count}_scaled{scaled}_k{ksize}_mol{moltype}.txt"), scaled=["200", "500", "1000","2000"], min_count=["2"], ksize=[5,7,11,13,15,17,19,21], moltype=["protein", "dayhoff", "hp"])
        #expand(os.path.join(out_dir, "hashclust", "all_hashes_above_{min_count}_plass100_scaled{scaled}_{moltype}_k{k}.sig"), scaled=["1"], min_count=["2"], moltype="dayhoff", k=["31"]), #,"50", "200", "500", "1000","2000"], min_count=["2"])
        #expand(os.path.join(out_dir, "hashclust", "sigs", "{sample}_plass100_scaled{scaled}_{moltype}_k{k}_mincount{min_count}.sig"), sample=SAMPLES, scaled=[1], min_count=[2], moltype="dayhoff", k=[31])
        #expand(os.path.join(out_dir, "hashclust", "sigs", "scaled{scaled}_{moltype}_k{k}", "{sample}_plass100_above{min_count}_renamed.sig"), sample=SAMPLES, scaled=[1], min_count=[2], moltype="dayhoff", k=[31])
        expand(os.path.join(out_dir, "plass", "plots", "plass_cdhit100_scaled{scaled}_{moltype}_k{k}_above{min_count}_{ctype}_compare.np.matrix.pdf"), scaled=[1], min_count=[2], moltype="dayhoff",k=ksizes, ctype=["cosine", "jaccard"])

### try dropping unique hashes before clustering. 
#First, build table of good remaining hashes

## would be best to do this with reads, but let's start with plass assemblies
rule drop_unique_hashes_plass:
# rule modified from @taylorreiter get_greater_than_1_filt_sigs at https://github.com/taylorreiter/ibd/blob/master/Snakefile
    input: expand(os.path.join(out_dir, "plass", "sigs", "{sample}_plass_cdhit100_scaled{{scaled}}_{{moltype}}_k{{ksize}}.sig"), sample=SAMPLES) 
    output: 
        hashes=os.path.join(out_dir, "hashclust", "all_hashes_above_{min_count}_plass100_scaled{scaled}_{moltype}_k{ksize}.hashes.txt"),
        sig=os.path.join(out_dir, "hashclust", "all_hashes_above_{min_count}_plass100_scaled{scaled}_{moltype}_k{ksize}.sig"),
    log: os.path.join(logs_dir, "hashclust", "allhashes_above_{min_count}_plass100_scaled{scaled}_{moltype}_k{ksize}.log")
    benchmark: os.path.join(logs_dir, "hashclust", "allhashes_above_{min_count}_plass100_scaled{scaled}_{moltype}_k{ksize}.benchmark")
    params:
        scaled= lambda w: w.scaled,
        ksize= lambda w: w.ksize,
        molecule= lambda w: w.moltype,
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "drop_unique_hashes.py")

#rule drop_unique_hashes_coding_reads:
# rule modified from @taylorreiter get_greater_than_1_filt_sigs at https://github.com/taylorreiter/ibd/blob/master/Snakefile
    #input: expand(os.path.join(out_dir, "plass", "{sample}_plass_cdhit100_scaled{{scaled}}.sig"), sample=SAMPLES) 
    #output: 
        #hashes=os.path.join(out_dir, "hashclust", "all_hashes_above_{min_count}_scaled{scaled}_k{ksize}_mol{moltype}.txt"),
    #    sig=os.path.join(out_dir, "hashclust", "all_hashes_above_{min_count}_scaled{scaled}.sig"),
    #log: os.path.join(logs_dir, "hashclust", "allhashes_plass_cdhit100_scaled{scaled}_drop_counts_under{min_count}.log")
    #benchmark: os.path.join(logs_dir, "hashclust", "allhashes_plass_cdhit100_scaled{scaled}_drop_counts_under{min_count}.benchmark")
    #params:
    #    scaled= lambda w: w.scaled,
    #    ksize= lambda w: w.ksize,
    #    molecule= lambda w: w.moltype,
    #conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    #script: os.path.join(wrappers_dir, "drop_unique_hashes.py")


rule intersect_to_drop_unique:
    input: 
        keep_hashes=os.path.join(out_dir, "hashclust", "all_hashes_above_{min_count}_plass100_scaled{scaled}_{moltype}_k{ksize}.sig"),
        sig=os.path.join(out_dir, "plass", "sigs", "{sample}_plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}.sig")
    output: 
        filt= os.path.join(out_dir, "hashclust", "sigs", "scaled{scaled}_{moltype}_k{ksize}", "{sample}_plass100_above{min_count}.sig"),
        filt_renamed=os.path.join(out_dir, "hashclust", "sigs", "{sample}_plass100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_renamed.sig") 
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash signature intersect -o {output.filt} -A {input.sig} -k {wildcards.ksize} {input.sig} {input.keep_hashes}
        sourmash signature rename -o {output.filt_renamed} -k {wildcards.ksize} {input.sig} {wildcards.sample}_above{wildcards.min_count}_renamed
        """

# build all compare matrices: np and csv output
rule sourmash_compare_cosine: # use abundances
    input: sigs=expand(os.path.join(out_dir, "hashclust", "sigs", "{sample}_plass100_scaled{{scaled}}_{{moltype}}_k{{ksize}}_above{{min_count}}_renamed.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "plass", "compare", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.np"),
        csv=os.path.join(out_dir, "plass", "compare", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.csv"),
    params:
        include_encodings = lambda w: w.moltype,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.ksize,
    log: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.log") 
    benchmark: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.benchmark") 
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule sourmash_compare_jaccard: #ignore abundances
    input: sigs=expand(os.path.join(out_dir, "hashclust", "sigs", "{sample}_plass100_scaled{{scaled}}_{{moltype}}_k{{ksize}}_above{{min_count}}_renamed.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "plass", "compare", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.np"),
        csv=os.path.join(out_dir, "plass", "compare", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.csv"),
    params:
        ignore_abundance = True,
        include_encodings = lambda w: w.moltype,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.ksize,
    log: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.log") 
    benchmark: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.benchmark") 
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule sourmash_plot_cosine:
    input: rules.sourmash_compare_cosine.output.np
    output: os.path.join(out_dir, "plass", "plots", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_compare.np.matrix.pdf"),
    params:
        plot_dir=os.path.join(out_dir, "plass", "plots")
    log: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_plot.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_cosine_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input} 2> {log}
        """

rule sourmash_plot_jaccard:
    input: rules.sourmash_compare_jaccard.output.np
    output: os.path.join(out_dir, "plass", "plots", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_compare.np.matrix.pdf"),
    params:
        plot_dir=os.path.join(out_dir, "plass", "plots")
    log: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_plot.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass_cdhit100_scaled{scaled}_{moltype}_k{ksize}_above{min_count}_jaccard_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input} 2> {log}
        """


#rule compare_tsne:
#    input:
#        sce="analysis/normalized.batch-removed.rds",
#        fits=expand("analysis/cellassign.{parent}.rds", parent=markers["parent"].unique())
#    output:
#        report("plots/celltype-tsne.seed={seed}.pdf",
#                   caption="../report/celltype-tsne.rst",
#                   category="Dimension Reduction")
#    log:
#        "logs/celltype-tsne/seed={seed}.log"
#    conda:
#        "../envs/eval.yaml"
#    wildcard_constraints:
#        seed="[0-9]+"
#    script:
#        "../scripts/celltype-tsne.R"
