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

pep_dir="/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"

rule all:
    input: 
        #expand(os.path.join(out_dir, "plass", "{sample}_plass.fa"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "transdecoder", "{sample}_trinity.transdecoder.pep"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "preprocess", "{sample}_1.polyAtrim.fq.gz"), sample=SAMPLES)
        #expand(os.path.join(out_dir, "plass", "{sample}_plass.cdhit100.fa"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "blast", "{sample}_plass.cdhit100.psq"),sample=SAMPLES)
        #expand(os.path.join(out_dir, "plass", "{sample}_split", "split_fasta.done"), sample=SAMPLES)
        #expand(os.path.join(out_dir, "cdhit", "{sample}_plass_cdhit{perc_id}.fa"), sample=SAMPLES, perc_id=[".95", ".97", ".85", ".80", ".70"])
        #expand(os.path.join(out_dir, "plass", "sigs", "{sample}_plass_cdhit100_scaled{scaled}_{encoding}_k{k}.sig"), sample=SAMPLES, scaled= ["1", "50", "200", "500", "1000","2000"], k=[5,7,11,13,15,17,19,21,25,31,35,41,45], encoding=["protein", "dayhoff", "hp"]),
        #expand(os.path.join(out_dir, "plass", "sourmash", "plots", "plass100_scaled{scaled}_{encoding}_k{k}_{type}_compare.np.matrix.pdf"),sample=SAMPLES, scaled= ["1", "50", "200", "500", "1000","2000"], k=[5,7,11,13,15,17,19,21,25,31,35,41], encoding=["protein", "dayhoff", "hp"], type=["jaccard", "cosine"])
        directory(os.path.join(out_dir, "sonicparanoid", "tenhapto"))


subworkflow preprocess:
    snakefile: "preprocess.snakefile"
    configfile: "dummy_config.yml"


rule trinity:
    input:
        left=preprocess(os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_1.fq.gz")),
        right=preprocess(os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_2.fq.gz")),
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
        fasta=rules.trinity.output.fasta
        #fasta=os.path.join(out_dir, "trinity", "{sample}_trinity.fa") 
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
    script: os.path.join(wrappers_dir, "transdecoder-longorfs.wrapper.py")
      
rule transdecoder_predict:
    input:
        fasta=os.path.join(out_dir, "trinity", "{sample}_trinity.fa"),
        longorfs=os.path.join(out_dir, "transdecoder", "{sample}.transdecoder_dir/longest_orfs.pep")
    output:
        bed=os.path.join(out_dir, "transdecoder", "{sample}_trinity.transdecoder.bed"),
        cds=os.path.join(out_dir, "transdecoder", "{sample}_trinity.transdecoder.cds"),
        pep=os.path.join(out_dir, "transdecoder", "{sample}_trinity.transdecoder.pep"),
        gff=os.path.join(out_dir, "transdecoder","{sample}_trinity.transdecoder.gff3")
    log:
        os.path.join(logs_dir, "transdecoder", "{sample}_trinity.transdecoder-predict.log")
    benchmark:
        os.path.join(logs_dir, "transdecoder", "{sample}_trinity.transdecoder-predict.benchmark")
    params:
        extra="" 
    threads: 10
    conda: os.path.join(wrappers_dir, "transdecoder-env.yml") 
    script: os.path.join(wrappers_dir, "transdecoder-predict.wrapper.py")


rule copy_pepfiles:
    input: rules.transdecoder_predict.output.pep
    output: os.path.join(out_dir, "sonicparanoid", "tenhapto_pepfiles", "{sample}_trinity.transdecoder.pep"),
    log: os.path.join(logs_dir, "utils", "copy_{sample}_trinity.transdecoder.pep.log")
    shell:
        """
        cp {input} {output} 2> {log}
        """

rule sonicparanoid:
    input: expand(os.path.join(out_dir, "sonicparanoid", "tenhapto_pepfiles", "{sample}_trinity.transdecoder.pep"), sample=SAMPLES)
    output: directory(os.path.join(out_dir, "sonicparanoid", "tenhapto"))
    params:
        name="test_ten_hapto",
        indir=os.path.join(out_dir, "sonicparanoid", "tenhapto_pepfiles")
    conda: os.path.join(wrappers_dir, "sonicparanoid-env.yml")
    shell:
        """
        sonicparanoid -i {params.indir} -o {output} -m fast -t 30 --project-id {params.name}
        """

rule pear_read_merging:
    """
    Merge PE reads with PEAR, for input into PLASS
    """
    input:
        r1=preprocess(os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_1.fq.gz")),
        r2=preprocess(os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_2.fq.gz")),
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
        left=preprocess(os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_1.fq.gz")),
        right=preprocess(os.path.join(out_dir, "preprocess", "khmer", "{sample}.cutpolyA.trim.abundtrim_2.fq.gz")),
    output: 
        os.path.join(out_dir, "plass", "{sample}_plass.fa")
    log: 
        os.path.join(logs_dir, "plass", "{sample}_plass.log")
    benchmark: 
        os.path.join(logs_dir, "plass", "{sample}_plass.benchmark")
    params:
        #tmpdir="tmp",
        min_len=30,
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
    input: os.path.join(out_dir, "plass", "{sample}_plass.fa")
    output: os.path.join(out_dir, "cdhit", "{sample}_plass_cdhit{perc_id}.fa")
    log: os.path.join(logs_dir, "cdhit", "{sample}_plass_cdhit{perc_id}.log")
    benchmark: os.path.join(logs_dir, "cdhit", "{sample}_plass_cdhit{perc_id}.benchmark")
    resources: 
        mem_mb=16000 # 16GB
    threads: 10
    params:
        coverage=60,
        perc_id=lambda w: {w.perc_id},
        word_size=5, # 5 for thresholds 0.7->1.0
        #prefix=lambda w: os.path.dirname(output),
        extra=""
    conda: os.path.join(wrappers_dir, "cdhit-env.yml")
    shell:
        """
         cd-hit -i {input} -o {output} -c {params.perc_id} -T {threads} -M {resources.mem_mb} {params.extra} &> {log}
        """
         #cd-hit -i {input} -o {output} -c {params.perc_id} -n {params.word_size}  -d 0 -aS {params.coverage} \
         #-aL {params.coverage} -T {threads} -M {resources.mem_mb} {params.extra} &> {log}
        #mv {params.prefix} {output[0]} 2>> {log}

rule plass_sourmash_compute:
    input: rules.plass_eliminate_identical_contigs.output 
    output: os.path.join(out_dir, "plass", "sigs", "{sample}_plass_cdhit100_scaled{scaled}_{encoding}_k{ksize}.sig")
    params:
        #k=[5,7,11,13,15,17,19,21,25,31,35,41],
        k= lambda w: w.ksize,
        scaled= lambda w: w.scaled,
        #compute_moltypes=["protein", "dayhoff", "hp"],
        compute_moltypes= lambda w: w.encoding,
        input_is_protein=True,
        track_abundance=True,
    log: os.path.join(logs_dir, "sourmash", "{sample}_plass_cdhit100_scaled{scaled}_{encoding}_k{ksize}_compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_plass_cdhit100_scaled{scaled}_{encoding}_k{ksize}_compute.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compute.wrapper.py")

# build all compare matrices: np and csv output
rule plass_sourmash_compare_cosine: # use abundances
    input: sigs=expand(os.path.join(out_dir, "plass", "sigs", "{sample}_plass_cdhit100_scaled{{scaled}}_{{encoding}}_k{{k}}.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "plass", "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_cosine_compare.np"),
        csv=os.path.join(out_dir, "plass", "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_cosine_compare.csv")
    params:
        include_encodings = lambda w: w.encoding,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.k,
    log: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_cosine_compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_cosine_compare.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule plass_sourmash_compare_jaccard: # ignore abundances
    input: sigs=expand(os.path.join(out_dir, "plass", "sigs", "{sample}_plass_cdhit100_scaled{{scaled}}_{{encoding}}_k{{k}}.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "plass", "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_jaccard_compare.np"),
        csv=os.path.join(out_dir, "plass", "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_jaccard_compare.csv")
    params:
        ignore_abundance=True,
        include_encodings = lambda w: w.encoding,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.k,
    log: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_jaccard_compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_jaccard_compare.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule plass_sourmash_plot_cosine:
    input: rules.plass_sourmash_compare_cosine.output.np
    output: os.path.join(out_dir, "plass", "sourmash", "plots", "plass100_scaled{scaled}_{encoding}_k{k}_cosine_compare.np.matrix.pdf"), 
    params:
        plot_dir=os.path.join(out_dir, "plass", "sourmash", "plots")
    log: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_cosine_plot.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_cosine_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """

rule plass_sourmash_plot_jaccard:
    input: rules.plass_sourmash_compare_jaccard.output.np
    output: os.path.join(out_dir, "plass", "sourmash", "plots", "plass100_scaled{scaled}_{encoding}_k{k}_jaccard_compare.np.matrix.pdf"), 
    params:
        plot_dir=os.path.join(out_dir, "plass", "sourmash", "plots")
    log: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_jaccard_plot.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_{encoding}_k{k}_jaccard_plot.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """


#rule compare_tsne:
    #input: rules.sourmash_compare_jaccard.output.csv
#    input: rules.plass_sourmash_compare_jaccard.output.csv
    #output: os.path.join(out_dir, "tsne", "codingpep_k{k}_{encoding}_jaccard_tsne.pdf")
#    output: os.path.join(out_dir, "tsne", "codingpep_k{k}_{encoding}_jaccard_tsne.pdf")
        #report("plots/celltype-tsne.seed={seed}.pdf", caption="../report/celltype-tsne.rst", category="Dimension Reduction")
#    log: os.path.join(logs_dir, "codingpep_k{k}_{encoding}_jaccard_tsne.log")
#    benchmark: os.path.join(logs_dir, "codingpep_k{k}_{encoding}_jaccard_tsne.benchmark")
#    conda: os.path.join(wrappers_dir, "sklearn.yml")
#    script: os.path.join(wrappers_dir, "tsne.py")
