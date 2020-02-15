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
        #expand(os.path.join(out_dir, "transdecoder", "{sample}_trinity.transdecoder.pep"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}_search_oh0/results.csv"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}_search_oh0/{sample}.cdbg_ids.contigs.fa.gz"), sample=["MMETSP0286"], k=ksizes, r=radiuses)
        #expand(os.path.join(sgc_configdir, "{sample}_k{k}_r{r}.yml"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "preprocess", "{sample}_1.polyAtrim.fq.gz"), sample=SAMPLES)
        #expand(os.path.join(out_dir, "plass", "{sample}_plass.cdhit100.fa"), sample=SAMPLES),
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}/first_doms.txt"), sample=SAMPLES, k=ksizes, r=radiuses),
        #expand(os.path.join(out_dir, "blast", "{sample}_plass.cdhit100.psq"),sample=SAMPLES)
        #expand(os.path.join(out_dir, "spacegraphcats", "{sample}_k{k}_r{r}", "{sample}_sgc_contigs_x_plass.paladin.sort.bam{end}"), sample=SAMPLES, k=ksizes, r=radiuses, end=[".bai", ".flagstat"]),
        #expand(os.path.join(out_dir, "plass", "{sample}_split", "split_fasta.done"), sample=SAMPLES)
        #expand(os.path.join(out_dir, "cdhit", "{sample}_plass_cdhit{perc_id}.fa"), sample=SAMPLES, perc_id=[".95", ".97", ".85", ".80", ".70"])
        expand(os.path.join(out_dir, "plass", "{sample}_plass_cdhit100_scaled{scaled}.sig"), sample=SAMPLES, scaled= ["1", "50", "200", "500", "1000","2000"])


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

rule sourmash_compute_plass:
    input: rules.plass_eliminate_identical_contigs.output 
    output: os.path.join(out_dir, "plass", "{sample}_plass_cdhit100_scaled{scaled}.sig")
    params:
        k=[5,7,11,13,15,17,19,21],
        scaled= lambda w: w.scaled,
        compute_moltypes=["protein", "dayhoff", "hp"],
        input_is_protein=True,
        track_abundance=True,
    log: os.path.join(logs_dir, "sourmash", "{sample}_plass_cdhit100_scaled{scaled}_compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_plass_cdhit100_scaled{scaled}_compute.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compute.wrapper.py")

# build all compare matrices: np and csv output
rule plass_sourmash_compare_cosine: # use abundances
    input: sigs=expand(os.path.join(out_dir, "plass", "{sample}_plass_cdhit100_scaled{{scaled}}.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "plass", "sourmash", "plass100_scaled{scaled}_k{k}_{encoding}_cosine_compare.np")
        csv=os.path.join(out_dir, "plass", "sourmash", "plass100_scaled{scaled}_k{k}_{encoding}_cosine_compare.csv")
    params:
        include_encodings = lambda w: f"{w.encoding}",
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: f"{w.k}",
    log: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_k{k}_{encoding}_cosine_compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_k{k}_{encoding}_cosine_compare.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule plass_sourmash_compare_jaccard: # ignore abundances
    input: sigs=expand(os.path.join(out_dir, "plass", "{sample}_plass_cdhit100_scaled{{scaled}}.sig"), sample=SAMPLES)
    output:
        np=os.path.join(out_dir, "plass", "sourmash", "plass100_scaled{scaled}_k{k}_{encoding}_jaccard_compare.np")
        csv=os.path.join(out_dir, "plass", "sourmash", "plass100_scaled{scaled}_k{k}_{encoding}_jaccard_compare.csv")
    params:
        ignore_abundance=True,
        include_encodings = lambda w: w.encoding,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.k,
    log: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_k{k}_{encoding}_jaccard_compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "plass100_scaled{scaled}_k{k}_{encoding}_jaccard_compare.benchmark")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "sourmash-compare.wrapper.py")

rule compare_tsne:
    #input: rules.sourmash_compare_jaccard.output.csv
    input: rules.plass_sourmash_compare_jaccard.output.csv
    #output: os.path.join(out_dir, "tsne", "codingpep_k{k}_{encoding}_jaccard_tsne.pdf")
    output: os.path.join(out_dir, "tsne", "codingpep_k{k}_{encoding}_jaccard_tsne.pdf")
        #report("plots/celltype-tsne.seed={seed}.pdf", caption="../report/celltype-tsne.rst", category="Dimension Reduction")
    log: os.path.join(logs_dir, "codingpep_k{k}_{encoding}_jaccard_tsne.log")
    benchmark: os.path.join(logs_dir, "codingpep_k{k}_{encoding}_jaccard_tsne.benchmark")
    conda: os.path.join(wrappers_dir, "sklearn.yml")
    script: os.path.join(wrappers_dir, "tsne.py")

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

