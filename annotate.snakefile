"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s annotate.snakefile --use-conda #--cores 26 --cluster "sbatch -t 10:00:00 -N 11 -n 26 -p bmm --mem=60gb" --jobs 5 -k --rerun-incomplete 
"""
## protein annotation
### mmseqs2
### hhseach2
###
## nucleotide annotation
import os
import pandas as pd

database_dir= "databases"
pep_dir="/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"

sample_namelist = config.get("samplelist", "ten_haptophytes.txt")
SAMPLES = [x.strip().split('\t')[0] for x in open(sample_namelist, "r")]
#SAMPLES.remove("MMETSP0009")
#SAMPLES = ["MMETSP0143"]
SAMPLES = ["MMETSP0091"]

out_dir = config.get("out_dir", "orthopep_out")
data_dir = os.path.join(out_dir, "input_data")
logs_dir = os.path.join(out_dir, "logs")
wrappers_dir = "wrappers"
#sgc_configdir = os.path.join(out_dir, "spacegraphcats", "sgc_config")

pep_dir="/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"

#ksizes= config.get("ksizes", ["31"])
#radiuses= config.get("radiuses", ["1"])


rule all:
    input: 
        #expand(os.path.join(out_dir, "plass", "{sample}_plass.cdhit100_x_MERC.metaeuk.fa"), sample=SAMPLES),
        expand(os.path.join(out_dir, "trinity", "{sample}_trinity_x_MERC.metaeuk.fa"), sample=SAMPLES)


rule anNOGtate:
    input: 
         expand(os.path.join(out_dir, "{sample}.emapper.annotations"), sample=SAMPLES),
         expand(os.path.join(out_dir, "{sample}.emapper.seed_orthologs"), sample=SAMPLES)


rule createdb_MERCref:
    input: "/home/ntpierce/2020-pep/khtools_testing/MERC.fasta.gz"
    output: "/home/ntpierce/2020-pep/khtools_testing/MERC.mmseqsDB"
    conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
    log: os.path.join(logs_dir, "mmseqs", "MERC.createDB.log")
    benchmark: os.path.join(logs_dir, "mmseqs", "MERC.createDB.benchmark")
    shell:
        """
        mmseqs createdb {input} {output} 2> {log}
        """
localrules: mmseqs_createdb_plass, mmseqs_createdb_trinity

rule mmseqs_createdb_plass:
    input: os.path.join(out_dir, "plass", "{sample}_plass.cdhit100.fa")
    output: os.path.join(out_dir, "plass", "{sample}_plass.cdhit100.mmseqsDB")
    log: os.path.join(logs_dir, "mmseqs", "{sample}_plass.createDB.log")
    benchmark: os.path.join(logs_dir, "mmseqs", "{sample}_plass.createDB.benchmark")
    conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
    shell:
        """
        mmseqs createdb {input} {output} 2> {log}
        """
rule mmseqs_createdb_trinity:
    input: os.path.join(out_dir, "trinity", "{sample}_trinity.fa")
    output: os.path.join(out_dir, "trinity", "{sample}_trinity.mmseqsDB")
    log: os.path.join(logs_dir, "mmseqs", "{sample}_trinity.createDB.log")
    benchmark: os.path.join(logs_dir, "mmseqs", "{sample}_trinity.createDB.benchmark")
    conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
    shell:
        """
        mmseqs createdb {input} {output} 2> {log}
        """

rule mmseqs_createdb_transdecoder:
    input: os.path.join(out_dir, "transdecoder", "{sample}_trinity.transdecoder.pep")
    output: os.path.join(out_dir, "linclust", "{sample}_transdecoder.mmseqsDB")
    log: os.path.join(logs_dir, "mmseqs", "{sample}_transdecoder.createDB.log")
    benchmark: os.path.join(logs_dir, "mmseqs", "{sample}_transdecoder.createDB.benchmark")
    conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
    shell:
        """
        mmseqs createdb {input} {output} 2> {log}
        """
#rule mmseqs_cluster:
#    input: os.path.join(out_dir, "trinity", "{sample}_trinity.fa")
#    output: os.path.join(out_dir, "trinity", "{sample}_trinity.mmseqsDB")
    #log: os.path.join(logs_dir, "mmseqs", "{sample}_trinity.createDB.log")
    #benchmark: os.path.join(logs_dir, "mmseqs", "{sample}_trinity.createDB.benchmark")
    #conda: os.path.join(wrappers_dir, "mmseqs2-env.yml")
    #shell:
        #mmseqs linclust inDB outDB tmp
    #    """
    #    mmseqs cluster {input} {output} tmp
    #    """

#rule mmseqs_cluster

### metaeuk annotation ###
rule metaeuk_easypredict_plass:
    input:
        plassDB= rules.mmseqs_createdb_plass.output,
        refDB=rules.createdb_MERCref.output
    output: os.path.join(out_dir, "plass", "{sample}_plass.cdhit100_x_MERC.metaeuk.fa")
    log: os.path.join(logs_dir, "metaeuk", "{sample}_plass.cdhit100_x_MERC.metaeuk.log")
    benchmark: os.path.join(logs_dir, "metaeuk", "{sample}_plass.cdhit100_x_MERC.metaeuk.benchmark")
    params:
        tmpdir=os.path.join(out_dir, "plass", "tmp")
    conda: os.path.join(wrappers_dir, "metaeuk-env.yml")
    shell:
        """
        metaeuk easy-predict {input.plassDB} {input.refDB} {output} {params.tmpdir} 2> {log}
        """

rule metaeuk_easypredict_trinity:
    input:
        plassDB= rules.mmseqs_createdb_trinity.output,
        refDB=rules.createdb_MERCref.output
    output: os.path.join(out_dir, "trinity", "{sample}_trinity_x_MERC.metaeuk.fa")
    log: os.path.join(logs_dir, "metaeuk", "{sample}_trinity_x_MERC.metaeuk.log")
    benchmark: os.path.join(logs_dir, "metaeuk", "{sample}_trinity_x_MERC.metaeuk.benchmark")
    conda: os.path.join(wrappers_dir, "metaeuk-env.yml")
    params:
        tmpdir=os.path.join(out_dir, "trinity", "tmp")
    shell:
        """
        metaeuk easy-predict {input.plassDB} {input.refDB} {output} {params.tmpdir} 2> {log}
        """

### eggnog annotation ###
rule get_eggnog_dbs:
    output:
        os.path.join(database_dir, "eggnog.db"),
        os.path.join(database_dir, "eggnog_proteins.dmnd"),
    conda: "eggnog.yml"
    params:
        data_dir = database_dir
    shell:
        """
        mkdir -p {params.data_dir}
        download_eggnog_data.py -y -f --data_dir {params.data_dir}
        """


### run 1. plass assembly 2. nbhd-plass 3. new trinity-transdecoder 4. old trinity-transdecoder

def get_pep(w): # get johnson pep files
#most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
#odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    if w.sample == "MMETSP0251":
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")
    else:
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")

## to optimize on cluster, run eggnog mapper steps separately
# https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2
# https://github.com/eggnogdb/eggnog-mapper/issues/80

# this part is cpu intensive
rule run_eggnog_mapper_dmnd:
    input:
        nog_db=os.path.join(database_dir, "eggnog.db"),
        nog_dmnd=os.path.join(database_dir, "eggnog_proteins.dmnd"),
        pep=get_pep
    output:
        os.path.join(out_dir, "{sample}.emapper.seed_orthologs")
    params:
        mode="diamond",
        data_dir=directory("eggnog_data"),
        out_dir=out_dir
    log: os.path.join(out_dir, "logs", "{sample}_emapper_dmnd.log")
    benchmark: os.path.join(out_dir, "logs", "{sample}_emapper_dmnd.benchmark")
    threads: 20
    shadow: "shallow" ## this means tmpdir will be on local scratch (using --shadow-prefix) --> faster? 
    conda: "eggnog.yml"
    shell:
        """
        emapper.py --no_annot --no_file_comments -i {input.pep} --output {params.out_dir}/{wildcards.sample} -m {params.mode} --cpu {threads} --data_dir {params.data_dir} > {log} 2>&1 
        """

# this part is I/O intensive --> USE LOCAL tmpdir! (use shadow, with --shadow-prefix /$SCRATCH/ntpierce)
rule run_eggnog_mapper_annotate:
    input:
        nog_db=os.path.join(database_dir, "eggnog.db"),
        seed_orthologs=os.path.join(out_dir, "{sample}.emapper.seed_orthologs")
    output:
        os.path.join(out_dir, "{sample}.emapper.annotations"),
    params:
        mode="diamond",
        data_dir=database_dir,
        out_dir=out_dir
    log: os.path.join(out_dir, "logs", "{sample}_emapper_annot.log")
    benchmark: os.path.join(out_dir, "logs", "{sample}_emapper_annot.benchmark")
    threads: 1
    conda: "eggnog.yml"
    shadow: "shallow" ## here, it helps to have local eggnog db, write locally, then copy out to final location
    shell:
        """
        emapper.py --annotate_hits_table {input.seed_orthologs} --no_file_comments -o {params.out_dir}/{wildcards.sample} --data_dir {params.data_dir} --cpu {threads} > {log} 2>&1
        """
