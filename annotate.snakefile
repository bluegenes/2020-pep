"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s eggnog.snakefile --use-conda #--cores 26 --cluster "sbatch -t 10:00:00 -N 11 -n 26 -p bmm --mem=60gb" --jobs 5 -k --rerun-incomplete 
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

rule anNOGtate:
    input: 
         expand(os.path.join(out_dir, "{sample}.emapper.annotations"), sample=SAMPLES),
         expand(os.path.join(out_dir, "{sample}.emapper.seed_orthologs"), sample=SAMPLES)

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
