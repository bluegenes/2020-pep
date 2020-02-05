"""Snakemake wrapper for Transdecoder Predict"""

__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2019, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

work_dir = os.path.dirname(str(snakemake.input.longorfs))

gff = snakemake.output.get("gff", "")
bed = snakemake.output.get("bed", "")
pep = snakemake.output.get("pep", "")
cds = snakemake.output.get("cds", "")
outfiles = [gff, bed, pep, cds]

assert (
    list(map(bool, outfiles)).count(True) >= 1
), "please specify at least one output file by using 'gff', 'bed', 'pep', or 'cds' keywords in the output field)"

addl_outputs = ""
pfam = snakemake.input.get("pfam_hits", "")
if pfam:
    addl_outputs += " --retain_pfam_hits " + pfam

blast = snakemake.input.get("blastp_hits", "")
if blast:
    addl_outputs += " --retain_blastp_hits " + blast

input_fasta = str(snakemake.input.fasta)
if input_fasta.endswith("gz"):
    input_fa = os.path.basename(input_fasta.rsplit(".gz", 1)[0])
#    shell("gunzip -c {input_fasta} > {input_fa}")
else:
    input_fa = os.path.basename(input_fasta)

shell("TransDecoder.Predict --output_dir {work_dir} -t {input_fasta} --final_output_dir {work_dir} {addl_outputs} {extra} {log}")

default_output_prefix = os.path.join(work_dir, f"{input_fa}.transdecoder.")

if gff:
    gff_default = default_output_prefix + "gff3"
    shell("mv {gff_default} {gff}")
if cds:
    cds_default = default_output_prefix + "cds"
    shell("mv {cds_default} {cds}")
if bed:
    bed_default = default_output_prefix + "bed"
    shell("mv {bed_default} {bed}")
if pep:
    pep_default = default_output_prefix + "pep"
    shell("mv {pep_default} {pep}")

#shell("mv pipeliner* {output_dir}")
