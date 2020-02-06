#! /usr/bin/env python

import os
import sys
import screed

log = str(snakemake.log)
fasta = str(snakemake.input)

default_prefix = os.path.basename(fasta).rsplit('.fa',1)[0]

prefix = snakemake.params.get("prefix", default_prefix)
output_touchfile = str(snakemake.output)
output_dirname = os.path.dirname(output_touchfile)
min_len=int(snakemake.params.get("min_len", 45))

with open(log, "w") as out_log:
    i=0
    out_log.write(f"Splitting {fasta} into one-file-per-contig. Writing files to {output_dirname} \n")
    for record in screed.open(fasta):
        if len(record.sequence) >= min_len:
            name = record.name.rsplit(" ", 1)[0]
            outfile = os.path.join(output_dirname, prefix + str(i) + ".fa")
            with open(outfile, "w") as out:
                out.write(f">{name}\n{record.sequence}\n")
            i+=1
    out_log.write(f"{str(i)} contigs written as individual fasta files\n")

with open(output_touchfile, 'w') as tf:
    tf.write("done!\n")
