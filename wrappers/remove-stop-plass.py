#! /usr/bin/env python
#https://raw.githubusercontent.com/spacegraphcats/2018-paper-spacegraphcats/master/pipeline-base/scripts/remove-stop-plass.py

import screed
import sys
import os

# read input files, make sure we're working with a list (even if single file)

log = snakemake.log

infiles = ([snakemake.input]
        if isinstance(snakemake.input, str)
        else snakemake.input)


with open(log, "w"):
    for filename in infiles:
        log.write("Removing stop (*) from plass contigs " + filename + "\n")
        if ".fa" in filename:
            outname = filename.rsplit(".fa", 1)[0] + ".nostop.fa"
        else:
            outname = filename + ".nostop.fa"
        with open(outname, 'wt') as fp:
            for record in screed.open(filename):
                record.sequence = record.sequence.replace('*', '')
                fp.write('>{}\n{}\n'.format(record.name, record.sequence))
