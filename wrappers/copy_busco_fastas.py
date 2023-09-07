import os
import sys
import argparse
import glob
from shutil import copyfile


# requires python >= 3.6
# run: python copy_busco_fastas.py --source orthopep_out/busco/MMETSP0091_transdecoder_busco_euk/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences  --destination orthopep_out/busco_fastas

log = str(snakemake.log)

basename=snakemake.params.get("basename")
assert basename is not None, "basename param is required"
source_dir=snakemake.params.get("source_dir")
assert source_dir is not None, "source dir is required"
dest_dir=snakemake.params.get("dest_dir")
assert dest_dir is not None, "dest dir is required"

def copy_busco_fastas(source, destination, logfile):
    with open(logfile, "w") as out_log:
        if not os.path.exists(destination):
            try:
                os.mkdir(destination)
            except Exception as e:
                out_log.write(f"\n\tError: cannot make {destination} destination folder. Please fix.\n\n")
                sys.exit(-1)

        # glob the fasta files
        buscofastas = glob.glob(os.path.join(source, "*"))
        if len(buscofastas) < 1:
            out_log.write(f"\n\tWarning: no fastas found at source dir {source} \n\n")
            #sys.exit(-1)

        # copy each file over
        else:
            for busco_fa in buscofastas:
                new_name = basename + "_" + os.path.basename(busco_fa)
                dest = os.path.join(destination, new_name)
                if dest.endswith('.faa'):
                    dest = dest.rsplit('.faa', 1)[0] + '.fasta'
                    #dest = dest + '.fasta'
                out_log.write(f"\n\t copying {busco_fa} to {dest} \n\n")
                copyfile(busco_fa, dest)


copy_busco_fastas(source_dir, dest_dir, log)

#if __name__ == '__main__':
#    p = argparse.ArgumentParser()
#    p.add_argument('--source', default = os.getcwd())
#    p.add_argument('--destination', required=True)
#    p.add_argument('--logfile', default= "copy_busco_fastas.log")
#    args = p.parse_args()
#    sys.exit(copy_busco_fastas(args.source, args.destination, args.logfile))

