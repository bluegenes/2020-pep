import os
import sys
import argparse
import sourmash
from sourmash import MinHash, SourmashSignature, signature
from sourmash import sourmash_args
from sourmash.logging import notify, error
from collections import Counter
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

sigs = snakemake.input
assert sigs is not None, "signatures are required as input"
min_count = int(snakemake.params.get("min_count", 2))
scaled = int(snakemake.params.get("scaled", 1)) # default scaled=1 keeps all hashes
# if these are not set, assume multiple ksizes/moltypes per file, do for each
ksize = int(snakemake.params.get("ksize"))
moltype = snakemake.params.get("molecule")

outhashes = str(snakemake.output.hashes)
outsig = str(snakemake.output.sig)

all_mins = {"protein":{},"dayhoff":{}, "hp":{}, "DNA":{} }
# structure: moltype : {ksize: mins}
# modified from https://github.com/taylorreiter/ibd/blob/master/Snakefile
for sigfile in sigs:
    if os.path.getsize(sigfile) > 0:
        sigfp = open(sigfile, 'rt')
        siglist = list(signature.load_signatures(sigfp))
        for sig in siglist:
            # ksize
            sig_ksize = sig.minhash.ksize
            if ksize: # if we set a single ksize in the params
                if sig_ksize != ksize:
                    continue
            # moltype
            if sig.minhash.is_protein:
                sig_moltype="protein"
            else:
                if sig.minhash.dayhoff:
                    sig_moltype="dayhoff"
                elif sig.minhash.hp:
                    sig_moltype="hp"
                else:
                    sig_moltype="DNA"
            if moltype:
                if sig_moltype != moltype:
                    continue
            print(f"molecule: {moltype}, ksize: {ksize}")
            mins = sig.minhash.get_mins() # Get the minhashes
            if sig_ksize in all_mins[sig_moltype].keys():
                prev_mins = all_mins[sig_moltype][sig_ksize]
                all_mins[sig_moltype][sig_ksize] = prev_mins + mins
            else:
                all_mins[sig_moltype][sig_ksize] = mins

sigobjs = []
for moltype, mins_by_k in all_mins.items():
    for ksize, mins in mins_by_k.items(): # here we can check if k, moltype are in user-defined criteria...?
        #counts = Counter(all_mins) # tally the number of hashes
        if len(mins) > 0: # anything to count?
            print(len(mins))
            #print(mins)
            counts = Counter(mins) # tally the number of hashes
            # remove hashes that occur only once
            for hashval, ct in counts.copy().items():
                print(f"{hashval}:{ct}")
                if ct < min_count:
                    counts.pop(hashval)
            # write out hashes

            # let's try building a sig. we will use this sig later to intersect with sample-specific sigs
            new_mins = set(counts.keys())
            print(len(new_mins))
            with open(outhashes, "w") as out:
                for hsh in new_mins:
                    out.write(str(hsh) + '\n')
            if len(new_mins) > 0:
                minhash = MinHash(n=0, ksize=ksize, scaled=scaled) # scaled=1 so we keep all (though these were previously at some other scaled val)
                minhash.add_many(set(counts.keys()))
                # write sig to file
                sigobj= sourmash.SourmashSignature(minhash, name=f"aggregated_hashvals_above_{min_count}", filename=f"generated with drop_unique_hashes.py")
                sigobjs +=[sigobj]

## this part only handles one output file -- doesn't take care of case with many ksizes/moltypes
with open(outsig, 'wt') as sigout:
    sourmash.save_signatures(sigobjs, sigout)
    #notify('wrote signature to {}', args.output)

# write out hashes to a text file


# this part is from
# https://github.com/dib-lab/sourmash/blob/7661087aa0b0e81bfec82a58002463d7c699528a/utils/hashvals-to-signature.py

#ksize = int(snakemake.params.get("ksize", 7))
#do some checking here?
#if scaled==0:
#    num=int(snakemake.params.get("num_hashes", 0))
#    if num==0:
#        notify('setting --num automatically from the number of hashes.')
#        num = len(counts.keys()) # can you access keys this was from Counter object?

#import pdb;pdb.set_trace()
# construct empty MinHash object according to args

# add hashes into!
#minhash.add_many(hashes)

#if len(minhash) < len(hashes):
    #notify("WARNING: loaded {} hashes, but only {} made it into MinHash.",
    #       len(hashes), len(minhash))
    #if scaled:
#        notify("This is probably because of the scaled argument.")
#    elif args.num:
#        notify("This is probably because your --num is set to {}",
#               args.num)

#if num > len(minhash):
#    notify("WARNING: --num set to {}, but only {} hashes in signature.",
#           num, len(minhash))



