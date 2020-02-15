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
#ksize = int(snakemake.params.get("ksize", 7))
#moltype = snakemake.params.get("molecule", "protein")

#outhashes = str(snakemake.output.hashes)
outsig = str(snakemake.output.sig)

all_mins = {"protein":{},"dayhoff":{}, "hp":{}, "DNA":{} }
# structure: moltype : {ksize: mins}
# modified from https://github.com/taylorreiter/ibd/blob/master/Snakefile
for sigfile in sigs:
    if os.path.getsize(sigfile) > 0:
        sigfp = open(sigfile, 'rt')
        siglist = list(signature.load_signatures(sigfp))  ### since I store mult k in same file --> work with single ksize here??
        # here, we need to make counters by ksize and molecule type. Set these as params, rather than doing all in this script
        #loaded_sig = siglist[1]
        for sig in siglist:
            ksize = sig.minhash.ksize
            is_protein = sig.minhash.is_protein
            if is_protein:
                moltype="protein"
            else:
                if sig.minhash.dayhoff:
                    moltype="dayhoff"
                elif sig.minhash.hp:
                    moltype="hp"
                else:
                    moltype="DNA"
            print(f"molecule: {moltype}, ksize: {ksize}")
            mins = sig.minhash.get_mins() # Get the minhashes
            if ksize in all_mins[moltype].keys():
                prev_mins = all_mins[moltype][ksize]
                all_mins[moltype][ksize] = prev_mins + mins
            else:
                all_mins[moltype][ksize] = mins
            #all_mins += mins

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
                print(f"{hashval}:ct")
                if ct < min_count:
                    counts.pop(hashval)

            # now need to store this somewhere again, or directly create a sig
            # let's try building a sig. we will use this sig later to intersect with sample-specific sigs
            new_mins = set(counts.keys())
            print(len(new_mins))
            if len(new_mins) > 0:
                import pdb;pdb.set_trace()
                minhash = MinHash(n=0, ksize=ksize, scaled=scaled) # scaled=1 so we keep all (though these were previously at some other scaled val)
                minhash.add_many(set(counts.keys()))
                # write sig to file
                sigobj= sourmash.SourmashSignature(minhash, name=f"aggregated_hashvals_above_{min_count}", filename=f"generated with drop_unique_hashes.py")
                import pdb;pdb.set_trace()
                sigobjs +=[sigobj]

with open(outsig, 'wt') as sigout:
    sourmash.save_signatures(sigobjs, sigout)
                    #notify('wrote signature to {}', args.output)

# write out hashes to a text file
#with open(outhashes, "w") as out:
#    for key in counts:
#        out.write(key + '\n')


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



