

snakemake -s orthopep.snakefile  --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -D {cluster.chdir} " --cluster-config /home/ntpierce/2019-protein-work/farm_config.yml --jobs 1 --latency-wait=25 --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete -k 


#snakemake -s orthopep.snakefile --latency-wait=25 --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete -k --nolock #--cores 10


