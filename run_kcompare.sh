
snakemake -s khtools_compare.snakefile --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete --nolock --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} -J {cluster.jobname} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -D {cluster.chdir} " --cluster-config /home/ntpierce/2020-pep/orthopep_sbatch_config.yml --jobs 1 --latency-wait=25 -k