

#snakemake -s orthopep.snakefile  --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} -J {cluster.jobname} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -D {cluster.chdir} " --cluster-config /home/ntpierce/2019-protein-work/farm_config.yml --jobs 5 --latency-wait=25 --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete -k 

#snakemake -s orthopep.snakefile  --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} -J {cluster.jobname} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -D {cluster.chdir} " --cluster-config /home/ntpierce/2020-pep/orthopep_sbatch_config.yml --jobs 10 --latency-wait=25 --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete -k --nolock


#snakemake -s orthopep.snakefile --latency-wait=25 --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete --nolock -n # -k --nolock #--cores 20
#snakemake -s orthopep.snakefile --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete --nolock --cores 1   # -k --nolock #--cores 20

#snakemake -s nhbds.snakefile --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete --nolock #-n #--cores 10 # -k --nolock #--cores 20
#snakemake -s nhbds.snakefile --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete --nolock --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} -J {cluster.jobname} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -D {cluster.chdir} " --cluster-config /home/ntpierce/2020-pep/orthopep_sbatch_config.yml --jobs 10 --latency-wait=25 -k

#snakemake -s preprocess.snakefile --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete --nolock --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} -J {cluster.jobname} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -D {cluster.chdir} " --cluster-config /home/ntpierce/2020-pep/orthopep_sbatch_config.yml --jobs 10 --latency-wait=25 -k

#snakemake -s annotate.snakefile --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete --nolock --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} -J {cluster.jobname} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -D {cluster.chdir} " --cluster-config /home/ntpierce/2020-pep/orthopep_sbatch_config.yml --jobs 10 --latency-wait=25 -k

snakemake -s hashclust-eukbusco.snakefile --use-conda -p --rerun-incomplete --nolock --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} -J {cluster.jobname} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -D {cluster.chdir} " --cluster-config /home/ntpierce/2020-pep/orthopep_sbatch_config.yml --jobs 10 --latency-wait=25 -k


