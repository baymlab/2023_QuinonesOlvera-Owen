snakemake \
-j 50 \
-p \
--use-conda \
--rerun-incomplete \
--conda-frontend mamba \
--cluster-config config.o2.json \
--cluster "sbatch -p {cluster.p} --mem {cluster.mem} -n {cluster.n} -J {cluster.J} -t {cluster.time} -o {cluster.o} -e {cluster.e}" \
--latency-wait 45 \
-k
