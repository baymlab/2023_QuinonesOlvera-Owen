snakemake \
--snakefile figure_ndiv.smk \
-k \
-p \
--use-conda \
--conda-frontend mamba \
--cluster-config config.o2.json \
--cluster "sbatch -p {cluster.p} --mem {cluster.mem} -n {cluster.n} -J {cluster.J} -t {cluster.time} -o {cluster.o} -e {cluster.e}" \
-j 50 \
--latency-wait 30