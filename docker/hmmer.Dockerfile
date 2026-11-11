# hmmsearch only, for x86 HPC via Singularity. Build: docker build --platform linux/amd64 -f docker/hmmer.Dockerfile -t kfpfam-hmmer:3.4-amd64 docker/
FROM debian:trixie-slim
# procps gives Nextflow its ps
RUN apt-get update && apt-get install -y --no-install-recommends hmmer procps && rm -rf /var/lib/apt/lists/*
