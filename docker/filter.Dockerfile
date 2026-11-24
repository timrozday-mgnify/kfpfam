# kfpfam's filter, for running the prefilter outside this repo (the raw-reads pipeline's kfpfam
# engine). Build from the repository root so the build sees filter/ and measure/:
#   docker build -f docker/filter.Dockerfile -t kfpfam-filter .
# On an HPC without Docker, convert it or build the image there:
#   singularity build kfpfam-filter.sif docker-daemon://kfpfam-filter:latest
# prep.py is included so the Pfam prep directory can be built with the same image, and pyhmmer and pyarrow so the
# pipeline can run HMMER on the filter's candidates as they come (its KFPFAM_HMMER).
FROM rust:1-slim-trixie AS build
WORKDIR /src
COPY filter/ /src/filter/
RUN cd filter && cargo build --release && cp target/release/filter /usr/local/bin/filter

FROM python:3.12-slim-trixie
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir numpy pyhmmer==0.12.3 pyarrow==21.0.0
COPY --from=build /usr/local/bin/filter /usr/local/bin/filter
COPY measure/prep.py /opt/kfpfam/prep.py
ARG KFPFAM_VERSION=unknown
RUN echo "$KFPFAM_VERSION" > /opt/kfpfam/VERSION
WORKDIR /work
