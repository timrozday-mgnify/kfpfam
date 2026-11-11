# hmmsearch (truth), MMseqs2 and nail (baselines). Build: docker build --platform linux/arm64 -t kfpfam-tools docker/
FROM debian:trixie-slim
ARG ARCH=arm64
ARG NAIL_ARCH=aarch64
RUN apt-get update && apt-get install -y --no-install-recommends hmmer ca-certificates curl time python3 && rm -rf /var/lib/apt/lists/*
# nail is tested against MMseqs2 15-6f452, so it gets its own copy; 18-8cc5c is the current release.
RUN for v in 15-6f452 18-8cc5c; do mkdir -p /opt/mmseqs-$v \
      && curl -fsSL https://github.com/soedinglab/MMseqs2/releases/download/$v/mmseqs-linux-$ARCH.tar.gz \
       | tar xz --strip-components=1 -C /opt/mmseqs-$v; done
RUN mkdir -p /opt/nail/bin \
    && curl -fsSL https://github.com/TravisWheelerLab/nail/releases/download/nail-v0.7.1/nail-$NAIL_ARCH-unknown-linux-gnu.tar.gz \
     | tar xz -C /opt/nail \
    && mv "$(find /opt/nail -type f \( -name nail -o -name 'nail-*' \) | head -1)" /opt/nail/bin/nail \
    && chmod +x /opt/nail/bin/nail
ENV PATH=/opt/mmseqs-18-8cc5c/bin:/opt/nail/bin:$PATH
WORKDIR /work
