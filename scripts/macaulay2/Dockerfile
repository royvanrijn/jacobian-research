FROM debian:trixie-slim@sha256:020c0d20b9880058cbe785a9db107156c3c75c2ac944a6aa7ab59f2add76a7bd

ARG TARGETARCH

RUN apt-get update \
    && case "$TARGETARCH" in \
         amd64) m2_version='1.24.11+ds-5+b1' ;; \
         arm64) m2_version='1.24.11+ds-5' ;; \
         *) echo "Unsupported Macaulay2 architecture: $TARGETARCH" >&2; exit 1 ;; \
       esac \
    && apt-get install -y --no-install-recommends "macaulay2=$m2_version" \
    && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["M2", "--script"]
