FROM debian:trixie-slim@sha256:020c0d20b9880058cbe785a9db107156c3c75c2ac944a6aa7ab59f2add76a7bd

RUN apt-get update \
    && apt-get install -y --no-install-recommends macaulay2=1.24.11+ds-5 \
    && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["M2", "--script"]
