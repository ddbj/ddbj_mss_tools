FROM --platform=linux/amd64 python:3.12-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    && rm -rf /var/lib/apt/lists/*

# libidn11 was removed from Debian 12; install from Debian 11 archive
RUN wget -q http://archive.debian.org/debian/pool/main/libi/libidn/libidn11_1.33-3_amd64.deb \
    && dpkg -i libidn11_1.33-3_amd64.deb \
    && rm libidn11_1.33-3_amd64.deb

# ── Layer 1: install third-party dependencies ─────────────────────────────
# Copied before source so this layer is cached on source-only changes.
COPY requirements.txt pyproject.toml ./
RUN pip install --no-cache-dir hatchling \
 && pip install --no-cache-dir -r requirements.txt

# ── Layer 2: copy source and install the package itself ───────────────────
COPY src/ src/
COPY examples/ examples/
COPY tests/ tests/

RUN pip install --no-cache-dir -e .

# ── Layer 3: download NCBI command-line tools ─────────────────────────────
# These binaries expire periodically; rebuild with --no-cache to refresh.
RUN mkdir -p bin \
 && wget -qO- https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/asn2gb.linux64.gz  | gunzip > bin/asn2gb \
 && wget -qO- https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/asn2fsa.linux64.gz | gunzip > bin/asn2fsa \
 && chmod +x bin/asn2gb bin/asn2fsa

CMD ["bash"]
