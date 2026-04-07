FROM --platform=linux/amd64 python:3.12-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# patch for linbidn11
RUN ln -s /usr/lib/x86_64-linux-gnu/libidn2.so /usr/lib/x86_64-linux-gnu/libidn.so.11

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

CMD ["bash"]
