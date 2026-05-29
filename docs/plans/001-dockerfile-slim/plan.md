# Dockerfile.slim 実装計画

**Goal:** 既存 `Dockerfile` を残したまま、最小構成の `Dockerfile.slim`（約325MB）を新規追加する。

**Architecture:** シングルステージ。`python:3.12-slim` をベースに、build-essential を入れず wheel 配布の依存のみ pip install。NCBI バイナリ（asn2gb/asn2fsa）は wget+libidn11 で `/app/bin` に焼き込み後、wget を purge。`src/` のみ COPY し `examples/`・`tests/` は含めない。`.dockerignore` でビルドコンテキストを軽量化する（既存 Dockerfile が COPY する `examples/`・`tests/` は除外しない）。

## Context & Constraints

- 設計ドキュメント: `docs/superpowers/specs/2026-05-29-dockerfile-slim-design.md`
- 変更/追加ファイル:
  - 新規: `Dockerfile.slim`
  - 新規: `.dockerignore`
  - 変更なし（保護対象）: `Dockerfile`
- 制約:
  - 全4ツール（egapx2mss / mss_builder / mss2ff / batch_wgs_builder）が動作すること
  - egapx2mss はオフラインで動作（NCBI バイナリ焼き込み必須）
  - NCBI バイナリの焼き込み先は **`/app/bin`**（`asn_tools.py:31` の `DEFAULT_BIN_DIR = Path(__file__).parent.parent.parent / "bin"`、editable install 時に `/app/bin` を指す）
  - WORKDIR は `/app` のまま（変えるとバイナリ位置がズレる）
  - `.dockerignore` に `examples/`・`tests/` を**含めない**（既存 `Dockerfile` の `COPY examples/`・`COPY tests/` を壊さないため）

## Phase 1: `.dockerignore` の作成

### Task 1.1: `.dockerignore` を新規作成

- **What:** ビルドコンテキストから VCS・ドキュメント・大容量データ・ルート直下の生成物を除外する
- **Where:** リポジトリルート `/.dockerignore`
- **How:** 以下の内容。ワイルドカードはルート限定（Docker の `*` は `/` を跨がないため `*.fa` は `examples/x.fa` に影響しない。これに依存する）:

```
.git/
.gitignore
deprecated/
docs/
BATCH/
.devcontainer/
__pycache__/
*.pyc
*.pyo

# root-level generated artifacts (do NOT use **/ — must not touch examples/)
/*.fna
/*.fna.gz
/*.fa
/*.ann
/*.ddbj
/*.ff
/test.json
```

- **重要な注意:** `examples/` と `tests/` は**書かない**（既存 Dockerfile が COPY するため）。`src/` も書かない。`/*.fa` のように先頭 `/` を付け、`examples/mss2ff/*.fa` 等が除外されないようにする。
- **Verify:**
  - `docker build -t ddbj-mss-tools:test -f Dockerfile .`（既存 Dockerfile）が成功し、`examples/` がイメージ内に存在すること:
    `docker run --rm ddbj-mss-tools:test ls examples/mss2ff/ | grep -q SAMD00000001_Example-1.fa`

## Phase 2: `Dockerfile.slim` の作成

### Task 2.1: `Dockerfile.slim` を新規作成

- **What:** 最小構成の Dockerfile を追加
- **Where:** リポジトリルート `/Dockerfile.slim`
- **How:** 以下の内容（設計ドキュメントの構成に準拠）:

```dockerfile
FROM --platform=linux/amd64 python:3.12-slim

WORKDIR /app

# Install NCBI command-line tools (asn2gb, asn2fsa) into /app/bin.
# wget + libidn11 are needed only for the download; wget is purged afterward.
# build-essential is intentionally omitted (all deps ship as manylinux wheels).
RUN apt-get update && apt-get install -y --no-install-recommends wget \
 && wget -q http://archive.debian.org/debian/pool/main/libi/libidn/libidn11_1.33-3_amd64.deb \
 && dpkg -i libidn11_1.33-3_amd64.deb \
 && rm libidn11_1.33-3_amd64.deb \
 && mkdir -p bin \
 && wget -qO- https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/asn2gb.linux64.gz  | gunzip > bin/asn2gb \
 && wget -qO- https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/asn2fsa.linux64.gz | gunzip > bin/asn2fsa \
 && chmod +x bin/asn2gb bin/asn2fsa \
 && apt-get purge -y wget && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*

# Install third-party dependencies (cached on source-only changes).
COPY requirements.txt pyproject.toml ./
RUN pip install --no-cache-dir hatchling \
 && pip install --no-cache-dir -r requirements.txt

# Copy source only (no examples/, no tests/) and install the package.
COPY src/ src/
RUN pip install --no-cache-dir -e .

CMD ["bash"]
```

- **Verify (build):** `docker build -t ddbj-mss-tools:slim -f Dockerfile.slim .` が成功すること。特に `pip install -r requirements.txt` が build-essential 無しで通ること（pandas/numpy/biopython/pydantic-core が wheel で入る）。

## Phase 3: イメージ検証

### Task 3.1: サイズ確認

- **What:** イメージサイズが目標（約325MB、少なくとも 719MB から大幅減）に収まること
- **How:** `docker images ddbj-mss-tools:slim --format '{{.Size}}'`
- **Verify:** 400MB 未満であること（目標 ~325MB）

### Task 3.2: NCBI バイナリの位置と実行可否

- **What:** バイナリが `/app/bin` にあり、`DEFAULT_BIN_DIR` がそこを指すこと
- **How:**
  ```bash
  docker run --rm ddbj-mss-tools:slim ls -l /app/bin/asn2gb /app/bin/asn2fsa
  docker run --rm ddbj-mss-tools:slim python3 -c \
    "from egapx2mss.asn_tools import DEFAULT_BIN_DIR; print(DEFAULT_BIN_DIR)"
  ```
- **Verify:** 両バイナリが存在し実行ビット付き、`DEFAULT_BIN_DIR` の出力が `/app/bin` であること

### Task 3.3: 全4ツールの起動確認

- **What:** 4 つの CLI が import エラー無しで `--help` を返すこと
- **How:**
  ```bash
  for t in egapx2mss mss_builder mss2ff batch_wgs_builder; do
    docker run --rm ddbj-mss-tools:slim $t --help >/dev/null && echo "$t OK" || echo "$t FAIL"
  done
  ```
- **Verify:** 4 つすべて `OK`

### Task 3.4: examples が含まれないこと（slim の確認）

- **What:** slim イメージに `examples/` が入っていないこと（サイズ削減の裏付け）
- **How:** `docker run --rm ddbj-mss-tools:slim sh -c 'ls examples 2>&1 || echo NO_EXAMPLES'`
- **Verify:** `NO_EXAMPLES`（または該当ディレクトリ無し）

## Phase 4: ドキュメント追記とコミット

### Task 4.1: README/CLAUDE.md に slim ビルド手順を追記（任意・軽微）

- **What:** slim イメージのビルド/実行方法を1ブロック追記
- **Where:** `CLAUDE.md` のセットアップ節（`## セットアップ` の Docker 部分）。README が別途あればそちらにも。
- **How:**
  ```bash
  # 軽量イメージ
  docker build -t ddbj-mss-tools:slim -f Dockerfile.slim .
  docker run --rm -it -v $(pwd):/app ddbj-mss-tools:slim
  ```
- **Verify:** 記述が既存のセットアップ節と整合していること
- **Note:** ボリュームマウント `-v $(pwd):/app` は `/app/bin` を上書きし焼き込みバイナリを隠す点に注意。egapx2mss を使う場合はマウント先を変える（例 `-v $(pwd):/work -w /work`）か、マウント無しで使う旨を一言添える。

### Task 4.2: コミット

- **What:** `Dockerfile.slim` と `.dockerignore`（+ 任意のドキュメント変更）をコミット
- **How:**
  ```
  git add Dockerfile.slim .dockerignore CLAUDE.md
  git commit -m "feat: add minimal Dockerfile.slim (~325MB)"
  ```
- **Verify:** `git status` がクリーン（対象ファイルのみコミット済み）、未追跡の大容量データは残ったまま

## Verification Strategy

- **ビルド検証:** 既存 `Dockerfile`（Phase 1.1 Verify）と新規 `Dockerfile.slim`（Phase 2.1 Verify）の両方がビルド成功すること。`.dockerignore` が既存ビルドを壊していないことを最優先で確認。
- **機能検証:** Phase 3.2/3.3 で NCBI バイナリ位置と全4 CLI 起動を確認。可能なら `examples/egapx2mss/` の入力で egapx2mss を 1 回実走し `.ann`/`.fa` が生成されることまで確認する（ただし examples はイメージ外なのでマウントして実行）:
  ```bash
  docker run --rm -v $(pwd)/examples/egapx2mss:/data ddbj-mss-tools:slim \
    egapx2mss /data/annotated_genome.asn -o /tmp/out
  ```
  （`/app/bin` を隠さないよう `/app` ではなく `/data` にマウントする点に注意）
- **サイズ検証:** Phase 3.1 で 400MB 未満を確認。
- **回帰:** 既存 `Dockerfile` は無変更のため機能回帰なし。`.dockerignore` のみが既存ビルドへの影響源 → Phase 1.1 で担保。

## Risks

1. **build-essential 無しで pip install が失敗する依存があった場合** → 失敗した依存だけを wheel 確認。最終手段として slim でも build-essential を一時導入→purge するマルチステージに切替（設計の方針②）。実装時に Phase 2.1 のビルドで判明する。
2. **`.dockerignore` のワイルドカードが examples を巻き込む** → 先頭 `/` 付きルート限定パターンで回避。Phase 1.1 Verify で実害が無いことを確認。
3. **libidn11 の Debian アーカイブ URL が将来失効** → 既存 Dockerfile と同一 URL のため現状はリスク同等。失効時は両 Dockerfile 共通の課題として別途対応。
