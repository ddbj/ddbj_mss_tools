# Dockerfile.slim Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 既存 `Dockerfile` を残したまま、最小構成の `Dockerfile.slim`（実測約329MB、現行719MBから約54%減）を新規追加する。

**Architecture:** シングルステージ。`python:3.12-slim` をベースに build-essential を入れず wheel 配布の依存のみ pip install。NCBI バイナリ（asn2gb/asn2fsa）は wget+libidn11 で `/app/bin` に焼き込み後 wget を purge。`src/` のみ COPY し `examples/`・`tests/` は含めない。`.dockerignore` でビルドコンテキストを軽量化するが、既存 `Dockerfile` が COPY する `examples/`・`tests/` は除外しない。

**Tech Stack:** Docker (BuildKit), python:3.12-slim, pip (manylinux wheels), NCBI cmdline tools (asn2gb/asn2fsa)

---

## Context & Constraints

- 設計ドキュメント: `docs/superpowers/specs/2026-05-29-dockerfile-slim-design.md`
- 制約:
  - 全4ツール（egapx2mss / mss_builder / mss2ff / batch_wgs_builder）が動作すること
  - egapx2mss はオフラインで動作（NCBI バイナリ焼き込み必須）
  - NCBI バイナリの焼き込み先は **`/app/bin`**（`src/egapx2mss/asn_tools.py:31` の `DEFAULT_BIN_DIR = Path(__file__).parent.parent.parent / "bin"`、editable install 時に `/app/bin` を指す）
  - WORKDIR は `/app` のまま（変えるとバイナリ位置がズレる）
  - `.dockerignore` に `examples/`・`tests/`・`src/` を**書かない**（既存 `Dockerfile` の `COPY examples/`・`COPY tests/` を壊さないため）
- **`.dockerignore` のパターン注意:** Docker のパターンは（先頭 `/` の有無に関わらず）任意の深さに適用されうる。`*.fa` と書くと `examples/mss2ff/test.fa` 等まで除外される。`examples/` には `.fa`/`.ann` が多数存在するため、**拡張子ワイルドカードは使わず、ルート直下の生成物はファイル名をフルで列挙する**。

### File Structure

| ファイル | 役割 | 操作 |
|---|---|---|
| `/.dockerignore` | ビルドコンテキストから VCS・docs・大容量データ・ルート生成物を除外 | 新規作成 |
| `/Dockerfile.slim` | 最小構成イメージ定義 | 新規作成 |
| `/Dockerfile` | 既存。**変更しない**（回帰防止の保護対象） | 不変 |
| `/CLAUDE.md` | セットアップ節に slim ビルド手順を追記 | 変更 |

---

## Task 1: `.dockerignore` の作成

**Files:**
- Create: `/.dockerignore`

- [ ] **Step 1: `.dockerignore` を作成**

ルート直下の生成物は（拡張子ワイルドカードを避け）実ファイル名で列挙する。`examples/`・`tests/`・`src/` は書かない。

```
.git
.gitignore
deprecated
docs
BATCH
.devcontainer
__pycache__
*.pyc
*.pyo

# root-level generated artifacts (listed by exact name to avoid matching examples/**)
GCA_003307255.1.fna.ann
GCA_003307255.1.fna.fa
GCF_020809275.1_ASM2080927v1_genomic.fna
GRU.ann
GRU.fa
SAMD888_COl.ann
SAMD888_COl.fa
SAMD888_COl.ddbj
test.json
```

- [ ] **Step 2: 既存 Dockerfile が壊れていないことを検証（examples が焼き込まれること）**

Run:
```bash
docker build -t ddbj-mss-tools:orig -f Dockerfile . \
 && docker run --rm ddbj-mss-tools:orig ls examples/mss2ff/SAMD00000001_Example-1.fa
```
Expected: ビルド成功し、最終行に `examples/mss2ff/SAMD00000001_Example-1.fa` が表示される（ファイルが存在）。エラー `No such file or directory` が出たら `.dockerignore` が examples を巻き込んでいる → Step 1 を修正。

- [ ] **Step 3: コミット**

```bash
git add .dockerignore
git commit -m "build: add .dockerignore to slim down build context"
```

---

## Task 2: `Dockerfile.slim` の作成

**Files:**
- Create: `/Dockerfile.slim`

- [ ] **Step 1: `Dockerfile.slim` を作成**

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

- [ ] **Step 2: ビルドが成功することを検証（build-essential 無しで pip が通ること）**

Run: `docker build -t ddbj-mss-tools:slim -f Dockerfile.slim .`
Expected: 全ステップ成功。特に `pip install -r requirements.txt`（pandas/numpy/biopython/pydantic-core）が manylinux wheel で入りコンパイル無しで完了。コンパイルエラー（`gcc: not found` 等）が出たら Task 5 のリスク手順へ。

- [ ] **Step 3: コミット**

```bash
git add Dockerfile.slim
git commit -m "feat: add minimal Dockerfile.slim (~325MB)"
```

---

## Task 3: イメージ検証

**Files:** （検証のみ、変更なし）

- [ ] **Step 1: サイズ確認**

Run: `docker images ddbj-mss-tools:slim --format '{{.Size}}'`
Expected: 400MB 未満（目標 ~325MB）。450MB を超える場合は `docker history ddbj-mss-tools:slim` で想定外の肥大レイヤーを調査。

- [ ] **Step 2: NCBI バイナリの位置と DEFAULT_BIN_DIR の一致確認**

Run:
```bash
docker run --rm ddbj-mss-tools:slim ls -l /app/bin/asn2gb /app/bin/asn2fsa
docker run --rm ddbj-mss-tools:slim python3 -c "from egapx2mss.asn_tools import DEFAULT_BIN_DIR; print(DEFAULT_BIN_DIR)"
```
Expected: 両バイナリが存在し実行ビット（`-rwxr-xr-x`）付き。2 コマンド目の出力が `/app/bin`。

- [ ] **Step 3: 全4ツールの起動確認（import エラーが無いこと）**

Run:
```bash
for t in egapx2mss mss_builder mss2ff batch_wgs_builder; do
  docker run --rm ddbj-mss-tools:slim $t --help >/dev/null 2>&1 && echo "$t OK" || echo "$t FAIL"
done
```
Expected: 4 行すべて `... OK`。`FAIL` が出たら該当ツールの import 依存（pandas 等）が欠けていないか確認。

- [ ] **Step 4: examples が slim イメージに含まれないことを確認**

Run: `docker run --rm ddbj-mss-tools:slim sh -c '[ -d examples ] && echo HAS_EXAMPLES || echo NO_EXAMPLES'`
Expected: `NO_EXAMPLES`

- [ ] **Step 5: egapx2mss の実走（オフライン動作・バイナリ焼き込みの実証）**

`examples/` はイメージ外なのでホストからマウントする。`/app/bin` を隠さないよう `/app` ではなく `/data` にマウントする。

Run:
```bash
docker run --rm --network none -v "$(pwd)/examples/egapx2mss:/data" -v /tmp/slimout:/out \
  ddbj-mss-tools:slim egapx2mss /data/annotated_genome.asn -o /out
```
Expected: `--network none`（ネット遮断）でも asn2gb/asn2fsa が動き、`/tmp/slimout` に `.ann` と `.fa` が生成される。`ls /tmp/slimout` に `annotated_genome.ann` と `annotated_genome.fa` が出ること。

---

## Task 4: ドキュメント追記とコミット

**Files:**
- Modify: `/CLAUDE.md`（`## セットアップ` 節の Docker 部分）

- [ ] **Step 1: CLAUDE.md のセットアップ節に slim 手順を追記**

`## セットアップ` 節の既存 Docker ブロックの直後に以下を追記する:

```markdown
# 軽量イメージ (Dockerfile.slim, 約325MB)
docker build -t ddbj-mss-tools:slim -f Dockerfile.slim .
# 注意: egapx2mss を使う場合、焼き込み済み NCBI バイナリ (/app/bin) を隠さないよう
#       作業ディレクトリは /app 以外にマウントする
docker run --rm -it -v $(pwd):/data -w /data ddbj-mss-tools:slim
```

- [ ] **Step 2: 追記が既存記述と整合していることを目視確認**

Run: `git --no-pager diff CLAUDE.md`
Expected: セットアップ節に上記ブロックのみが追加されている。

- [ ] **Step 3: コミット**

```bash
git add CLAUDE.md
git commit -m "docs: document Dockerfile.slim build/run in CLAUDE.md"
```

---

## Task 5: リスク対応（条件付き — Task 2 Step 2 が失敗した場合のみ）

build-essential 無しで `pip install` が失敗した場合のフォールバック。**Task 2 Step 2 が成功したら本タスクはスキップする。**

**Files:**
- Modify: `/Dockerfile.slim`

- [ ] **Step 1: 失敗した依存を特定**

ビルドログで「どのパッケージが」「ソースからビルドしようとしたか」を確認（`Building wheel for <pkg>` → `gcc: not found` 等）。

- [ ] **Step 2: build-essential を一時導入→purge する形に修正**

依存インストールの RUN を以下に置き換える（コンパイラを入れて wheel をビルドし、同一レイヤーで撤去）:

```dockerfile
COPY requirements.txt pyproject.toml ./
RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
 && pip install --no-cache-dir hatchling \
 && pip install --no-cache-dir -r requirements.txt \
 && apt-get purge -y build-essential && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*
```

- [ ] **Step 3: 再ビルドして Task 2 Step 2 / Task 3 を再検証**

Run: `docker build -t ddbj-mss-tools:slim -f Dockerfile.slim .`
Expected: 成功。サイズは build-essential 撤去により増加が最小限（同一レイヤー purge のため）。Task 3 の全 Step を再実行して合格を確認。

- [ ] **Step 4: コミット**

```bash
git add Dockerfile.slim
git commit -m "fix: build wheels with transient build-essential in Dockerfile.slim"
```

---

## Verification Strategy

- **ビルド回帰:** Task 1 Step 2 で既存 `Dockerfile` がビルド成功し examples を焼き込むことを確認（`.dockerignore` が既存ビルドを壊さない保証）。`Dockerfile` 自体は不変。
- **slim ビルド:** Task 2 Step 2 でコンパイラ無しビルド成功を確認。
- **機能:** Task 3 Step 2/3/5 で NCBI バイナリ位置・全4 CLI 起動・egapx2mss のオフライン実走を確認。
- **サイズ:** Task 3 Step 1 で 400MB 未満を確認。

## Risks

1. **build-essential 無しで pip が失敗する依存** → Task 5 で一時導入→purge に切替。
2. **`.dockerignore` のワイルドカードが examples を巻き込む** → 拡張子ワイルドカードを使わずルート生成物を実名列挙（Task 1 Step 1）。Task 1 Step 2 で実害無しを確認。
3. **libidn11 の Debian アーカイブ URL 失効** → 既存 Dockerfile と同一 URL のためリスク同等。失効時は両 Dockerfile 共通課題として別途対応。
4. **ボリュームマウントが `/app/bin` を隠す** → ドキュメント（Task 4）と実走検証（Task 3 Step 5）で `/data` マウントを徹底。
