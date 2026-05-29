# Dockerfile.slim 設計

- 日付: 2026-05-29
- 対象: `ddbj_mss_tools` の軽量コンテナイメージ
- 方針: 既存 `Dockerfile` を残したまま、最小構成の `Dockerfile.slim` を新規追加する

## 背景・課題

現行 `nigyta/ddbj_mss_tools:latest` は **719MB** と大きい。`docker history` による内訳:

| レイヤー | サイズ | 備考 |
|---|---:|---|
| `apt install build-essential wget` | 315MB | 削減対象（依存は全て wheel 配布のためコンパイラ不要） |
| python:3.12-slim ベース | 144MB | ベース（削減不可） |
| `pip install`（pandas/numpy/biopython 等） | 181MB | 全4ツールに必要なため維持 |
| `COPY examples/` | 78.6MB | ランタイム不要（削減対象） |
| NCBI バイナリ（asn2gb/asn2fsa） | 23.5MB | egapx2mss に必要なため維持 |
| libidn11 + wget DL | 0.7MB | NCBI バイナリの実行時依存 |

最大の無駄は **build-essential (315MB)** と **examples/ (78.6MB)**。

## 要件

- 全4ツール（egapx2mss / mss_builder / mss2ff / batch_wgs_builder）がそのまま動作すること
- egapx2mss はオフラインでも即動く（NCBI バイナリをイメージに焼き込む）
- 既存 `Dockerfile` は変更しない（後方互換）

## 設計

### `Dockerfile.slim`（シングルステージ）

```dockerfile
FROM --platform=linux/amd64 python:3.12-slim

WORKDIR /app

# NCBI バイナリ取得用に wget と libidn11 を入れ、焼き込み後に wget を撤去。
# build-essential は入れない（依存は全て wheel 配布でコンパイラ不要）。
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

# 依存インストール（wheel 配布なのでコンパイラ不要）
COPY requirements.txt pyproject.toml ./
RUN pip install --no-cache-dir hatchling \
 && pip install --no-cache-dir -r requirements.txt

# ソースのみコピー（examples/ tests/ は含めない）
COPY src/ src/
RUN pip install --no-cache-dir -e .

CMD ["bash"]
```

### NCBI バイナリの焼き込み先（重要な制約）

`asn_tools.py` の既定ディレクトリは:

```python
DEFAULT_BIN_DIR = Path(__file__).parent.parent.parent / "bin"
```

`pip install -e .`（editable）では `__file__` = `/app/src/egapx2mss/asn_tools.py` のため、
3階層上 = `/app` → `DEFAULT_BIN_DIR = /app/bin`。
よって NCBI バイナリは **`/app/bin`** に焼き込む必要がある（現行 Dockerfile と同じ位置）。
WORKDIR を `/app` から変えない限り一致が保たれる。

### `.dockerignore`（新規）

ビルドコンテキストを軽量化する。**既存 `Dockerfile` は `examples/`・`tests/` を `COPY` するため、
これらは除外しない**（除外すると既存ビルドが壊れる）。`Dockerfile.slim` は `src/` 等のみを
`COPY` するので、examples を除外しなくても slim イメージには焼き込まれない。

除外対象は「どちらの Dockerfile も COPY しない大容量データ・生成物・VCS」に限定する:

```
.git/
deprecated/
docs/
BATCH/
*.fna
*.fna.gz
SAMD*.ddbj
*.ff
test.json
__pycache__/
*.pyc
.devcontainer/
```

※ ルート直下の生成物（`GRU.ann`/`GRU.fa`/`GCA_*.fna.*` など）も対象に含めるが、
`src/`・`examples/` 配下の正規ファイルを誤って除外しないよう、ワイルドカードは
ルート限定（先頭 `/`）で指定するか、拡張子パターンの影響範囲を実装時に確認する。

## 期待される効果

| 削減項目 | 削減量 |
|---|---:|
| build-essential 削除 | −315MB |
| examples/ 非コピー | −78.6MB |
| tests/ 非コピー | 0（元々 0B） |
| wget purge | わずか |

**719MB → 約 329MB（約 54% 減）**（実測値）

## リスク・検証項目

1. **build-essential 無しで `pip install` が通るか**
   pandas / numpy / biopython / pydantic-core はいずれも manylinux wheel を配布しているため
   コンパイラ不要のはず。ビルド時に実際に成功するか要検証。
2. **NCBI バイナリの位置一致**
   ビルド後に `egapx2mss --help` 等ではなく、`DEFAULT_BIN_DIR` が `/app/bin` を指すこと、
   および `bin/asn2gb`・`bin/asn2fsa` が存在し実行可能であることを確認。
3. **`.dockerignore` と既存 `Dockerfile` の整合**
   `.dockerignore` には `examples/`・`tests/` を**含めない**（既存 `Dockerfile` がそれらを COPY する
   ため、除外すると既存ビルドが壊れる）。除外は両 Dockerfile が COPY しない大容量データ・生成物・
   `.git/`・`docs/` 等に限定する。ワイルドカード（`*.fa` 等）が `src/`・`examples/` 配下の正規
   ファイルを巻き込まないよう、ルート限定パターンを使うか影響範囲をビルド時に確認する。
4. **全4ツールの起動確認**
   ビルドしたイメージ内で 4 つの CLI が `--help` を返すこと（import エラーが無いこと）を確認。

## スコープ外

- マルチステージ化（方針②、追加削減 ~25MB だが複雑化のため不採用）
- NCBI バイナリの実行時 DL 委譲（方針③、オフライン要件と相反するため不採用）
- 既存 `Dockerfile` のリファクタ
