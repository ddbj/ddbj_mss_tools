# ddbj_mss_tools [beta version, under development]

DDBJ MSS (Mass Submission System) 登録ファイルを生成するPythonツール群。

- **egapx2mss** — NCBI EGAPx の出力 (ASN.1形式) を DDBJ MSS形式 (.ann / .fa) に変換
- **wgs_maker** — FASTAファイルと Excel/TSV のメタデータから DDBJ MSS形式ファイルを一括生成 (作成中)

> English documentation is available in the [second half of this page](#english).

---

## 目次

- [インストール](#インストール)
- [egapx2mss の使い方](#egapx2mss-の使い方)
  - [基本的な使い方](#基本的な使い方)
  - [オプション一覧](#オプション一覧)
  - [common JSON ファイル](#common-json-ファイル)
  - [染色体テーブル (--chromosomes)](#染色体テーブル---chromosomes)
  - [注意点](#注意点)
- [wgs_maker の使い方](#wgs_maker-の使い方)

---

## インストール

```bash
pip install -e .
```

Docker を使う場合:

```bash
docker build -t ddbj-mss-tools .
docker run --rm -it -v $(pwd):/app ddbj-mss-tools
```

---

## egapx2mss の使い方

詳細は [src/egapx2mss/](src/egapx2mss/) を参照。

### 基本的な使い方

```bash
egapx2mss input.asn \
  --common examples/egapx2mss/common_example.json \
  --output output_prefix
```

デフォルトでは入力ファイル名（拡張子なし）が出力プレフィックスになります。
上記の例では `output_prefix.ann` と `output_prefix.fa` が生成されます。

### オプション一覧

| オプション | 説明 |
|---|---|
| `input` | 入力 ASN.1 ファイル (.asn) |
| `-o`, `--output` | 出力ファイルのプレフィックス（省略時は入力ファイル名） |
| `--common` | 共通メタデータ JSON ファイル（DBLINK, SUBMITTER, REFERENCE 等） |
| `--chromosomes` | 染色体テーブル TSV ファイル（後述） |
| `--bin-dir` | asn2gb / asn2fsa バイナリの保存ディレクトリ（デフォルト: `bin/`） |
| `--keep-tmp` | 中間ファイル (.tbl, raw FASTA) を削除せず保持する |

### common JSON ファイル

`--common` で指定する JSON ファイルには、登録者情報・文献情報・BioProject/BioSample リンク等を記述します。
`DBLINK.project` と `DBLINK.sample` は必須です。

```json
{
    "DBLINK": {
        "project": "PRJD000001",
        "sample": "SAMD000001",
        "DRA": ["DRA000001"]
    },
    "SUBMITTER": {
        "ab_name": ["Tanizawa,Y."],
        "contact": "Yuki Tanizawa",
        "email": "xxx@ddbj.nig.ac.jp",
        "institute": "National Institute of Genetics",
        "country": "Japan"
    },
    "REFERENCE": [{
        "title": "Genome sequencing of ...",
        "ab_name": ["Tanizawa,Y."],
        "status": "Unpublished",
        "year": 2025
    }],
    "SOURCE": {
        "organism": "Brassica rapa",
        "mol_type": "genomic DNA",
        "cultivar": "NAPPA"
    },
    "SOURCE_MODIFIER": "cultivar",
    "ASSEMBLY_GAP": {
        "linkage_evidence": "proximity ligation",
        "min_gap_length": 10
    }
}
```

- JSON5 スタイルの末尾カンマ (trailing comma) が使えます
- サンプルファイル: [examples/egapx2mss/common_example.json](examples/egapx2mss/common_example.json)

#### SOURCE セクション

`SOURCE` に記載した qualifier がそのまま source フィーチャーに書き込まれます。
`SOURCE_MODIFIER` に qualifier 名を指定すると、`ff_definition` の生物名に続く識別子として使われます（例: `"cultivar"` → `Brassica rapa NAPPA DNA, ...`）。

#### ASSEMBLY_GAP セクション

連続する N 塩基（デフォルト10塩基以上）を自動検出し `assembly_gap` フィーチャーを挿入します。
`linkage_evidence` に指定できる値は以下の3種類です:

| 値 | gap_type | estimated_length |
|---|---|---|
| `paired-ends` | within scaffolds | known |
| `proximity ligation` | within scaffolds | unknown |
| `align genus` | within scaffolds | unknown |

### 染色体テーブル (--chromosomes)

ゲノムアセンブリの配列を染色体・オルガネラ・unplaced に分類するための5列タブ区切りファイルです。

```
# seq_id    type          seq_name       status    topology
Chr1        chromosome    1              complete  linear
ChrM        organelle     mitochondrion  complete  circular
scaffold001 unplaced                     partial   linear
```

| 列 | 内容 |
|---|---|
| seq_id | FASTA ヘッダーの配列 ID |
| type | `chromosome` / `organelle` / `unplaced` |
| seq_name | 染色体番号やオルガネラ名（unplaced の場合は空でも可） |
| status | `complete` / `partial` |
| topology | `linear` / `circular` |

- 省略した場合、全配列が unplaced として扱われ WGS モードで出力されます
- `#` で始まる行はコメントとして無視されます

### 注意点

#### asn2gb / asn2fsa バイナリの自動ダウンロード

`egapx2mss` は内部で NCBI の `asn2gb` および `asn2fsa` コマンドを使用します。
これらのバイナリは初回実行時に NCBI のFTPサーバーから自動でダウンロードされ、
`--bin-dir` で指定したディレクトリ（デフォルト: `bin/`）にキャッシュされます。

- macOS: `asn2gb.mac`, `asn2fsa.mac`
- Linux: `asn2gb.linux64`, `asn2fsa.linux64`

#### バイナリの有効期限

NCBI の `asn2gb` / `asn2fsa` には **利用期限**があります。
期限切れのバイナリを実行すると空の出力が返ります。
`egapx2mss` はこれを検出した場合、自動で最新版を再ダウンロードして1回リトライします。
リトライ後も失敗する場合は入力ファイルが正しい ASN.1 形式か確認してください。

#### catenated ASN.1 ファイル

EGAPx が出力する ASN.1 ファイルには複数の `Seq-entry ::=` ブロックが連結されています。
`asn2gb` はこの形式に対応していますが、`asn2fsa` は非対応のため、
`egapx2mss` が内部でブロックごとに一時ファイルへ分割して処理します。

---

## wgs_maker の使い方

作成中。

---

<a id="english"></a>

---

# ddbj_mss_tools (English)

A set of Python tools for generating DDBJ MSS (Mass Submission System) submission files.

- **egapx2mss** — Converts NCBI EGAPx output (ASN.1 format) to DDBJ MSS format (.ann / .fa)
- **wgs_maker** — Generates DDBJ MSS files in bulk from FASTA files and Excel/TSV metadata *(under development)*

---

## Table of Contents

- [Installation](#installation)
- [egapx2mss Usage](#egapx2mss-usage)
  - [Basic Usage](#basic-usage)
  - [Options](#options)
  - [Common JSON File](#common-json-file)
  - [Chromosome Table (--chromosomes)](#chromosome-table---chromosomes)
  - [Important Notes](#important-notes)
- [wgs_maker Usage](#wgs_maker-usage)

---

## Installation

```bash
pip install -e .
```

Using Docker:

```bash
docker build -t ddbj-mss-tools .
docker run --rm -it -v $(pwd):/app ddbj-mss-tools
```

---

## egapx2mss Usage

See [src/egapx2mss/](src/egapx2mss/) for source code.

### Basic Usage

```bash
egapx2mss input.asn \
  --common examples/egapx2mss/common_example.json \
  --output output_prefix
```

By default, the output prefix is the input filename without extension.
The above example produces `output_prefix.ann` and `output_prefix.fa`.

### Options

| Option | Description |
|---|---|
| `input` | Input ASN.1 file (.asn) |
| `-o`, `--output` | Output file prefix (default: input basename) |
| `--common` | Common metadata JSON file (DBLINK, SUBMITTER, REFERENCE, etc.) |
| `--chromosomes` | Chromosome table TSV file (see below) |
| `--bin-dir` | Directory for asn2gb / asn2fsa binaries (default: `bin/`) |
| `--keep-tmp` | Keep intermediate files (.tbl, raw FASTA) |

### Common JSON File

The JSON file specified with `--common` describes submitter information, references, and BioProject/BioSample links.
`DBLINK.project` and `DBLINK.sample` are required.

```json
{
    "DBLINK": {
        "project": "PRJD000001",
        "sample": "SAMD000001",
        "DRA": ["DRA000001"]
    },
    "SUBMITTER": {
        "ab_name": ["Tanizawa,Y."],
        "contact": "Yuki Tanizawa",
        "email": "xxx@ddbj.nig.ac.jp",
        "institute": "National Institute of Genetics",
        "country": "Japan"
    },
    "REFERENCE": [{
        "title": "Genome sequencing of ...",
        "ab_name": ["Tanizawa,Y."],
        "status": "Unpublished",
        "year": 2025
    }],
    "SOURCE": {
        "organism": "Brassica rapa",
        "mol_type": "genomic DNA",
        "cultivar": "NAPPA"
    },
    "SOURCE_MODIFIER": "cultivar",
    "ASSEMBLY_GAP": {
        "linkage_evidence": "proximity ligation",
        "min_gap_length": 10
    }
}
```

- Trailing commas (JSON5-style) are accepted.
- Sample file: [examples/egapx2mss/common_example.json](examples/egapx2mss/common_example.json)

#### SOURCE Section

Qualifiers listed under `SOURCE` are written directly into the source feature.
`SOURCE_MODIFIER` specifies which qualifier name is appended to the organism name in `ff_definition`
(e.g. `"cultivar"` → `Brassica rapa NAPPA DNA, ...`).

#### ASSEMBLY_GAP Section

Runs of consecutive N bases (10 or more by default) are automatically detected and written as `assembly_gap` features.
Valid values for `linkage_evidence` are:

| Value | gap_type | estimated_length |
|---|---|---|
| `paired-ends` | within scaffolds | known |
| `proximity ligation` | within scaffolds | unknown |
| `align genus` | within scaffolds | unknown |

### Chromosome Table (--chromosomes)

A 5-column tab-separated file that classifies sequences into chromosomes, organelles, or unplaced scaffolds.

```
# seq_id    type          seq_name       status    topology
Chr1        chromosome    1              complete  linear
ChrM        organelle     mitochondrion  complete  circular
scaffold001 unplaced                     partial   linear
```

| Column | Description |
|---|---|
| seq_id | Sequence ID from the FASTA header |
| type | `chromosome` / `organelle` / `unplaced` |
| seq_name | Chromosome number or organelle name (may be empty for unplaced) |
| status | `complete` / `partial` |
| topology | `linear` / `circular` |

- If omitted, all sequences are treated as unplaced and output in WGS mode.
- Lines beginning with `#` are treated as comments.

### Important Notes

#### Automatic Download of asn2gb / asn2fsa

`egapx2mss` uses NCBI's `asn2gb` and `asn2fsa` commands internally.
These binaries are automatically downloaded from the NCBI FTP server on first run
and cached in the directory specified by `--bin-dir` (default: `bin/`).

- macOS: `asn2gb.mac`, `asn2fsa.mac`
- Linux: `asn2gb.linux64`, `asn2fsa.linux64`

#### Binary Expiration

NCBI's `asn2gb` / `asn2fsa` binaries have an **expiration date**.
An expired binary returns empty output.
`egapx2mss` detects this and automatically re-downloads the latest version and retries once.
If the retry also fails, verify that the input file is a valid ASN.1 file.

#### Catenated ASN.1 Files

ASN.1 files output by EGAPx contain multiple concatenated `Seq-entry ::=` blocks.
While `asn2gb` supports this format natively, `asn2fsa` does not.
`egapx2mss` handles this by splitting the blocks into temporary files internally.

---

## wgs_maker Usage

Under development.
