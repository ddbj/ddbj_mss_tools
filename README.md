# ddbj_mss_tools

> **本ツールは beta 版です。不明点があれば dfast[@]ddbj.nig.ac.jp までお問い合わせください。**

DDBJ MSS (Mass Submission System) 登録ファイルを生成するPythonツール群。

| ツール | 概要 |
|---|---|
| **egapx2mss** | NCBI EGAPx の出力 (ASN.1形式) を DDBJ MSS形式 (.ann / .fa) に変換 |
| **mss_builder** | FASTA ファイルから DDBJ MSS形式ファイル (.ann / .fa) を生成 |
| **batch_wgs_builder** | 複数ゲノムの DDBJ MSS WGS/MAG-WGS 登録ファイルを一括生成 |
| **mss2ff** | MSS アノテーションファイルから DDBJ Flat File を生成 |

\* **mss_builder** と **batch_wgs_builder** には生物学的注釈 (CDSなどの遺伝子アノテーション情報) を行う機能はありません。塩基配列のみを登録するためのファイルを生成します。

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
- [mss_builder の使い方](#mss_builder-の使い方)
  - [基本的な使い方](#基本的な使い方-1)
  - [オプション一覧](#オプション一覧-1)
  - [WGS モードと染色体モード](#wgs-モードと染色体モード)
- [batch_wgs_builder の使い方](#batch_wgs_builder-の使い方)
  - [基本的な使い方](#基本的な使い方-2)
  - [オプション一覧](#オプション一覧-2)
  - [TSV ファイルの形式](#tsv-ファイルの形式)
  - [common JSON ファイル](#common-json-ファイル-1)
  - [WGS と MAG-WGS](#wgs-と-mag-wgs)
- [mss2ff の使い方](#mss2ff-の使い方)
  - [基本的な使い方](#基本的な使い方-3)
  - [オプション一覧](#オプション一覧-3)

---

## インストール

```bash
git clone https://github.com/ddbj/ddbj_mss_tools.git
cd ddbj_mss_tools
pip install -e .
```

Docker を使う場合:

```bash
docker build -t ddbj-mss-tools .
docker run --rm -it -v $(pwd):/app ddbj-mss-tools
```

---

## egapx2mss の使い方

NCBI EGAPx が出力する ASN.1 ファイルを DDBJ MSS 形式の `.ann` / `.fa` ファイルに変換します。

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
        "ab_name": ["Tanizawa,Y.", "Mishima,H.", "Smith,J."],
        "contact": "Yasuhiro Tanizawa",
        "email": "xxx@ddbj.nig.ac.jp",
        "institute": "National Institute of Genetics",
        "country": "Japan",
        "city": "Mishima",
        "street": "1111 Yata",
        "zip": "411-8540"
    },
    "REFERENCE": [{
        "title": "Genome sequencing of ...",
        "ab_name": ["Tanizawa,Y.", "Mishima,H."],
        "status": "Unpublished",
    }],
    "SOURCE": {
        "organism": "Brassica rapa",
        "mol_type": "genomic DNA",
        "cultivar": "NAPPA",
        "collection_date": "2025-05-13",
        "geo_loc_name": "Japan:Shizuoka, Mishima"
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

`SOURCE_MODIFIER` には、SOURCE 内に記載した qualifier のうち、**種内での個体を識別する名称**として用いるものを指定します。
登録後の公開ファイルの DEFINITION 行に反映され、たとえば `"SOURCE_MODIFIER": "cultivar"` と指定した場合、
`Brassica rapa NAPPA DNA, chromosome 1, complete sequence.` のような形式で生物名に続けて記載されます。

#### ASSEMBLY_GAP セクション

連続する N 塩基を自動検出し、`assembly_gap` フィーチャーとしてアノテーションファイルに記載します。
`min_gap_length` 以上の長さのギャップが対象です。DDBJ では長さ10以上のギャップへの記載を求めているため、
デフォルト値 (`10`) から変更する必要は通常ありません。

`linkage_evidence` にはギャップが導入された理由を記載します。指定できる値:

| 値 | 用途 |
|---|---|
| `paired-ends` | ペアエンドリードによるスキャフォールディング |
| `proximity ligation` | Hi-C を用いたスキャフォールディング |
| `align genus` | 同属の近縁種ゲノムへのアラインメントで決定 |
| `align xgenus` | 属をまたぐ近縁種ゲノムへのアラインメントで決定 |
| `align trnscpt` | トランスクリプトのアラインメントで決定 |
| `map` | 遺伝地図・物理マップ等をもとに決定 |
| `within clone` | クローン内配列から決定 |
| `clone contig` | クローンコンティグから決定 |
| `strobe` | ストローブリードによるスキャフォールディング |
| `unspecified` | 上記以外 / 不明 |

gap_type と estimated_length の対応:

| linkage_evidence | gap_type | estimated_length |
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

## mss_builder の使い方

FASTA ファイルから DDBJ MSS 形式の `.ann` / `.fa` ファイルを生成します。
遺伝子アノテーションを含まない WGS コンティグ登録や、染色体・オルガネラを指定した complete genome 登録に使用します。

### 基本的な使い方

```bash
# WGS 登録（全配列を unplaced コンティグとして扱う）
mss_builder genome.fa \
  --common examples/mss_builder/common_example.json \
  --output output_prefix

# Complete genome 登録（染色体テーブルを指定）
mss_builder genome.fa \
  --common examples/mss_builder/common_example.json \
  --chromosomes chromosomes.tsv \
  --output output_prefix
```

### オプション一覧

| オプション | 説明 |
|---|---|
| `input` | 入力 FASTA ファイル (.fa / .fasta) |
| `-o`, `--output` | 出力ファイルのプレフィックス（省略時は入力ファイル名） |
| `--common` | 共通メタデータ JSON ファイル（egapx2mss と同形式） |
| `--chromosomes` | 染色体テーブル TSV（省略時は WGS モード） |

common JSON の形式は [egapx2mss と同じ](#common-json-ファイル)です。
`DBLINK.project` と `DBLINK.sample` が必須です。

### WGS モードと染色体モード

**WGS モード**（`--chromosomes` 省略時）:

- source フィーチャーを COMMON ブロックに `@@[entry]@@` メタ記法で記載
- 各エントリには `assembly_gap` フィーチャーのみ記載

**染色体モード**（`--chromosomes` 指定時）:

- 各エントリに独立した source フィーチャーを記載
- 染色体テーブルの `type` / `seq_name` / `status` / `topology` を反映

---

## batch_wgs_builder の使い方

サンプルリスト TSV と common JSON を入力とし、複数ゲノムの DDBJ MSS WGS/MAG-WGS 登録ファイル（`.ann` / `.fa`）を一括生成します。

### 基本的な使い方

```bash
# WGS 一括生成
batch_wgs_builder sample_list_WGS.tsv \
  --common examples/batch_wgs_builder/common_example.json \
  --out-dir output_dir

# MAG-WGS を含む場合（TSV 内の _trad_submission_category 列で指定）
batch_wgs_builder sample_list_MAG-WGS.tsv \
  --common examples/batch_wgs_builder/common_example.json \
  --out-dir output_dir
```

出力ファイルは `{biosample}_{strain_or_isolate}.ann` / `.fa` という名前になります。

### オプション一覧

| オプション | 省略形 | 説明 |
|---|---|---|
| `tsv` | — | サンプルリスト TSV（必須） |
| `--common` | `-m` | 共通メタデータ JSON ファイル |
| `--out-dir` | `-o` | 出力ディレクトリ（デフォルト: `.`） |
| `--hold-date` | `-H` | 公開保留日（YYYYMMDD 形式） |

### TSV ファイルの形式

ヘッダー行が2行あります。1行目がフィーチャー名、2行目が qualifier 名に対応します。

```
_           DBLINK    DBLINK      DBLINK                  ST_COMMENT        ...  source    source  ...  COMMENT
_file_path  project   biosample   sequence read archive   Assembly Method   ...  organism  strain  ...  line
path/to/genome.fa.gz  PRJDB99999  SAMD999997  DRR999997  Skesa v. 1.0      ...  Homo sapiens  HG001  ...  Comment text
```

| 列ヘッダー (行1 / 行2) | 内容 |
|---|---|
| `_` / `_file_path` | FASTA ファイルへのパス（必須） |
| `_` / `_trad_submission_category` | 登録カテゴリ（`MAG-WGS` を指定、省略時は `WGS`） |
| `DBLINK` / `project` | BioProject ID |
| `DBLINK` / `biosample` | BioSample ID |
| `DBLINK` / `sequence read archive` | DRA アクセッション（`;` 区切りで複数指定可） |
| `ST_COMMENT` / `Assembly Method` 等 | Genome-Assembly-Data の qualifier |
| `source` / `organism` 等 | source フィーチャーの qualifier |
| `COMMENT` / `line` | COMMENT ブロックの内容（`;` 区切りで複数行） |

- サンプルファイル: [examples/batch_wgs_builder/sample_list_WGS.tsv](examples/batch_wgs_builder/sample_list_WGS.tsv)

### common JSON ファイル

SUBMITTER, REFERENCE, ASSEMBLY_GAP, SOURCE_MODIFIER など、全サンプルに共通するメタデータを記載します。
DBLINK や source フィーチャーの情報は TSV で指定するため **DBLINK は不要**ですが、
共通値を書いておくことも可能で、その場合は TSV の値で上書きされます。

```json
{
    "SUBMITTER": {
        "ab_name": ["Suzuki,K.", "Doe,J."],
        "contact": "Jane Doe",
        "email": "xxx@ddbj.nig.ac.jp",
        "institute": "National Institute of Genetics",
        "country": "Japan",
        "city": "Mishima",
        "street": "Yata 1111",
        "zip": "411-8540"
    },
    "REFERENCE": [{
        "ab_name": ["Suzuki,K.", "Doe,J."],
        "status": "Unpublished",
        "title": "Genome sequences for ..."
    }],
    "SOURCE_MODIFIER": "strain",
    "ASSEMBLY_GAP": {
        "linkage_evidence": "paired-ends",
        "min_gap_length": 10
    }
}
```

- サンプルファイル: [examples/batch_wgs_builder/common_example.json](examples/batch_wgs_builder/common_example.json)

### WGS と MAG-WGS

TSV の `_trad_submission_category` 列で登録カテゴリを指定します。

| カテゴリ | DATATYPE | DIVISION | KEYWORD |
|---|---|---|---|
| `WGS`（デフォルト） | WGS | — | WGS, STANDARD_DRAFT |
| `MAG-WGS` | WGS | ENV | ENV, WGS, STANDARD_DRAFT, Metagenome Assembled Genome, MAG |

MAG-WGS の場合、source フィーチャーに `environmental_sample`（値なし）が自動付加されます。
`metagenome_source` は TSV の `source` 列として記載します。

---

## mss2ff の使い方

MSS アノテーションファイル（`.ann` または `.annt.tsv`）と FASTA ファイルから DDBJ Flat File を生成します。

### 基本的な使い方

```bash
# 基本
mss2ff annotation.ann --fasta genome.fa --division BCT --output output.ff

# 全オプション指定
mss2ff annotation.ann \
    --fasta genome.fa \
    --output output.ff \
    --division BCT \
    --submission-date 2025-04-01 \
    --file-date 2025-04-01 \
    --email your@email.com \
    --accession AP000001
```

### オプション一覧

| オプション | 省略形 | デフォルト | 説明 |
|---|---|---|---|
| `ANN` | — | 必須 | MSS アノテーションファイル (.ann または .annt.tsv) |
| `--fasta` | `-f` | なし | FASTA シーケンスファイル（CDS 翻訳に必要） |
| `--output` | `-o` | 標準出力 | 出力ファイルパス |
| `--division` | `-d` | `UNK` | DDBJ division コード（BCT, VRL, PLN 等） |
| `--submission-date` | `-s` | 今日 | REFERENCE 1 の投稿日（YYYY-MM-DD） |
| `--file-date` | — | 今日 | LOCUS 行のファイル作成日（YYYY-MM-DD） |
| `--email` | — | `mss2ff@ddbj.nig.ac.jp` | NCBI Entrez API 用メールアドレス |
| `--accession` | `-a` | なし | 開始アクセッション番号 |
| `--no-taxonomy` | — | false | NCBI 分類情報の取得をスキップ |

#### アクセッション番号の形式

| 形式 | 例 |
|---|---|
| 2文字プレフィックス + 6桁以上 | `AP000001` |
| 4文字プレフィックス + 2桁バージョン + 6桁以上 | `AAXJ010000001` |
| 6文字プレフィックス + 2桁バージョン + 6桁以上 | `AAXJEM010000001` |

エントリーはアノテーションファイルの順に serial を +1 ずつ割り当てます。

---

<a id="english"></a>

---

# ddbj_mss_tools (English)

A set of Python tools for generating DDBJ MSS (Mass Submission System) submission files.

| Tool | Description |
|---|---|
| **egapx2mss** | Converts NCBI EGAPx output (ASN.1 format) to DDBJ MSS format (.ann / .fa) |
| **mss_builder** | Generates DDBJ MSS files (.ann / .fa) from a single FASTA file |
| **batch_wgs_builder** | Batch-generates DDBJ MSS WGS/MAG-WGS submission files for multiple genomes |
| **mss2ff** | Generates DDBJ Flat Files from MSS annotation files |

---

## Table of Contents

- [Installation](#installation)
- [egapx2mss Usage](#egapx2mss-usage)
  - [Basic Usage](#basic-usage)
  - [Options](#options)
  - [Common JSON File](#common-json-file)
  - [Chromosome Table (--chromosomes)](#chromosome-table---chromosomes)
  - [Important Notes](#important-notes)
- [mss_builder Usage](#mss_builder-usage)
  - [Basic Usage](#basic-usage-1)
  - [Options](#options-1)
  - [WGS Mode and Chromosome Mode](#wgs-mode-and-chromosome-mode)
- [batch_wgs_builder Usage](#batch_wgs_builder-usage)
  - [Basic Usage](#basic-usage-2)
  - [Options](#options-2)
  - [TSV File Format](#tsv-file-format)
  - [Common JSON File](#common-json-file-1)
  - [WGS and MAG-WGS](#wgs-and-mag-wgs)
- [mss2ff Usage](#mss2ff-usage)
  - [Basic Usage](#basic-usage-3)
  - [Options](#options-3)

---

## Installation

```bash
git clone https://github.com/ddbj/ddbj_mss_tools.git
cd ddbj_mss_tools
pip install -e .
```

Using Docker:

```bash
docker build -t ddbj-mss-tools .
docker run --rm -it -v $(pwd):/app ddbj-mss-tools
```

---

## egapx2mss Usage

Converts ASN.1 files produced by NCBI EGAPx into DDBJ MSS format `.ann` / `.fa` files.

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

`SOURCE_MODIFIER` specifies which qualifier in `SOURCE` is used as the **intraspecific identifier** for the organism.
It is reflected in the DEFINITION line of the published flat file — for example, `"SOURCE_MODIFIER": "cultivar"` produces a definition like
`Brassica rapa NAPPA DNA, chromosome 1, complete sequence.`

#### ASSEMBLY_GAP Section

Runs of consecutive N bases are automatically detected and written as `assembly_gap` features in the annotation file.
Only gaps of `min_gap_length` bases or longer are annotated. DDBJ requires annotation of gaps ≥ 10 bases, so the default value (`10`) does not normally need to be changed.

`linkage_evidence` describes how the gap was introduced. Valid values:

| Value | When to use |
|---|---|
| `paired-ends` | Scaffolding with paired-end reads |
| `proximity ligation` | Scaffolding with Hi-C |
| `align genus` | Determined by alignment to a congeneric genome |
| `align xgenus` | Determined by alignment to a genome from another genus |
| `align trnscpt` | Determined by transcript alignment |
| `map` | Determined from a genetic or physical map |
| `within clone` | Determined from within-clone sequence |
| `clone contig` | Determined from a clone contig |
| `strobe` | Scaffolding with strobe reads |
| `unspecified` | Other / unknown |

Corresponding gap_type and estimated_length:

| linkage_evidence | gap_type | estimated_length |
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

## mss_builder Usage

Generates DDBJ MSS format `.ann` / `.fa` files from a FASTA file.
Use this for WGS contig submissions without gene annotations, or for complete genome submissions with chromosome/organelle assignments.

### Basic Usage

```bash
# WGS submission (all sequences treated as unplaced contigs)
mss_builder genome.fa \
  --common examples/mss_builder/common_example.json \
  --output output_prefix

# Complete genome submission (with chromosome table)
mss_builder genome.fa \
  --common examples/mss_builder/common_example.json \
  --chromosomes chromosomes.tsv \
  --output output_prefix
```

### Options

| Option | Description |
|---|---|
| `input` | Input FASTA file (.fa / .fasta) |
| `-o`, `--output` | Output file prefix (default: input basename) |
| `--common` | Common metadata JSON file (same format as egapx2mss) |
| `--chromosomes` | Chromosome table TSV (if omitted, WGS mode is used) |

The common JSON format is the [same as egapx2mss](#common-json-file).
`DBLINK.project` and `DBLINK.sample` are required.

### WGS Mode and Chromosome Mode

**WGS mode** (no `--chromosomes`):

- The source feature is written in the COMMON block using `@@[entry]@@` meta-notation.
- Each entry contains only `assembly_gap` features (if `ASSEMBLY_GAP` is configured).

**Chromosome mode** (`--chromosomes` specified):

- A separate source feature is written per entry.
- Chromosome/organelle names, topology, and completeness are derived from the chromosome table.

---

## batch_wgs_builder Usage

Batch-generates DDBJ MSS WGS/MAG-WGS submission files (`.ann` / `.fa`) for multiple genomes from a sample list TSV and a common JSON file.

### Basic Usage

```bash
# WGS batch generation
batch_wgs_builder sample_list_WGS.tsv \
  --common examples/batch_wgs_builder/common_example.json \
  --out-dir output_dir
```

Output files are named `{biosample}_{strain_or_isolate}.ann` / `.fa`.

### Options

| Option | Short | Description |
|---|---|---|
| `tsv` | — | Sample list TSV (required) |
| `--common` | `-m` | Common metadata JSON file |
| `--out-dir` | `-o` | Output directory (default: `.`) |
| `--hold-date` | `-H` | Public release hold date (YYYYMMDD) |

### TSV File Format

The TSV has two header rows: row 1 contains feature names, row 2 contains qualifier names.

```
_           DBLINK    DBLINK      DBLINK                  ST_COMMENT        ...  source    source  ...  COMMENT
_file_path  project   biosample   sequence read archive   Assembly Method   ...  organism  strain  ...  line
path/to/genome.fa.gz  PRJDB99999  SAMD999997  DRR999997  Skesa v. 1.0      ...  Homo sapiens  HG001  ...  Comment text
```

| Header (row 1 / row 2) | Description |
|---|---|
| `_` / `_file_path` | Path to the FASTA file (required) |
| `_` / `_trad_submission_category` | Submission category (`MAG-WGS`; default is `WGS` if omitted) |
| `DBLINK` / `project` | BioProject ID |
| `DBLINK` / `biosample` | BioSample ID |
| `DBLINK` / `sequence read archive` | DRA accession(s) (semicolon-separated for multiple) |
| `ST_COMMENT` / `Assembly Method` etc. | Genome-Assembly-Data qualifiers |
| `source` / `organism` etc. | source feature qualifiers |
| `COMMENT` / `line` | COMMENT block content (semicolon-separated for multiple lines) |

Sample file: [examples/batch_wgs_builder/sample_list_WGS.tsv](examples/batch_wgs_builder/sample_list_WGS.tsv)

### Common JSON File

Describes metadata common to all samples: SUBMITTER, REFERENCE, ASSEMBLY_GAP, SOURCE_MODIFIER, etc.
**DBLINK is not required** here (it is specified per sample in the TSV), but common DBLINK or SOURCE values may be included and will be overridden by TSV values.

```json
{
    "SUBMITTER": {
        "ab_name": ["Suzuki,K.", "Doe,J."],
        "contact": "Jane Doe",
        "email": "xxx@ddbj.nig.ac.jp",
        "institute": "National Institute of Genetics",
        "country": "Japan",
        "city": "Mishima",
        "street": "Yata 1111",
        "zip": "411-8540"
    },
    "REFERENCE": [{
        "ab_name": ["Suzuki,K.", "Doe,J."],
        "status": "Unpublished",
        "title": "Genome sequences for ..."
    }],
    "SOURCE_MODIFIER": "strain",
    "ASSEMBLY_GAP": {
        "linkage_evidence": "paired-ends",
        "min_gap_length": 10
    }
}
```

Sample file: [examples/batch_wgs_builder/common_example.json](examples/batch_wgs_builder/common_example.json)

### WGS and MAG-WGS

The submission category is specified in the `_trad_submission_category` column of the TSV.

| Category | DATATYPE | DIVISION | KEYWORD |
|---|---|---|---|
| `WGS` (default) | WGS | — | WGS, STANDARD_DRAFT |
| `MAG-WGS` | WGS | ENV | ENV, WGS, STANDARD_DRAFT, Metagenome Assembled Genome, MAG |

For MAG-WGS, `environmental_sample` (no value) is automatically added to the source feature.
`metagenome_source` should be specified as a `source` column in the TSV.

---

## mss2ff Usage

Generates a DDBJ Flat File from an MSS annotation file (`.ann` or `.annt.tsv`) and a FASTA file.

### Basic Usage

```bash
# Basic
mss2ff annotation.ann --fasta genome.fa --division BCT --output output.ff

# All options
mss2ff annotation.ann \
    --fasta genome.fa \
    --output output.ff \
    --division BCT \
    --submission-date 2025-04-01 \
    --file-date 2025-04-01 \
    --email your@email.com \
    --accession AP000001
```

### Options

| Option | Short | Default | Description |
|---|---|---|---|
| `ANN` | — | required | MSS annotation file (.ann or .annt.tsv) |
| `--fasta` | `-f` | none | FASTA sequence file (required for CDS translation) |
| `--output` | `-o` | stdout | Output file path |
| `--division` | `-d` | `UNK` | DDBJ division code (BCT, VRL, PLN, etc.) |
| `--submission-date` | `-s` | today | Submission date for Reference 1 (YYYY-MM-DD) |
| `--file-date` | — | today | File creation date for LOCUS line (YYYY-MM-DD) |
| `--email` | — | `mss2ff@ddbj.nig.ac.jp` | Email address for NCBI Entrez API calls |
| `--accession` | `-a` | none | Starting accession number |
| `--no-taxonomy` | — | false | Skip NCBI taxonomy lookup |

#### Accession Number Formats

| Format | Example |
|---|---|
| 2-letter prefix + ≥6 digits | `AP000001` |
| 4-letter prefix + 2-digit version + ≥6 digits | `AAXJ010000001` |
| 6-letter prefix + 2-digit version + ≥6 digits | `AAXJEM010000001` |

Entries are assigned serial numbers sequentially in the order they appear in the annotation file.
