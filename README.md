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
  - [.tbl / .fa から直接変換する](#tbl--fa-から直接変換する)
  - [common JSON ファイル](#common-json-ファイル)
  - [Sequence role ファイル (--sequence_roles)](#sequence-role-ファイル---sequence_roles)
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
  - [登録カテゴリ (_submission_category)](#登録カテゴリ-_submission_category)
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

`gff2mss` サブツールは追加の `ddbj-gff` パッケージを必要とします（他ツールは不要）。gff2mss を使う場合のみ、extra を指定してインストールしてください:

```bash
pip install -e ".[gff2mss]"
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
# 通常変換（出力は input.ann / input.fa）
egapx2mss input.asn --common examples/egapx2mss/common_example.json

# 出力ディレクトリとファイル名を指定
egapx2mss input.asn \
  --common examples/egapx2mss/common_example.json \
  --outdir results/ \
  --prefix output
```

デフォルトでは入力ファイルと同じディレクトリに、入力ファイル名（拡張子なし）をプレフィックスとして出力します。
上記2番目の例では `results/output.ann` と `results/output.fa` が生成されます。

### オプション一覧

| オプション | 説明 |
|---|---|
| `input` | 入力 ASN.1 ファイル (.asn)。`--tbl` と `--fsa` を両方指定する場合は省略可 |
| `-o`, `--outdir` | 出力先ディレクトリ（存在しない場合は自動作成。デフォルト: 入力ファイルと同じディレクトリ） |
| `-p`, `--prefix` | 出力ファイルのベースネーム（ディレクトリ区切り文字不可。デフォルト: 入力ファイルのベースネーム） |
| `--common` | 共通メタデータ JSON ファイル（DBLINK, SUBMITTER, REFERENCE 等） |
| `--submission_category` | 登録カテゴリ（`WGS`, `GNM`, `MAG-WGS` 等）。JSON の `_submission_category` を上書き（後述） |
| `--sequence_roles` | Sequence role ファイル TSV（後述）。旧名 `--chromosomes` も互換のため受け付けます |
| `--bin-dir` | asn2gb / asn2fsa バイナリの保存ディレクトリ（デフォルト: `~/.local/share/ddbj_mss_tools/bin`） |
| `--keep-tmp` | 中間ファイル (.tbl, raw FASTA) を削除せず保持する |
| `--tbl` | 既存の NCBI feature table (.tbl) を直接指定（step 1/3 をスキップ） |
| `--fsa` | 既存の FASTA ファイル (.fa/.fsa) を直接指定（step 2/3 をスキップ） |
| `--preconvert-only` | step 1/3・2/3 のみ実行して終了（.tbl と .fa を生成） |

### .tbl / .fa から直接変換する

`--tbl` と `--fsa` を両方指定すると、ASN.1 ファイルを省略して既存の中間ファイルから step 3/3（MSS アノテーション変換）のみ実行できます。
この場合 `asn2gb` / `asn2fsa` のダウンロード・実行は不要です。

```bash
# step 3/3 だけ実行（既存の .tbl と .fa を使用）
egapx2mss --tbl input.tbl --fsa input.fa \
  --common common_example.json \
  --outdir results/ --prefix output
```

また、`--preconvert-only` を使うと step 1/3・2/3 だけ実行して止めることができます。

```bash
# step 1/3・2/3 だけ実行（.tbl と .fa を生成して終了）
egapx2mss input.asn --preconvert-only --outdir tmp/
```

同じ入力に対して再実行した場合、既に生成済みの `.tbl` や `.fa` が存在するステップは自動的にスキップされます。

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
    "SOURCE_IDENTIFIER": "cultivar",
    "ASSEMBLY_GAP": [
        {
            "enabled": true,
            "linkage_evidence": "proximity ligation",
            "min_gap_length": 100,
            "max_gap_length": 100,
            "gap_type": "within scaffold",
            "estimated_length": "unknown"
        },
        {
            "enabled": true,
            "linkage_evidence": "paired-ends",
            "min_gap_length": 10,
            "gap_type": "within scaffold",
            "estimated_length": "known"
        }
    ]
}
```

- JSON5 スタイルの末尾カンマ (trailing comma) が使えます
- サンプルファイル: [examples/egapx2mss/common_example.json](examples/egapx2mss/common_example.json)

#### _submission_category

`_submission_category` キーで登録カテゴリを指定すると、そのカテゴリに応じた DATATYPE / DIVISION / KEYWORD が自動注入され、必要な source qualifier が不足している場合には警告が表示されます。

```json
{
    "_submission_category": "GNM",
    ...
}
```

コマンドラインオプション `--submission_category` を指定すると、JSON の値を上書きできます。
対応カテゴリの一覧は [登録カテゴリ一覧](#登録カテゴリ-_submission_category) を参照してください。

#### SOURCE セクション

`SOURCE` に記載した qualifier がそのまま source フィーチャーに書き込まれます。

`SOURCE_IDENTIFIER` には、SOURCE 内に記載した qualifier のうち、**種内での個体を識別する名称**として用いるものを指定します。
登録後の公開ファイルの DEFINITION 行に反映され、たとえば `"SOURCE_IDENTIFIER": "cultivar"` と指定した場合、
`Brassica rapa NAPPA DNA, chromosome 1, complete sequence.` のような形式で生物名に続けて記載されます。
旧名 `INFRASPECIFIC_NAME_MODIFIER` も後方互換性のため引き続き受け付けます。

#### ASSEMBLY_GAP セクション

連続する N 塩基を自動検出し、`assembly_gap` フィーチャーとしてアノテーションファイルに記載します。
`ASSEMBLY_GAP` は **配列（配列形式）** で記載し、複数のルールを優先順に指定できます。
各 N-run に対して、リストの先頭から順に条件を照合し、最初にマッチしたルールが適用されます。
結果は座標順に出力されます。

各ルールのフィールド:

| フィールド | 必須 | デフォルト | 説明 |
|---|---|---|---|
| `enabled` | — | `true` | `false` にするとこのルールをスキップ |
| `linkage_evidence` | ✓ | — | ギャップ導入の根拠（下表参照） |
| `min_gap_length` | — | `10` | アノテーション対象の最小ギャップ長 |
| `max_gap_length` | — | 上限なし | アノテーション対象の最大ギャップ長 |
| `gap_type` | — | 推奨値 | `gap_type` qualifier の値 |
| `estimated_length` | — | 推奨値 | `estimated_length` qualifier の値（`known` または `unknown`） |

`gap_type` と `estimated_length` を省略した場合は `linkage_evidence` に応じた推奨値が自動設定されます:

| linkage_evidence | gap_type | estimated_length |
|---|---|---|
| `paired-ends` | within scaffold | known |
| `proximity ligation` | within scaffold | unknown |
| `align genus` | within scaffold | unknown |

`linkage_evidence` に指定できる値:

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

**記載例（Hi-C スキャフォールディング + ペアエンドの2ルール）:**

```json
"ASSEMBLY_GAP": [
    {
        "enabled": true,
        "linkage_evidence": "proximity ligation",
        "min_gap_length": 100,
        "max_gap_length": 100,
        "gap_type": "within scaffold",
        "estimated_length": "unknown"
    },
    {
        "enabled": true,
        "linkage_evidence": "paired-ends",
        "min_gap_length": 10,
        "gap_type": "within scaffold",
        "estimated_length": "known"
    }
]
```

この例では、長さ 100 の N-run には最初のルール (`unknown`) が適用され、
長さ 10〜99 の N-run には2番目のルール (`known`) が適用されます。
長さ 9 以下の N-run はどのルールにもマッチしないためアノテーションされません。

### Sequence role ファイル (--sequence_roles)

ゲノムアセンブリの配列を染色体・オルガネラ・プラスミド・セグメント・unplaced に分類するための5列タブ区切りファイル（通称 `sequence_roles.tsv`）です。
以前は「染色体テーブル」と呼ばれており、オプションも `--chromosomes` でした。旧名は互換のため引き続き受け付けます。

```
# seq_id    type          seq_name       status    topology
Chr1        chromosome    1              complete  linear
ChrM        organelle     mitochondrion  complete  circular
scaffold001 unplaced                     partial   linear
```

| 列 | 内容 |
|---|---|
| seq_id | FASTA ヘッダーの配列 ID |
| type | `chromosome` / `organelle` / `plasmid` / `segment` / `unplaced` |
| seq_name | 染色体番号やオルガネラ名。`plasmid`、および submission 全体で2件以上ある `chromosome` / `segment` では必須（空欄はエラーになります）。単一 chromosome・単一 segment・unplaced では空でも可 |
| status | `complete` / `partial` |
| topology | `linear` / `circular` |

- 省略した場合、全配列が unplaced として扱われ WGS モードで出力されます
- `#` で始まる行はコメントとして無視されます
- source の ff_definition（DEFINITION 行）は DDBJ MSS のメタ記法（`@@[qualifier_name]@@`）を使ったテンプレートとして出力され、登録時に同じ source フィーチャーの qualifier 実値に展開されます
- `segment`（分節ゲノムのセグメント）は、submission 全体で `segment` が1件のみの場合は `complete genome` / `partial genome`（source に `/segment` は付与されない）、複数件ある場合は `segment @@[segment]@@, complete sequence` / `segment @@[segment]@@, unlocalized sequence @@[entry]@@`（source に `/segment` を付与）として出力されます
- type / count / status ごとの ff_definition の完全な決定テーブルは以下のとおりです。`{P}` = `@@[organism]@@ @@[{source_identifier}]@@`（`source_identifier`＝`SOURCE_IDENTIFIER` の qualifier 名。空/None なら `@@[organism]@@` のみ）、`{mol}` は mol_type 由来の具体値（DNA / RNA / tRNA / rRNA / mRNA）:

| type | count | status | ff_definition |
|------|-------|--------|---------------|
| unplaced（entry=None） | — | is_wgs=true | `{P} {mol}, @@[submitter_seqid]@@` |
| unplaced（entry=None） | — | is_wgs=false | `{P} {mol}, unplaced sequence @@[entry]@@` |
| chromosome | count==1 | complete | `{P} {mol}, chromosome, complete genome` |
| chromosome | count≥2 | complete | `{P} {mol}, chromosome @@[chromosome]@@, complete sequence` |
| chromosome | count==1 | partial | `{P} {mol}, chromosome, partial genome` |
| chromosome | count≥2 | partial | `{P} {mol}, chromosome @@[chromosome]@@, unlocalized sequence @@[entry]@@` |
| organelle | — | complete | `{P} {organelle_code} {mol}, complete genome` |
| organelle | — | partial | `{P} {organelle_code} {mol}, partial genome` |
| plasmid | — | complete | `{P} plasmid @@[plasmid]@@ {mol}, complete sequence` |
| plasmid | — | partial | `{P} plasmid @@[plasmid]@@ {mol}, partial sequence` |
| segment | count==1 | complete | `{P} {mol}, complete genome` |
| segment | count==1 | partial | `{P} {mol}, partial genome` |
| segment | count≥2 | complete | `{P} {mol}, segment @@[segment]@@, complete sequence` |
| segment | count≥2 | partial | `{P} {mol}, segment @@[segment]@@, unlocalized sequence @@[entry]@@` |
| その他（未知 type） | — | — | `{P} {mol}, @@[entry]@@` |

  - `@@[chromosome]@@` / `@@[plasmid]@@` / `@@[segment]@@` を出力する分岐（chromosome count≥2 / plasmid / segment count≥2）では `seq_name` が必須で、空だと `ValueError` になります。
  - `{organelle_code}` は `/organelle` 値を DEFINITION 用の形容詞形に変換した値（`mitochondrion`→`mitochondrial`、`plastid:chloroplast`→`chloroplast` 等）で、これはメタ記法化しません。

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
mss_builder genome.fa --common examples/mss_builder/common_example.json

# 出力ディレクトリとファイル名を指定
mss_builder genome.fa \
  --common examples/mss_builder/common_example.json \
  --outdir results/ --prefix output

# Complete genome 登録（sequence role ファイルを指定）
mss_builder genome.fa \
  --common examples/mss_builder/common_example.json \
  --sequence_roles sequence_roles.tsv \
  --outdir results/ --prefix output
```

### オプション一覧

| オプション | 説明 |
|---|---|
| `input` | 入力 FASTA ファイル (.fa / .fasta) |
| `-o`, `--outdir` | 出力先ディレクトリ（存在しない場合は自動作成。デフォルト: 入力ファイルと同じディレクトリ） |
| `-p`, `--prefix` | 出力ファイルのベースネーム（ディレクトリ区切り文字不可。デフォルト: 入力ファイルのベースネーム） |
| `--common` | 共通メタデータ JSON ファイル（egapx2mss と同形式） |
| `--submission_category` | 登録カテゴリ（`WGS`, `GNM`, `MAG-WGS` 等）。JSON の `_submission_category` を上書き（後述） |
| `--sequence_roles` | Sequence role ファイル TSV（省略時は WGS モード）。旧名 `--chromosomes` も互換のため受け付けます |

common JSON の形式は [egapx2mss と同じ](#common-json-ファイル)です。
`DBLINK.project` と `DBLINK.sample` が必須です。

### WGS モードと染色体モード

**WGS モード**（`--sequence_roles` 省略時。旧名 `--chromosomes` も同義）:

- source フィーチャーを COMMON ブロックに `@@[entry]@@` メタ記法で記載
- 各エントリには `assembly_gap` フィーチャーのみ記載

**染色体モード**（`--sequence_roles` 指定時。旧名 `--chromosomes` も同義）:

- 各エントリに独立した source フィーチャーを記載
- Sequence role ファイルの `type` / `seq_name` / `status` / `topology` を反映

---

## batch_wgs_builder の使い方

サンプルリスト TSV と common JSON を入力とし、複数ゲノムの DDBJ MSS WGS/MAG-WGS 登録ファイル（`.ann` / `.fa`）を一括生成します。

### 基本的な使い方

```bash
# WGS 一括生成
batch_wgs_builder sample_list_WGS.tsv \
  --common examples/batch_wgs_builder/common_example.json \
  --out-dir output_dir

# MAG-WGS を含む場合（TSV 内の _submission_category 列または --submission_category で指定）
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
| `--submission_category` | — | 登録カテゴリ（`WGS`, `MAG-WGS` 等）。JSON および TSV 各行の `_submission_category` を一括上書き |

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
| `_` / `_submission_category` | 登録カテゴリ（`WGS`, `MAG-WGS` 等。省略時は `WGS`。後述） |
| `DBLINK` / `project` | BioProject ID |
| `DBLINK` / `biosample` | BioSample ID |
| `DBLINK` / `sequence read archive` | DRA アクセッション（`;` 区切りで複数指定可） |
| `ST_COMMENT` / `Assembly Method` 等 | Genome-Assembly-Data の qualifier |
| `source` / `organism` 等 | source フィーチャーの qualifier |
| `COMMENT` / `line` | COMMENT ブロックの内容（`;` 区切りで複数行） |

- サンプルファイル: [examples/batch_wgs_builder/sample_list_WGS.tsv](examples/batch_wgs_builder/sample_list_WGS.tsv)

### common JSON ファイル

SUBMITTER, REFERENCE, ASSEMBLY_GAP, SOURCE_IDENTIFIER など、全サンプルに共通するメタデータを記載します。
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
    "SOURCE_IDENTIFIER": "strain",
    "ASSEMBLY_GAP": [
        {
            "enabled": true,
            "linkage_evidence": "paired-ends",
            "min_gap_length": 10,
            "gap_type": "within scaffold",
            "estimated_length": "known"
        }
    ]
}
```

`ASSEMBLY_GAP` の詳細は [egapx2mss の ASSEMBLY_GAP セクション](#assembly_gap-セクション)を参照してください。

- サンプルファイル: [examples/batch_wgs_builder/common_example.json](examples/batch_wgs_builder/common_example.json)

### 登録カテゴリ (_submission_category)

登録カテゴリは以下の3か所で指定できます（優先度は上から高い順）:

1. **CLI オプション** `--submission_category CATEGORY` — TSV・JSON の値を一括上書き
2. **TSV** の `_submission_category` 列 — 行ごとに指定（`batch_wgs_builder` のみ）
3. **common JSON** の `_submission_category` キー — JSON 内のデフォルト値

カテゴリを指定すると、DATATYPE / DIVISION / KEYWORD が自動注入されます。
必須フィールドが不足している場合は実行時に警告が表示され、空値で補完されます。

| カテゴリ | DATATYPE | DIVISION | 主な KEYWORD | source_identifier | 主な用途 |
|---|---|---|---|---|---|
| （未指定） | — | — | — | — | カテゴリなし |
| `GNM` | — | — | — | — | コンプリートゲノム |
| `WGS` | WGS | — | WGS, STANDARD_DRAFT | — | ドラフトゲノム（WGS） |
| `ENV` | — | ENV | ENV | — | 環境 DNA |
| `MAG` | — | ENV | ENV, MAG, Metagenome Assembled Genome | `isolate` | MAG（コンプリートゲノム） |
| `MAG-WGS` | WGS | ENV | ENV, MAG, Metagenome Assembled Genome, WGS, STANDARD_DRAFT | `isolate` | MAG（ドラフトゲノム） |
| `TSA` | — | TSA | TSA, Transcriptome Shotgun Assembly | — | Transcriptome Shotgun Assembly |
| `TPA` | TPA | — | TPA, Third Party Data, TPA:assembly | — | Third Party Data |
| `TPA-WGS` | TPA-WGS | — | TPA, Third Party Data, TPA:assembly, WGS, STANDARD_DRAFT | — | TPA ドラフトゲノム |
| `TPA-GNM` | — | — | TPA, Third Party Data, TPA:assembly | — | TPA コンプリートゲノム |

**source_identifier** 列が空欄のカテゴリでは、ff_definition（DEFINITION 行）の生物名修飾子は `SOURCE_IDENTIFIER` キーで個別に指定してください。
`MAG` / `MAG-WGS` は `isolate` が自動的に使用されます（`SOURCE_IDENTIFIER` が指定されている場合はそちらが優先されます）。

**WGS の品質区分 KEYWORD** について: `STANDARD_DRAFT` はデフォルト値です。以下のうち1つを明示的に指定することもできます:

| KEYWORD | 説明 |
|---|---|
| `STANDARD_DRAFT` | 標準的なドラフトゲノム |
| `HIGH_QUALITY_DRAFT` | 高品質ドラフト |
| `IMPROVED_HIGH_QUALITY_DRAFT` | 改善された高品質ドラフト |
| `ANNOTATION_GRADE` | アノテーション向け品質 |
| `NON_CONTIGUOUS_FINISHED` | 未連結フィニッシュ |

**`MAG-WGS` の必須 source qualifier:**

| qualifier | 説明 |
|---|---|
| `isolate` | ゲノムの識別子（例: ゲノムビン名） |
| `metagenome_source` | 由来するメタゲノムの種類（例: `soil metagenome`） |
| `isolation_source` | サンプリング環境の説明 |
| `environmental_sample` | 自動付加（値なし） |

また `DBLINK` に `sequence read archive` が必要です。

---

## mss2ff の使い方

MSS アノテーションファイル（`.ann` または `.annt.tsv`）と FASTA ファイルから DDBJ Flat File を生成します。

### 基本的な使い方

```bash
# 基本
mss2ff annotation.ann genome.fa --division BCT --output output.ff

# 全オプション指定
mss2ff annotation.ann genome.fa \
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
| `FASTA` | — | なし | FASTA シーケンスファイル（CDS 翻訳に必要。source が COMMON に定義されている場合も必須） |
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
  - [Converting from .tbl / .fa directly](#converting-from-tbl--fa-directly)
  - [Common JSON File](#common-json-file)
  - [Sequence Role File (--sequence_roles)](#sequence-role-file---sequence_roles)
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
  - [Submission Categories (_submission_category)](#submission-categories-_submission_category)
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

The `gff2mss` subtool needs the extra `ddbj-gff` package (the other tools do not). Install with the extra only if you use gff2mss:

```bash
pip install -e ".[gff2mss]"
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
# Basic conversion (output: input.ann / input.fa)
egapx2mss input.asn --common examples/egapx2mss/common_example.json

# Specify output directory and filename
egapx2mss input.asn \
  --common examples/egapx2mss/common_example.json \
  --outdir results/ \
  --prefix output
```

By default, output files are written to the same directory as the input file, using the input basename as the prefix.
The second example above produces `results/output.ann` and `results/output.fa`.

### Options

| Option | Description |
|---|---|
| `input` | Input ASN.1 file (.asn). Can be omitted when both `--tbl` and `--fsa` are provided. |
| `-o`, `--outdir` | Output directory (created if absent; default: same directory as input file) |
| `-p`, `--prefix` | Output filename prefix, basename only — no directory separators (default: input basename) |
| `--common` | Common metadata JSON file (DBLINK, SUBMITTER, REFERENCE, etc.) |
| `--submission_category` | Submission category (`WGS`, `GNM`, `MAG-WGS`, etc.). Overrides `_submission_category` in the JSON (see below). |
| `--sequence_roles` | Sequence role file (TSV; see below). The legacy name `--chromosomes` is still accepted. |
| `--bin-dir` | Directory for asn2gb / asn2fsa binaries (default: `~/.local/share/ddbj_mss_tools/bin`) |
| `--keep-tmp` | Keep intermediate files (.tbl, raw FASTA) |
| `--tbl` | Pre-existing NCBI feature table (.tbl); skips step 1/3 |
| `--fsa` | Pre-existing FASTA file (.fa/.fsa); skips step 2/3 |
| `--preconvert-only` | Run steps 1/3 and 2/3 only (generate .tbl and .fa, then stop) |

### Converting from .tbl / .fa directly

When both `--tbl` and `--fsa` are provided, the ASN.1 input file can be omitted.
Only step 3/3 (MSS annotation conversion) is executed, and `asn2gb` / `asn2fsa` are not needed.

```bash
# Run step 3/3 only using existing .tbl and .fa files
egapx2mss --tbl input.tbl --fsa input.fa \
  --common common_example.json \
  --outdir results/ --prefix output
```

You can also run only steps 1/3 and 2/3 using `--preconvert-only`:

```bash
# Generate .tbl and .fa without MSS conversion
egapx2mss input.asn --preconvert-only --outdir tmp/
```

When re-running on the same input, any step whose output file already exists is automatically skipped.

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
    "SOURCE_IDENTIFIER": "cultivar",
    "ASSEMBLY_GAP": [
        {
            "enabled": true,
            "linkage_evidence": "proximity ligation",
            "min_gap_length": 100,
            "max_gap_length": 100,
            "gap_type": "within scaffold",
            "estimated_length": "unknown"
        },
        {
            "enabled": true,
            "linkage_evidence": "paired-ends",
            "min_gap_length": 10,
            "gap_type": "within scaffold",
            "estimated_length": "known"
        }
    ]
}
```

- Trailing commas (JSON5-style) are accepted.
- Sample file: [examples/egapx2mss/common_example.json](examples/egapx2mss/common_example.json)

#### _submission_category

Setting the `_submission_category` key automatically injects the appropriate DATATYPE, DIVISION, and KEYWORD values, and warns when required source qualifiers or DBLINK fields are missing.

```json
{
    "_submission_category": "GNM",
    ...
}
```

The `--submission_category` command-line option overrides the JSON value.
See [Submission Categories](#submission-categories-_submission_category) for a full list of supported categories.

#### SOURCE Section

Qualifiers listed under `SOURCE` are written directly into the source feature.

`SOURCE_IDENTIFIER` specifies which qualifier in `SOURCE` is used as the **intraspecific identifier** for the organism.
It is reflected in the DEFINITION line of the published flat file — for example, `"SOURCE_IDENTIFIER": "cultivar"` produces a definition like
`Brassica rapa NAPPA DNA, chromosome 1, complete sequence.`
The legacy key `INFRASPECIFIC_NAME_MODIFIER` is still accepted for backward compatibility.

#### ASSEMBLY_GAP Section

Runs of consecutive N bases are automatically detected and written as `assembly_gap` features in the annotation file.
`ASSEMBLY_GAP` is specified as an **array**, allowing multiple rules with different length ranges to be applied in priority order.
For each N-run, rules are evaluated from the first entry; the first matching rule is applied.
Output rows are in ascending coordinate order.

Fields for each rule:

| Field | Required | Default | Description |
|---|---|---|---|
| `enabled` | — | `true` | Set to `false` to skip this rule |
| `linkage_evidence` | ✓ | — | Evidence for how the gap was introduced (see table below) |
| `min_gap_length` | — | `10` | Minimum gap length to annotate |
| `max_gap_length` | — | no limit | Maximum gap length to annotate |
| `gap_type` | — | recommended | Value for the `gap_type` qualifier |
| `estimated_length` | — | recommended | Value for the `estimated_length` qualifier (`known` or `unknown`) |

If `gap_type` and `estimated_length` are omitted, recommended values are applied based on `linkage_evidence`:

| linkage_evidence | gap_type | estimated_length |
|---|---|---|
| `paired-ends` | within scaffold | known |
| `proximity ligation` | within scaffold | unknown |
| `align genus` | within scaffold | unknown |

Valid values for `linkage_evidence`:

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

**Example (Hi-C scaffolding + paired-ends, two rules):**

```json
"ASSEMBLY_GAP": [
    {
        "enabled": true,
        "linkage_evidence": "proximity ligation",
        "min_gap_length": 100,
        "max_gap_length": 100,
        "gap_type": "within scaffold",
        "estimated_length": "unknown"
    },
    {
        "enabled": true,
        "linkage_evidence": "paired-ends",
        "min_gap_length": 10,
        "gap_type": "within scaffold",
        "estimated_length": "known"
    }
]
```

In this example, N-runs of exactly 100 bases are annotated with the first rule (`unknown`),
while N-runs of 10–99 bases are annotated with the second rule (`known`).
N-runs shorter than 10 bases match no rule and are not annotated.

### Sequence Role File (--sequence_roles)

A 5-column tab-separated file (conventionally `sequence_roles.tsv`) that classifies sequences into chromosomes, organelles, plasmids, segments, or unplaced scaffolds.
Formerly called the *chromosome table* and passed via `--chromosomes`; the legacy option name is still accepted for backward compatibility.

```
# seq_id    type          seq_name       status    topology
Chr1        chromosome    1              complete  linear
ChrM        organelle     mitochondrion  complete  circular
scaffold001 unplaced                     partial   linear
```

| Column | Description |
|---|---|
| seq_id | Sequence ID from the FASTA header |
| type | `chromosome` / `organelle` / `plasmid` / `segment` / `unplaced` |
| seq_name | Chromosome number or organelle name. Required for `plasmid`, and for `chromosome` / `segment` when there are 2 or more entries of that type across the submission (an empty value raises an error). May be empty for a single chromosome, a single segment, or unplaced |
| status | `complete` / `partial` |
| topology | `linear` / `circular` |

- If omitted, all sequences are treated as unplaced and output in WGS mode.
- Lines beginning with `#` are treated as comments.
- The source feature's ff_definition (the DEFINITION line) is emitted as a DDBJ MSS meta-notation template (`@@[qualifier_name]@@`), which is expanded to the actual qualifier value on the same source feature at registration time.
- `segment` (a segment of a segmented/multipartite genome) is output as `complete genome` / `partial genome` (no `/segment` qualifier on source) when there is only one `segment` entry across the whole submission, or as `segment @@[segment]@@, complete sequence` / `segment @@[segment]@@, unlocalized sequence @@[entry]@@` (with `/segment` on source) when there are multiple.
- The full ff_definition decision table by type / count / status is below. `{P}` = `@@[organism]@@ @@[{source_identifier}]@@` (`source_identifier` = the name of the `SOURCE_IDENTIFIER` qualifier; `@@[organism]@@` alone when empty/None); `{mol}` is the mol_type-derived token (DNA / RNA / tRNA / rRNA / mRNA):

| type | count | status | ff_definition |
|------|-------|--------|---------------|
| unplaced (entry=None) | — | is_wgs=true | `{P} {mol}, @@[submitter_seqid]@@` |
| unplaced (entry=None) | — | is_wgs=false | `{P} {mol}, unplaced sequence @@[entry]@@` |
| chromosome | count==1 | complete | `{P} {mol}, chromosome, complete genome` |
| chromosome | count≥2 | complete | `{P} {mol}, chromosome @@[chromosome]@@, complete sequence` |
| chromosome | count==1 | partial | `{P} {mol}, chromosome, partial genome` |
| chromosome | count≥2 | partial | `{P} {mol}, chromosome @@[chromosome]@@, unlocalized sequence @@[entry]@@` |
| organelle | — | complete | `{P} {organelle_code} {mol}, complete genome` |
| organelle | — | partial | `{P} {organelle_code} {mol}, partial genome` |
| plasmid | — | complete | `{P} plasmid @@[plasmid]@@ {mol}, complete sequence` |
| plasmid | — | partial | `{P} plasmid @@[plasmid]@@ {mol}, partial sequence` |
| segment | count==1 | complete | `{P} {mol}, complete genome` |
| segment | count==1 | partial | `{P} {mol}, partial genome` |
| segment | count≥2 | complete | `{P} {mol}, segment @@[segment]@@, complete sequence` |
| segment | count≥2 | partial | `{P} {mol}, segment @@[segment]@@, unlocalized sequence @@[entry]@@` |
| other (unknown type) | — | — | `{P} {mol}, @@[entry]@@` |

  - The branches that emit `@@[chromosome]@@` / `@@[plasmid]@@` / `@@[segment]@@` (chromosome count≥2 / plasmid / segment count≥2) require a non-empty `seq_name`; an empty one raises `ValueError`.
  - `{organelle_code}` is the `/organelle` value converted to its DEFINITION adjectival form (`mitochondrion`→`mitochondrial`, `plastid:chloroplast`→`chloroplast`, etc.); it is not meta-notation.

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
mss_builder genome.fa --common examples/mss_builder/common_example.json

# Specify output directory and filename
mss_builder genome.fa \
  --common examples/mss_builder/common_example.json \
  --outdir results/ --prefix output

# Complete genome submission (with a sequence role file)
mss_builder genome.fa \
  --common examples/mss_builder/common_example.json \
  --sequence_roles sequence_roles.tsv \
  --outdir results/ --prefix output
```

### Options

| Option | Description |
|---|---|
| `input` | Input FASTA file (.fa / .fasta) |
| `-o`, `--outdir` | Output directory (created if absent; default: same directory as input file) |
| `-p`, `--prefix` | Output filename prefix, basename only — no directory separators (default: input basename) |
| `--common` | Common metadata JSON file (same format as egapx2mss) |
| `--submission_category` | Submission category (`WGS`, `GNM`, `MAG-WGS`, etc.). Overrides `_submission_category` in the JSON (see below). |
| `--sequence_roles` | Sequence role file TSV (if omitted, WGS mode is used). The legacy name `--chromosomes` is still accepted. |

The common JSON format is the [same as egapx2mss](#common-json-file).
`DBLINK.project` and `DBLINK.sample` are required.

### WGS Mode and Chromosome Mode

**WGS mode** (no `--sequence_roles`; the legacy `--chromosomes` option is equivalent):

- The source feature is written in the COMMON block using `@@[entry]@@` meta-notation.
- Each entry contains only `assembly_gap` features (if `ASSEMBLY_GAP` is configured).

**Chromosome mode** (`--sequence_roles` specified; the legacy `--chromosomes` option is equivalent):

- A separate source feature is written per entry.
- Chromosome/organelle names, topology, and completeness are derived from the sequence role file.

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
| `--submission_category` | — | Submission category (`WGS`, `MAG-WGS`, etc.). Overrides `_submission_category` in the JSON and in every TSV row. |

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
| `_` / `_submission_category` | Submission category (`WGS`, `MAG-WGS`, etc.; default is `WGS` if omitted; see below) |
| `DBLINK` / `project` | BioProject ID |
| `DBLINK` / `biosample` | BioSample ID |
| `DBLINK` / `sequence read archive` | DRA accession(s) (semicolon-separated for multiple) |
| `ST_COMMENT` / `Assembly Method` etc. | Genome-Assembly-Data qualifiers |
| `source` / `organism` etc. | source feature qualifiers |
| `COMMENT` / `line` | COMMENT block content (semicolon-separated for multiple lines) |

Sample file: [examples/batch_wgs_builder/sample_list_WGS.tsv](examples/batch_wgs_builder/sample_list_WGS.tsv)

### Common JSON File

Describes metadata common to all samples: SUBMITTER, REFERENCE, ASSEMBLY_GAP, SOURCE_IDENTIFIER, etc.
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
    "SOURCE_IDENTIFIER": "strain",
    "ASSEMBLY_GAP": [
        {
            "enabled": true,
            "linkage_evidence": "paired-ends",
            "min_gap_length": 10,
            "gap_type": "within scaffold",
            "estimated_length": "known"
        }
    ]
}
```

For details on `ASSEMBLY_GAP`, see the [ASSEMBLY_GAP Section](#assembly_gap-section) in the egapx2mss documentation.

Sample file: [examples/batch_wgs_builder/common_example.json](examples/batch_wgs_builder/common_example.json)

### Submission Categories (_submission_category)

The submission category can be specified in three ways (in decreasing priority):

1. **CLI option** `--submission_category CATEGORY` — overrides the TSV and JSON for all samples
2. **TSV** `_submission_category` column — set per row (`batch_wgs_builder` only)
3. **Common JSON** `_submission_category` key — default applied to all samples

When a category is set, DATATYPE / DIVISION / KEYWORD are injected automatically.
Missing required fields trigger a warning and are filled with empty values.

| Category | DATATYPE | DIVISION | Main KEYWORD | source_identifier | Typical use |
|---|---|---|---|---|---|
| (none) | — | — | — | — | No category |
| `GNM` | — | — | — | — | Complete genome |
| `WGS` | WGS | — | WGS, STANDARD_DRAFT | — | Draft genome (WGS) |
| `ENV` | — | ENV | ENV | — | Environmental DNA |
| `MAG` | — | ENV | ENV, MAG, Metagenome Assembled Genome | `isolate` | MAG (complete genome) |
| `MAG-WGS` | WGS | ENV | ENV, MAG, Metagenome Assembled Genome, WGS, STANDARD_DRAFT | `isolate` | MAG (draft genome) |
| `TSA` | — | TSA | TSA, Transcriptome Shotgun Assembly | — | Transcriptome Shotgun Assembly |
| `TPA` | TPA | — | TPA, Third Party Data, TPA:assembly | — | Third Party Data |
| `TPA-WGS` | TPA-WGS | — | TPA, Third Party Data, TPA:assembly, WGS, STANDARD_DRAFT | — | TPA draft genome |
| `TPA-GNM` | — | — | TPA, Third Party Data, TPA:assembly | — | TPA complete genome |

The **source_identifier** column indicates the qualifier automatically used as the intraspecific modifier in the DEFINITION line.
For `MAG` / `MAG-WGS`, `isolate` is used by default; an explicit `SOURCE_IDENTIFIER` in the JSON takes precedence.

**WGS draft quality keywords:** `STANDARD_DRAFT` is the default. One of the following may be specified explicitly:

| KEYWORD | Description |
|---|---|
| `STANDARD_DRAFT` | Standard draft |
| `HIGH_QUALITY_DRAFT` | High-quality draft |
| `IMPROVED_HIGH_QUALITY_DRAFT` | Improved high-quality draft |
| `ANNOTATION_GRADE` | Annotation-grade quality |
| `NON_CONTIGUOUS_FINISHED` | Non-contiguous finished |

**Required source qualifiers for `MAG-WGS`:**

| Qualifier | Description |
|---|---|
| `isolate` | Genome identifier (e.g. bin name) |
| `metagenome_source` | Type of source metagenome (e.g. `soil metagenome`) |
| `isolation_source` | Description of the sampling environment |
| `environmental_sample` | Added automatically (no value) |

`DBLINK` must also include a `sequence read archive` entry.

---

## mss2ff Usage

Generates a DDBJ Flat File from an MSS annotation file (`.ann` or `.annt.tsv`) and a FASTA file.

### Basic Usage

```bash
# Basic
mss2ff annotation.ann genome.fa --division BCT --output output.ff

# All options
mss2ff annotation.ann genome.fa \
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
| `FASTA` | — | none | FASTA sequence file (required for CDS translation; also required when source is defined only in COMMON) |
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
