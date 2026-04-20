# CLAUDE.md

## プロジェクト概要

DDBJ MSS (Mass Submission System) 登録ファイルを生成するPythonツール群。
2つのコマンドラインツールを含む。

- **egapx2mss**: NCBI EGAPx の出力 (ASN.1形式) を DDBJ MSS形式 (.ann / .fa) に変換
- **wgs_maker**: FASTAファイルと Excel/TSV のメタデータから DDBJ MSS形式ファイルを一括生成
- **mss2ff**: MSS アノテーションファイル (.ann / .annt.tsv) と FASTA から DDBJ Flat File を生成
- **mss_builder**: MSS 登録ファイルのビルド支援ツール

## リポジトリ構成

```
src/
├── common/
│   ├── json2mss.py          # MSS COMMON行生成ユーティリティ（両ツール共通）
│   ├── models.py            # 共通 pydantic モデル
│   ├── fasta.py             # FASTA読み書きユーティリティ
│   ├── gap_annotator.py     # assembly_gap フィーチャー生成
│   └── source_builder.py   # source フィーチャー生成
├── egapx2mss/
│   ├── cli.py               # エントリーポイント (egapx2mss コマンド)
│   ├── models.py            # pydantic モデル・load_common_json
│   ├── asn_tools.py         # asn2gb/asn2fsa のダウンロード・実行・期限切れリトライ
│   ├── tbl_parser.py        # NCBI feature table (.tbl) パーサー・ロケーション変換
│   └── ann_writer.py        # DDBJ MSS アノテーションファイル (.ann) 書き出し
├── wgs_maker/
│   ├── cli.py               # エントリーポイント (wgs_maker コマンド)
│   ├── core.py              # メインロジック (Excel/TSV → MSS変換)
│   └── schema_util.py       # JSON Schema バリデーション
├── mss_builder/
│   └── cli.py               # エントリーポイント (mss_builder コマンド)
└── mss2ff/
    ├── __init__.py
    ├── cli.py               # エントリーポイント (mss2ff コマンド)
    ├── ann_parser.py        # MSS .ann/.annt.tsv パーサー
    ├── ff_writer.py         # DDBJ Flat File ライター
    ├── location.py          # MSS ロケーション文字列パーサー・BioPython変換
    ├── taxonomy.py          # NCBI Entrez 分類情報取得・ST_COMMENT タグセットID変換
    └── translate_with_transl_except.py  # transl_except 対応 CDS翻訳

examples/
├── egapx2mss/               # brapa.asn, common_example.json
├── wgs_maker/               # sample_list_*.tsv, sample_list.xlsx, common_example.json
└── mss2ff/                  # DDBJ.annt.tsv, DDBJ.seq.fa (入力例)

tests/                       # pytest テスト (今後追加予定)
```

**旧ディレクトリ**: `egapx2mss/` と `wgs_maker/` は移行前の旧コードが残っている。新規開発は `src/` 以下で行う。

## セットアップ

```bash
# 通常インストール
pip install -e .

# Docker
docker build -t ddbj-mss-tools .
docker run --rm -it -v $(pwd):/app ddbj-mss-tools
```

## 依存パッケージ

`requirements.txt` で管理。`pyproject.toml` の `dependencies` と内容を一致させること。

```
pydantic>=2.0, biopython, pandas, openpyxl, jsonschema
```

## コマンド使用例

```bash
# egapx2mss (organism は --common の JSON 内 SOURCE.organism で指定)
egapx2mss input.asn --common examples/egapx2mss/common_example.json

# wgs_maker
wgs_maker --tsv examples/wgs_maker/sample_list_WGS.tsv -m examples/wgs_maker/common_example.json -o OUT

# mss2ff (基本)
mss2ff examples/mss2ff/DDBJ.annt.tsv --fasta examples/mss2ff/DDBJ.seq.fa --division BCT -o output.ff

# mss2ff (全オプション)
mss2ff DDBJ.annt.tsv \
    --fasta DDBJ.seq.fa \
    --output output.ff \
    --division BCT \
    --submission-date 2025-04-01 \
    --file-date 2025-04-01 \
    --email your@email.com \
    --accession AP000001 \
    --no-taxonomy
```

## DDBJ MSSファイル形式

アノテーションファイル (.ann) は **5列タブ区切り** テキスト。

```
列1: entry名 (最初のフィーチャーのみ記載、以降は空)
列2: フィーチャー名 (source, gene, mRNA, CDS, assembly_gap, ...)
列3: ロケーション (1..100, join(1..50,80..100), complement(...))
列4: qualifier名
列5: qualifier値
```

- ファイル先頭に `COMMON` ブロック（submitter, reference, DBLINK等）
- FASTA ファイル (.fsa) は各エントリの末尾に `//` セパレーター

## common_example.json の形式

`DBLINK.project` と `DBLINK.sample` が必須。pydantic でバリデーションされる。

```json
{
    "DBLINK": {
        "project": "PRJD000001",
        "sample": "SAMD000001",
        "DRA": ["DRA000001"]
    },
    "SUBMITTER": { "ab_name": ["Tanizawa,Y."], ... },
    "REFERENCE": [{ "title": "...", "ab_name": [...], "status": "Unpublished", "year": 2025 }]
}
```

末尾カンマ (JSON5スタイル) は許容される。

## egapx2mss 固有の注意点

### asn2gb / asn2fsa バイナリ
- NCBI ツール。プロジェクトルートの `bin/` にキャッシュ (`--bin-dir` で変更可)
- **利用期限あり**。出力が空のとき自動で最新版を再ダウンロードして1回リトライする
- ダウンロード元: `https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/`
- macOS: `asn2gb.mac.gz` / `asn2fsa.mac.gz`、Linux: `*.linux64.gz`

### ASN.1 ファイル
- 1ファイルに複数の `Seq-entry ::=` ブロックが連結されたフォーマット (catenated)
- `asn2gb` は `-a q` オプションで catenated ファイルに対応
- `asn2fsa` は catenated 非対応のため、ブロックごとに一時ファイルへ分割して実行する (`asn_tools.py`)

### ロケーション変換 (`tbl_parser.py`)
- NCBI tbl のマイナス鎖: `start > end` → DDBJ MSS: `complement(start..end)`
- 部分配列マーカー `<` (5'端) / `>` (3'端) を DDBJ 形式に変換する

## wgs_maker 固有の注意点

- 登録カテゴリ (`_trad_submission_category`): `WGS`, `MAG`, `GNM`, `MAG-WGS`
- `assembly_gap` フィーチャーは連続する N (デフォルト10塩基以上) を自動検出
- JSON Schema バリデーションは `src/wgs_maker/MSS_COMMON_template.json` を使用
  (ネットワーク接続時は GitHub から最新版を取得)

## mss2ff 固有の注意点

### CLIオプション

| オプション | 省略形 | デフォルト | 説明 |
|---|---|---|---|
| `ANN` | — | (必須) | MSS アノテーションファイル (.ann または .annt.tsv) |
| `--fasta` | `-f` | なし | FASTA/FSA シーケンスファイル (CDS翻訳に必要) |
| `--output` | `-o` | stdout | 出力ファイルパス |
| `--division` | `-d` | `UNK` | DDBJ division コード (BCT, VRL, PLN 等) |
| `--submission-date` | `-s` | 今日 | REFERENCE 1 の投稿日 (YYYY-MM-DD または DD-Mon-YYYY) |
| `--file-date` | — | 今日 | LOCUS 行のファイル作成日 |
| `--email` | — | `mss2ff@ddbj.nig.ac.jp` | NCBI Entrez API 用メールアドレス |
| `--accession` | `-a` | なし | 開始アクセッション番号 (後述) |
| `--no-taxonomy` | — | false | NCBI 分類情報の取得をスキップ |

### アクセッション番号フォーマット (`--accession`)

3種類のフォーマットをサポート。serial 部は6桁以上の可変長。

| フォーマット | 例 | 内訳 |
|---|---|---|
| 2文字 + serial | `AP000001` | プレフィックス2文字 + serial6桁以上 |
| 4文字 + 2桁バージョン + serial | `AAXJ010000001` | プレフィックス4文字 + バージョン2桁 + serial6桁以上 |
| 6文字 + 2桁バージョン + serial | `AAXJEM010000001` | プレフィックス6文字 + バージョン2桁 + serial6桁以上 |

エントリーはアノテーションファイルの順番に serial を +1 ずつ割り当てる。

### REFERENCE JOURNAL フォーマット

`status` 値により JOURNAL 行のフォーマットが変わる:

| status | JOURNAL 出力形式 |
|---|---|
| `Unpublished` (year なし) | `Unpublished.` |
| `Unpublished` (year あり) | `Unpublished. (year)` |
| `In press` | `{journal} ({year}) In press` |
| `Published` | `{journal} {volume}, {from_page}-{to_page} ({year})` |

### REFERENCE 1 (Direct Submission)

COMMON ブロックの SUBMITTER 情報から自動生成。JOURNAL 行の形式:

```
Submitted (DD-Mon-YYYY) to the DDBJ/EBI/GenBank databases.
Contact:contact_name
[department, ]institute;
street, city, state zip, country
URL    :url
```

- `department` が空の場合は省略
- `URL` フィールドは7文字分のアライメント (`URL    :`)

### CONSRTM (コンソーシアム著者)

REFERENCE および SUBMITTER ブロックで `consrtm` qualifier をサポート。
AUTHORS 行に個人著者名が存在しない場合は AUTHORS 行を省略し、CONSRTM 行のみ出力する。

### assembly_gap の estimated_length 展開

アノテーションファイルに `estimated_length=known` と記載されている場合、
ロケーション文字列から実際のギャップ長を計算して数値に置き換える。
(例: ロケーション `100..109` → `estimated_length=10`)

### ST_COMMENT タグセットID

`tagset_id` の値をそのままヘッダ/フッタに使用する:

| tagset_id | ヘッダ/フッタ |
|---|---|
| `Genome-Assembly-Data` | `##Genome-Assembly-Data-START##` / `##Genome-Assembly-Data-END##` |
| `Assembly-Data` | `##Assembly-Data-START##` / `##Assembly-Data-END##` |

### FASTA/FSA フォーマット

DDBJ FSA フォーマット (エントリ末尾に `//` セパレーター) を正しく処理する。
`//` 行はシーケンスデータとして読み込まず、スキップする。

### CDS 翻訳

- `pseudo` または `pseudogene` qualifier がある CDS は翻訳しない
- `transl_except` qualifier がある CDS は `translate_with_transl_except.py` で処理し、
  終止コドン位置をアミノ酸に変換してから `translation` qualifier を生成する

## テスト

```bash
pytest
```

テストは `tests/` に追加する。`examples/` のファイルを入力として使う統合テストを優先。

## 開発上の留意点

- Python 3.10 以上必須 (`match` 文や `X | Y` 型ヒントを使用)
- `src/common/` は `egapx2mss` と `wgs_maker` 両方から使う共有ライブラリ
- 新しいバリデーションルールは `src/egapx2mss/models.py` の pydantic モデルに追加する
- `requirements.txt` を更新したら `pyproject.toml` の `dependencies` も合わせて更新する
