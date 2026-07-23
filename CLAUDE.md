# CLAUDE.md

## プロジェクト概要

DDBJ MSS (Mass Submission System) 登録ファイルを生成するPythonツール群。

- **egapx2mss**: NCBI EGAPx の出力 (ASN.1形式) を DDBJ MSS形式 (.ann / .fa) に変換。既存の .tbl / .fa ファイルから直接変換することも可能
- **batch_wgs_builder**: FASTAファイルと TSV のメタデータから DDBJ MSS WGS/MAG-WGS形式ファイルを一括生成
- **mss2ff**: MSS アノテーションファイル (.ann / .annt.tsv) と FASTA から DDBJ Flat File を生成
- **mss_builder**: MSS 登録ファイルのビルド支援ツール

## リポジトリ構成

```
src/
├── common/
│   ├── cli_args.py          # 共通CLIオプション (-o/--outdir, -p/--prefix) のargparseヘルパー
│   ├── common_builder.py    # MSS COMMON行生成ユーティリティ（各ツール共通）
│   ├── models.py            # 共通 pydantic モデル
│   ├── fasta.py             # FASTA読み書きユーティリティ
│   ├── gap_annotator.py     # assembly_gap フィーチャー生成
│   └── source_builder.py   # source フィーチャー生成
├── egapx2mss/
│   ├── cli.py               # エントリーポイント (egapx2mss コマンド)
│   ├── asn_tools.py         # asn2gb/asn2fsa のダウンロード・実行・期限切れリトライ
│   ├── tbl_parser.py        # NCBI feature table (.tbl) パーサー・ロケーション変換
│   └── ann_writer.py        # DDBJ MSS アノテーションファイル (.ann) 書き出し
├── batch_wgs_builder/
│   ├── cli.py               # エントリーポイント (batch_wgs_builder コマンド)
│   └── core.py              # メインロジック (TSV → MSS変換)
├── mss_builder/
│   ├── cli.py               # エントリーポイント (mss_builder コマンド)
│   └── ann_writer.py        # アノテーションファイル書き出し
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
├── batch_wgs_builder/       # sample_list.xlsx, *.fna.gz
└── mss2ff/                  # DDBJ.annt.tsv, DDBJ.seq.fa (入力例)

deprecated/                  # 移行前の旧コード (egapx2mss/, wgs_maker/)
tests/                       # pytest テスト (今後追加予定)
```

新規開発は `src/` 以下で行う。`deprecated/` は参照のみ。

## セットアップ

```bash
# 通常インストール
pip install -e .

# Docker
docker build -t ddbj-mss-tools .
docker run --rm -it -v $(pwd):/app ddbj-mss-tools

# 軽量イメージ (Dockerfile.slim, 約329MB / 通常版は約719MB)
# build-essential と examples/ を含めず、NCBI バイナリ (asn2gb/asn2fsa) は /app/bin に焼き込み済み
# gff2mss 用の ddbj-gff wheel を ../gff_submission から先に生成する（コミットされない。ddbj-gff 更新時は再実行）
scripts/build-ddbj-gff-wheel.sh
docker build -t ddbj-mss-tools:slim -f Dockerfile.slim .
# 注意: egapx2mss を使う場合、焼き込み済みバイナリ (/app/bin) を隠さないよう
#       作業ディレクトリは /app 以外にマウントする
docker run --rm -it -v $(pwd):/data -w /data ddbj-mss-tools:slim
```

## 依存パッケージ

`requirements.txt` で管理。`pyproject.toml` の `dependencies`（コア）と内容を一致させること。

```
pydantic>=2.0, biopython, pandas, openpyxl, jsonschema
```

`gff2mss` サブツールだけが必要とする `ddbj-gff` は**コア依存ではなく optional extra**（`pyproject.toml` の `[project.optional-dependencies]` の `gff2mss`）。他ツール（egapx2mss / mss_builder / mss2ff / batch_wgs_builder）は `ddbj-gff` 無しで install・実行できる。gff2mss を使う場合のみ `pip install ".[gff2mss]"`（または別途 `ddbj-gff` を導入）。`ddbj-gff` は PyPI 非公開のため、コンテナ/CI では wheel を導入する。`gff2mss` 実行時に未導入なら親切なエラーで終了する（`src/gff2mss/cli.py`。パッケージ import 時に `ddbj-gff` を要求しないよう `src/gff2mss/__init__.py` は遅延 re-export）。


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

### source の flag 型 qualifier（値なし qualifier）

`environmental_sample` など値を持たない qualifier（INSDC flag 型）は、真偽値で付与を制御できる。

- 付与しない: `false` / `no`（大文字小文字問わず）、JSON boolean `false`、または **TSV の空欄**
- 付与する（値なし行 `/qualifier` として出力）: 上記以外。**推奨は `yes`**。JSON では `true`（boolean）や空文字 `""`（後方互換）も付与扱い。

対象 flag: `environmental_sample`, `transgenic`, `germline`, `rearranged`, `proviral`, `macronuclear`, `metagenomic`, `focus`。
flag 型でない通常 qualifier（`strain` 等）はこの判定の対象外で、値はそのまま出力される（例: `strain` の値 `No` は文字列 `No` のまま）。

注意: ENV/MAG/MAG-WGS など `environmental_sample` を自動付与するカテゴリでは、`false` を指定しても自動付与が優先される。

### sequence role (`--sequence_roles` TSV) と ff_definition

5列 TSV: `seq_id <TAB> type <TAB> seq_name <TAB> status <TAB> topology`。`type` は
`chromosome` / `organelle` / `plasmid` / `segment` / `unplaced` のいずれか。`status` は `complete` / `partial`。
type に応じて source の ff_definition（DDBJ Flat File の DEFINITION 行）が下記のように構築される。

ff_definition は DDBJ MSS のメタ記法（`@@[qualifier_name]@@`）を使ったテンプレート文字列として出力され、
MSS 側の登録処理で同じ source フィーチャーが持つ qualifier の実値に展開される。

- **`{P}`（prefix）**: `source_identifier`（`SOURCE_IDENTIFIER` の qualifier 名。例 `cultivar`）が
  非空なら `@@[organism]@@ @@[{source_identifier}]@@`、空/None なら `@@[organism]@@` のみ。
- **`{mol}`**: 従来どおり mol_type 由来の具体値（DNA/RNA/tRNA/rRNA/mRNA）。メタ化しない。

| type | count | status | ff_definition |
|------|-------|--------|---------------|
| unplaced（entry=None） | — | is_wgs=true | `{P} {mol}, @@[submitter_seqid]@@` |
| unplaced（entry=None） | — | is_wgs=false | `{P} {mol}, unplaced sequence @@[entry]@@` |
| chromosome | count==1 | complete | `{P} {mol}, chromosome, complete genome` |
| chromosome | count==1 | partial | `{P} {mol}, chromosome` |
| chromosome | count≥2 | complete | `{P} {mol}, chromosome @@[chromosome]@@, complete sequence` |
| chromosome | count≥2 | partial | `{P} {mol}, chromosome @@[chromosome]@@` |
| organelle | — | complete | `{P} {organelle_code} {mol}, complete genome` |
| organelle | — | partial | `{P} {organelle_code} {mol}, partial genome` |
| plasmid | — | complete | `{P} plasmid @@[plasmid]@@ {mol}, complete sequence` |
| plasmid | — | partial | `{P} plasmid @@[plasmid]@@ {mol}, partial sequence` |
| segment | count==1 | complete | `{P} {mol}, complete genome` |
| segment | count==1 | partial | `{P} {mol}, partial genome` |
| segment | count≥2 | complete | `{P} {mol}, segment @@[segment]@@, complete sequence` |
| segment | count≥2 | partial | `{P} {mol}, segment @@[segment]@@` |
| その他（未知 type） | — | — | `{P} {mol}, @@[entry]@@` |

`@@[chromosome]@@` / `@@[plasmid]@@` / `@@[segment]@@` を出力する分岐（chromosome count≥2、
plasmid、segment count≥2）では entry の `seq_name` が必須で、空文字の場合は `ff_definition()` が
`ValueError` を送出する。単一 chromosome（count==1）・単一 segment（count==1）は `seq_name` を
参照しないため空文字でも許容される。

organelle の `{organelle_code}` は INSDC `/organelle` 値（`mitochondrion`, `plastid:chloroplast` 等）を
DEFINITION 用の形容詞形（`mitochondrial`, `chloroplast` 等）に変換した値で、こちらはメタ記法化せず
従来どおり変換済みの具体値を出力する（`@@[organelle]@@` にすると生値 `mitochondrion` が展開され
不整合になるため）。変換表に無い値はそのまま使う。source の `/organelle` qualifier には変換前の
生の値が出力される。

## mss_builder コマンドオプション

```
mss_builder input.fa [オプション]
```

| オプション | 説明 |
|-----------|------|
| `input` | 入力 FASTA ファイル (.fa / .fasta) |
| `-o/--outdir DIR` | 出力先ディレクトリ（存在しない場合は自動作成）。デフォルト: 入力ファイルと同じディレクトリ |
| `-p/--prefix NAME` | 出力ファイルのベースネーム（ディレクトリ区切り文字不可）。デフォルト: 入力ファイルのベースネーム |
| `--common JSON` | 共通メタデータ JSON ファイル (DBLINK, SUBMITTER, REFERENCE, ASSEMBLY_GAP 等) |
| `--sequence_roles TSV` | seq_id → 染色体/organelle マッピング sequence role ファイル (5列 TSV)。旧名 `--chromosomes` も互換のため受け付ける |

**典型的な使い方:**
```bash
# 通常変換 (出力は input.ann / input.fa)
mss_builder input.fa

# 出力ディレクトリとファイル名を指定
mss_builder input.fa -o results/ -p submission --common common.json
```

## egapx2mss コマンドオプション

```
egapx2mss [input.asn] [オプション]
```

| オプション | 説明 |
|-----------|------|
| `input` | 入力 ASN.1 ファイル (.asn)。`--tbl` と `--fsa` を両方指定する場合は省略可 |
| `-o/--outdir DIR` | 出力先ディレクトリ（存在しない場合は自動作成）。デフォルト: 入力ファイルと同じディレクトリ |
| `-p/--prefix NAME` | 出力ファイルのベースネーム（ディレクトリ区切り文字不可）。デフォルト: 入力ファイルのベースネーム |
| `--tbl FILE` | 既存の NCBI feature table (.tbl) を直接指定。step 1/3 をスキップ |
| `--fsa FILE` | 既存の FASTA ファイル (.fa/.fsa) を直接指定。step 2/3 をスキップ |
| `--common JSON` | 共通メタデータ JSON ファイル (DBLINK, SUBMITTER, REFERENCE 等) |
| `--sequence_roles TSV` | seq_id → 染色体/organelle マッピング sequence role ファイル (5列 TSV)。旧名 `--chromosomes` も互換のため受け付ける |
| `--keep-tmp` | 中間ファイル (.tbl, _raw.fa) を削除せずに保持 |
| `--preconvert-only` | step 1/3・2/3 のみ実行して終了 (.tbl と .fa を生成) |
| `--bin-dir DIR` | asn2gb/asn2fsa バイナリの置き場所 |

**`--tbl` と `--fsa` の組み合わせルール:**
- 両方同時に指定する必要がある（片方のみは不可）
- `input` (.asn) との同時指定は不可
- 両方指定時は asn2gb/asn2fsa のダウンロード・実行をスキップ
- `-o/--outdir` と `-p/--prefix` を省略した場合は `--tbl` ファイルと同じディレクトリ・ベースネームをデフォルトとする

**典型的な使い方:**
```bash
# 通常変換 (出力は brapa.ann / brapa.fa)
egapx2mss brapa.asn

# 出力ディレクトリとファイル名を指定
egapx2mss brapa.asn -o results/ -p output

# step 1/2 だけ実行して .tbl と .fa を生成
egapx2mss brapa.asn --preconvert-only -o tmp/

# 既存の .tbl / .fa から step 3 だけ実行
egapx2mss --tbl brapa.tbl --fsa brapa.fa -o results/ -p output
```

## テスト

```bash
pytest
```

テストは `tests/` に追加する。`examples/` のファイルを入力として使う統合テストを優先。

## 開発上の留意点

- Python 3.10 以上必須 (`match` 文や `X | Y` 型ヒントを使用)
- `src/common/` は各ツールから使う共有ライブラリ
- 新しいバリデーションルールは `src/common/models.py` の pydantic モデルに追加する
- `requirements.txt` を更新したら `pyproject.toml` の `dependencies` も合わせて更新する
- `-o/--outdir` と `-p/--prefix` オプションは `src/common/cli_args.py` の共通ヘルパー (`add_output_args`, `validate_prefix`, `resolve_output`) を使う。新ツールで同パターンを追加する場合はここから import すること
