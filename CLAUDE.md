# CLAUDE.md

## プロジェクト概要

DDBJ MSS (Mass Submission System) 登録ファイルを生成するPythonツール群。
2つのコマンドラインツールを含む。

- **egapx2mss**: NCBI EGAPx の出力 (ASN.1形式) を DDBJ MSS形式 (.ann / .fa) に変換
- **wgs_maker**: FASTAファイルと Excel/TSV のメタデータから DDBJ MSS形式ファイルを一括生成

## リポジトリ構成

```
src/
├── common/
│   └── json2mss.py          # MSS COMMON行生成ユーティリティ（両ツール共通）
├── egapx2mss/
│   ├── cli.py               # エントリーポイント (egapx2mss コマンド)
│   ├── models.py            # pydantic モデル・load_common_json
│   ├── asn_tools.py         # asn2gb/asn2fsa のダウンロード・実行・期限切れリトライ
│   ├── fasta.py             # FASTA読み書きユーティリティ
│   ├── tbl_parser.py        # NCBI feature table (.tbl) パーサー・ロケーション変換
│   └── ann_writer.py        # DDBJ MSS アノテーションファイル (.ann) 書き出し
└── wgs_maker/
    ├── cli.py               # エントリーポイント (wgs_maker コマンド)
    ├── core.py              # メインロジック (Excel/TSV → MSS変換)
    ├── gap_annotator.py     # assembly_gap フィーチャー生成
    ├── schema_util.py       # JSON Schema バリデーション
    └── seq_util.py          # FASTA読み込み・source フィーチャー生成

examples/
├── egapx2mss/               # brapa.asn, common_example.json
└── wgs_maker/               # sample_list_*.tsv, sample_list.xlsx, common_example.json

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
# egapx2mss
egapx2mss input.asn --organism "Brassica rapa" --common examples/egapx2mss/common_example.json

# wgs_maker
wgs_maker --tsv examples/wgs_maker/sample_list_WGS.tsv -m examples/wgs_maker/common_example.json -o OUT
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
