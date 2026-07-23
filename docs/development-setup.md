# 開発環境の構築（別マシンへの移行手順）

`ddbj_mss_tools` の開発環境を新しいマシンで再構築するための手順。

## 前提

| 項目 | 要件 |
|---|---|
| Python | **3.11 以上**（`ddbj_mss_tools` 自体は 3.10+ だが、`gff2mss` が使う `ddbj-gff` が 3.11+ を要求するため、揃えて 3.11+ にする） |
| git | 必須 |
| ネットワーク | 依存パッケージの取得、および egapx2mss が使う NCBI バイナリの初回自動ダウンロードに必要 |
| Docker（任意） | slim イメージをビルドする場合のみ |
| uv（任意・推奨） | `gff_submission` が uv 管理（`uv.lock`）。gff2mss を開発するなら uv があると楽 |

## リポジトリ配置（重要: 隣接させる）

`gff2mss` サブツールは外部ライブラリ `ddbj-gff`（別リポジトリ `gff_submission`）に依存し、開発時は
`ddbj_mss_tools/pyproject.toml` の `[tool.uv.sources]` が **`../gff_submission`** を editable 参照します。
そのため2リポジトリを**同じ親ディレクトリに並べて**チェックアウトします。

```
<任意の親ディレクトリ>/
├── ddbj_mss_tools/     # このリポジトリ（GitHub: ddbj/ddbj_mss_tools）
└── gff_submission/     # ddbj-gff 本体（gff2mss を触る場合のみ必要）
```

- `ddbj_mss_tools`: `git clone git@github.com:ddbj/ddbj_mss_tools.git`
- `gff_submission`: **現状 git remote が無い**。旧マシンからディレクトリごとコピーするか、リモートを設定して転送する（`gff_submission/docs/development-setup.md` 参照）。`gff2mss` を使わないなら不要。

## セットアップ

### A. 4ツール（egapx2mss / mss_builder / mss2ff / batch_wgs_builder）だけ使う場合

`ddbj-gff` は不要。仮想環境を作って editable install するだけ:

```bash
cd ddbj_mss_tools
python3.11 -m venv .venv && source .venv/bin/activate   # または conda create -n mss python=3.11
pip install -e .          # コア依存のみ（pydantic, biopython, pandas, openpyxl, jsonschema）
pip install pytest        # テスト用
```

`pip install -e .` で `ddbj-gff` は入りません（optional extra 化済み）。これら4ツールは `ddbj-gff` 無しで動作します。

### B. gff2mss も開発する場合

上記に加えて、隣接する `gff_submission` を **editable** で入れます（これが `ddbj_gff` を提供）:

```bash
# 上の A を実施済みの環境で
pip install -e ../gff_submission      # ddbj_gff を editable 提供（gff_submission を編集すれば即反映）
```

> 注意: `pip install -e ".[gff2mss]"` は pip が `[tool.uv.sources]` を解釈せず `ddbj-gff` を PyPI に探しに行って失敗します（PyPI 未公開）。pip では上記のように `../gff_submission` を直接 editable install してください。uv を使う場合は `uv pip install -e .` で `[tool.uv.sources]` の path 参照が効きます。

## テストの実行

```bash
# editable install 済みなら
pytest                                   # tests/ を実行

# ソースから直接動かす場合（install せず）
PYTHONPATH=src pytest -m "not slow"
```

- `gff2mss` 関連テストは `ddbj_gff` が必要です。B のセットアップをしていないと、それらは collection error になります（他ツールのテストは通ります）。
- `-m "not slow"` で大きな example を使う統合テストを除外できます。

## egapx2mss の NCBI バイナリ

`egapx2mss` は `asn2gb` / `asn2fsa`（NCBI）を**初回実行時に自動ダウンロード**して `ddbj_mss_tools/bin/`（既定）にキャッシュします。手動導入は不要ですが、初回はネットワークが必要です。`--bin-dir` で保存先を変更可。

## Docker（任意）

```bash
# 通常イメージ
docker build -t ddbj-mss-tools .

# 軽量イメージ（gff2mss も動く）: ddbj-gff wheel を ../gff_submission から先に生成してから build
scripts/build-ddbj-gff-wheel.sh
docker build -t ddbj-mss-tools:slim -f Dockerfile.slim .
```

生成した wheel（`ddbj_gff-*.whl`）はコミットされません（`.gitignore` 済み）。`ddbj-gff` を更新したら
スクリプトを再実行してから build し直します。詳細は `gff_submission/docs/mss-tools-integration.md`。

## 関連リポジトリ

- **gff_submission** — `ddbj-gff`（GFF 正規化・パースライブラリ）。gff2mss が依存。セットアップは同リポジトリの `docs/development-setup.md`。
- **mss_tools_web**（別リポジトリ）— これらの CLI ツールをコンテナで呼び出す Web アプリ。本リポジトリの開発には必須ではない。
