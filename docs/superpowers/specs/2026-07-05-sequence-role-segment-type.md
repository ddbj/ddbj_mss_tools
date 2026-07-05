# sequence role `type=segment` 対応 設計

## 背景と目的

分節ゲノム（influenza などの multipartite ウイルス）を DDBJ MSS で登録する際、
エントリごとに INSDC `/segment` qualifier を付け、DEFINITION 行に segment を反映したい。

現状の sequence role ファイル (`--sequence_roles` TSV) の `type` 列は
`chromosome` / `organelle` / `plasmid` / `unplaced` の 4 種のみを扱い、segment を表現できない。
`type` にホワイトリスト検証は無く、未知 type は `source_qualifier()` / `ff_definition()` の
フォールバックに落ちるため、`segment` と書いても専用の qualifier / DEFINITION は生成されない。

本設計で `type=segment` を追加し、submission 内の segment 数と `status` に応じて
source qualifier と ff_definition を生成する。

対象ツール: **mss_builder** と **egapx2mss** の両方（実装は共有モジュール
`src/common/source_builder.py` に置き、両ツールの ann_writer から利用する）。

## 挙動仕様

判定軸は 2 つ:

- `segment_count` = submission 内の `type == "segment"` エントリ数
- `status` 列（`complete` / それ以外は partial 扱い）

これは既存の `chromosome`（`chromosome_count` による genome/sequence 切替）と
`organelle` / `plasmid`（status による complete/partial 切替）を組み合わせた挙動に相当する。

### source フィーチャーの `/segment` qualifier

| 条件 | `/segment` |
|---|---|
| segment が 1 つだけ (`segment_count == 1`) | 付与しない（`seq_name` は不要・あっても無視） |
| segment が複数 (`segment_count >= 2`) | `/segment = {seq_name}` を付与 |

`seq_name` が空文字の場合は既存の chromosome / plasmid と同様に qualifier を省略する。

### ff_definition

`{prefix}` = `{organism} {識別子の値}`（識別子は `SOURCE_IDENTIFIER` の qualifier 値。空なら organism のみ）、
`{mol}` = `mol_type` 由来トークン（DNA / RNA / tRNA / rRNA / mRNA。既存の `_molecule_token()`）。

| segment 数 | status | ff_definition |
|---|---|---|
| 単一 (`segment_count == 1`) | complete | `{prefix} {mol}, complete genome` |
| 単一 | partial 等 | `{prefix} {mol}, partial genome` |
| 複数 (`segment_count >= 2`) | complete | `{prefix} {mol}, segment {seq_name}, complete sequence` |
| 複数 | partial 等 | `{prefix} {mol}, segment {seq_name}, partial sequence` |

複数かつ `seq_name` が空の場合は segment 語のみ（`{prefix} {mol}, segment, complete sequence`）とする
（chromosome の `chr_part` と同じフォールバック方針）。

### その他の挙動

- topology 列が `circular` の場合は既存どおり `TOPOLOGY /circular` 行が付く（変更なし）。
- segment は「配置済み」type として扱う。既存の `_is_unplaced()` は
  `type == "unplaced"` のみを unplaced とみなすため、segment があると
  `is_wgs = False` になる（変更不要）。

## コード変更

### `src/common/source_builder.py`（共有）

1. `SequenceRoleEntry` の `type` コメントおよび `load_sequence_roles()` docstring に
   `segment` を追記（type にホワイトリスト検証は無いため検証コードの追加は不要）。

2. `source_qualifier(entry, seq_id, is_wgs=False, segment_count=0)` に
   `segment_count` 引数を追加。`entry.type == "segment"` 分岐を追加し:
   - `segment_count <= 1` → `{}`（`/segment` を出さない）
   - `segment_count >= 2` → `{"segment": entry.seq_name}`（`seq_name` が空なら `{}`）

3. `ff_definition(..., chromosome_count=0, segment_count=0)` に `segment_count` を追加。
   `entry.type == "segment"` 分岐を追加し上表の 4 パターンを実装。

いずれの新引数もデフォルト値を持たせ、既存呼び出し（テスト含む）を壊さない。

### `src/mss_builder/ann_writer.py` / `src/egapx2mss/ann_writer.py`（両方）

既存の `chromosome_count` 計算の直後に同型の集計を追加:

```python
segment_count = 0
if sequence_roles:
    segment_count = sum(1 for e in sequence_roles.values() if e.type == "segment")
```

`source_qualifier(...)` と `ff_definition(...)` の呼び出しに `segment_count=segment_count` を渡す。

## テスト

`tests/test_ff_definition_molecule.py` に追加（既存の chromosome count テストのパターン踏襲）:

- `ff_definition` segment 単一 complete → `... {mol}, complete genome`
- `ff_definition` segment 単一 partial → `... {mol}, partial genome`
- `ff_definition` segment 複数 complete → `... {mol}, segment {name}, complete sequence`
- `ff_definition` segment 複数 partial → `... {mol}, segment {name}, partial sequence`
- `source_qualifier` segment 単一 → `segment` キーを含まない
- `source_qualifier` segment 複数 → `{"segment": name}` を返す

可能であれば mss_builder のエンドツーエンド（segment 行を含む `sequence_roles.tsv` → `.ann`）で
`/segment` 行と DEFINITION が期待どおり出ることを確認する統合テストを 1 本追加する。

## ドキュメント更新

- `CLAUDE.md`: sequence role の `type` 一覧に `segment` を追加し、ff_definition 表に segment 行を追記。
- `README.md`: 同様に sequence role の説明へ `segment` を追記。
- `docs/design/submission_category.md` は本件の対象外（カテゴリ定義であり sequence role とは別）。

## スコープ外（YAGNI）

- エントリごとに異なる mol_type / segment 名の高度なパターン（現行の per-entry 仕組みで足りる）。
- `type` のホワイトリスト検証追加（既存でも未検証であり本件の目的ではない）。
- egapx2mss の入力 (.tbl/.asn) 側で segment を自動推定する仕組み（sequence role で明示指定する前提）。
