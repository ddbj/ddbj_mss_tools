# ff_definition のメタ記法化 設計

## 背景と目的

`src/common/source_builder.py` の `ff_definition()`（sequence role による per-entry パス）は現在、
`{organism} {識別子の値}` や `{seq_name}` を**具体値に展開した文字列**を生成している。

これを DDBJ MSS のメタ記法（`@@[qualifier_name]@@`）を使ったテンプレート文字列に変更する。
source フィーチャーが持つ qualifier（`/organism`, `/cultivar`, `/chromosome`, `/plasmid`, `/segment`,
`/submitter_seqid`）や MSS 提供の `@@[entry]@@` を参照させ、MSS 側で登録時に展開させる。

対象ツール: **mss_builder / egapx2mss / gff2mss** の3つ（共有 `ff_definition()` を使う全経路）。

WGS / COMMON ブロック用の `create_source_feature()`（`_create_source_with_meta`）は既にメタ記法を使っており、
本設計はそれと整合する形で per-entry パスもメタ記法にする。

## シグネチャ変更

現状:
```python
def ff_definition(entry, seq_id, organism, infraspecific_name_modifier, mol_type,
                  is_wgs=False, chromosome_count=0, segment_count=0) -> str:
```

変更後:
```python
def ff_definition(entry, source_identifier, mol_type,
                  is_wgs=False, chromosome_count=0, segment_count=0) -> str:
```

- **削除**: `seq_id`（新ルールで未使用）、`organism`（値。→ リテラル `@@[organism]@@`）、
  `infraspecific_name_modifier`（値。→ 名前を使う `source_identifier` に置換）
- **追加**: `source_identifier`（`SOURCE_IDENTIFIER` の qualifier **名**。例 `"cultivar"`。空/None 可）
- **継続**: `entry`（`Optional[SequenceRoleEntry]`）, `mol_type`, `is_wgs`, `chromosome_count`, `segment_count`

## prefix と mol

- **`{P}`（prefix）**:
  - `source_identifier` が非空 → `@@[organism]@@ @@[{source_identifier}]@@`
  - `source_identifier` が空/None → `@@[organism]@@`
- **`{mol}`**: 従来どおり `_molecule_token(mol_type)` の具体値（DNA/RNA/tRNA/rRNA/mRNA）

## 決定テーブル（上から順に判定、最初に一致した分岐）

| # | type | 条件 | ff_definition |
|---|------|------|---------------|
| 1 | unplaced / entry=None | is_wgs=true | `{P} {mol}, @@[submitter_seqid]@@` |
| 2 | unplaced / entry=None | is_wgs=false | `{P} {mol}, unplaced sequence @@[entry]@@` |
| 3 | chromosome | count==1・complete | `{P} {mol}, chromosome, complete genome` |
| 4 | chromosome | count≥2・complete | `{P} {mol}, chromosome @@[chromosome]@@, complete sequence`（seq_name空→`chromosome`）|
| 5a | chromosome | count==1・partial | `{P} {mol}, chromosome` |
| 5b | chromosome | count≥2・partial | `{P} {mol}, chromosome @@[chromosome]@@`（seq_name空→`chromosome`）|
| 6 | organelle | complete | `{P} {organelle_code} {mol}, complete genome` |
| 7 | organelle | partial | `{P} {organelle_code} {mol}, partial genome` |
| 8 | plasmid | complete | `{P} plasmid @@[plasmid]@@ {mol}, complete sequence`（seq_name空→ValueError）|
| 9 | plasmid | partial | `{P} plasmid @@[plasmid]@@ {mol}, partial sequence`（seq_name空→ValueError）|
| 10 | segment | count==1・complete | `{P} {mol}, complete genome` |
| 11 | segment | count==1・partial | `{P} {mol}, partial genome` |
| 12 | segment | count≥2・complete | `{P} {mol}, segment @@[segment]@@, complete sequence`（seq_name空→`segment, complete sequence`）|
| 13 | segment | count≥2・partial | `{P} {mol}, segment @@[segment]@@`（seq_name空→`segment`）|
| 14 | その他（未知 type） | — | `{P} {mol}, @@[entry]@@` |

### ルールの要点

- **chromosome は count で分岐**（従来 #3/#4 は complete のみ count 分岐だったが、partial も分岐する）:
  - 単一（count==1）は complete/partial とも `@@[chromosome]@@` を出さず `chromosome`（番号なし）
  - 複数（count≥2）は `chromosome @@[chromosome]@@`
- **organelle** は prefix のみメタ化。`{organelle_code}`（`mitochondrial` 等の形容詞形）は
  従来どおり `_organelle_code(entry.seq_name)` の変換値のまま（`@@[organelle]@@` だと生値
  `mitochondrion` になり不整合のため）。
- **plasmid の seq_name 空は `ValueError`**（`ff_definition()` 内で送出）。メッセージ例:
  `"plasmid entry requires a non-empty seq_name"`。
- **seq_name 空時のベア語フォールバック**（chromosome #4/#5b, segment #12/#13）は、
  展開先 qualifier が無いのに `@@[...]@@` を残さないための措置。

## メタ参照と source qualifier の対応（整合性）

各メタ参照は、同じ source フィーチャーが持つ qualifier（または MSS 提供値）に解決される。
`source_qualifier()` は本設計では変更しない。

| メタ参照 | 解決先 | 供給元 |
|---|---|---|
| `@@[organism]@@` | `/organism` | base_source（common.SOURCE） |
| `@@[{source_identifier}]@@` | `/cultivar` 等 | base_source（SOURCE_IDENTIFIER の qualifier） |
| `@@[chromosome]@@` | `/chromosome` | `source_qualifier`（seq_name 非空時） |
| `@@[plasmid]@@` | `/plasmid` | `source_qualifier`（plasmid は空不可） |
| `@@[segment]@@` | `/segment` | `source_qualifier`（count≥2 かつ seq_name 非空時） |
| `@@[submitter_seqid]@@` | `/submitter_seqid` | `source_qualifier`（unplaced+is_wgs 時） |
| `@@[entry]@@` | エントリ名 | MSS 提供 |

注意（本設計のスコープ外だが記録）: 単一 chromosome（count==1, seq_name 非空）では ff_definition は
`@@[chromosome]@@` を参照しないが、`source_qualifier` は `/chromosome` を出力し続ける（余分な
qualifier だが害はない）。`source_qualifier` の変更は本タスクでは行わない。

## 呼び出し側の変更

`ff_definition(...)` を呼ぶ3箇所を新シグネチャに合わせる:
- `src/mss_builder/ann_writer.py`（`if not is_wgs:` 内。source_id_key を渡す）
- `src/egapx2mss/ann_writer.py`（is_wgs を渡す。source_id_key を渡す）
- `src/gff2mss/assemble.py`（src_id_key を渡す）

各 writer で ff_definition 専用に算出していた `organism`（値）/`infraspecific_name_modifier`（値）は、
他で使われていなければ削除する（要 grep 確認）。`source_id_key`（名前）は既に算出済み。

## テスト

- `tests/test_ff_definition_molecule.py`: `ff_definition` の全アサーションを新シグネチャ・メタ記法へ
  書き換え。テーブルの各行（#1〜#14、空 seq_name のフォールバック、plasmid 空名の `ValueError`）を網羅。
  `_molecule_token` / `_organelle_code` の単体テストは無変更。`create_source_feature`（メタ経路）の
  テストも無変更。
- `tests/test_mss_segment.py`: 統合アサーションをメタ記法へ更新
  （例: `Influenza A virus RNA, segment 4, complete sequence` →
  `@@[organism]@@ {mol}, segment @@[segment]@@, complete sequence` 相当。source_identifier 無しの
  common なら prefix は `@@[organism]@@` のみ）。`/segment` qualifier 行の検証は維持。
- 3 writer の呼び出し変更後、既存の他テストが壊れないことを全体回帰で確認。

## ドキュメント更新

- `CLAUDE.md`: 「sequence role と ff_definition」節の表をメタ記法版に全面差し替え、prefix/mol の説明追記。
- `README.md`: sequence role 節の ff_definition 記述をメタ記法へ更新（日・英）。

## スコープ外（YAGNI）

- `source_qualifier()` の挙動変更（単一 chromosome の `/chromosome` 抑制など）。
- WGS/COMMON メタ経路 `create_source_feature()` の変更（既にメタ記法）。
- MSS が per-entry source フィーチャー内の `@@[...]@@` を展開する仕様の検証（DDBJ 側仕様に依拠。
  本設計はテンプレート生成のみ担当）。
