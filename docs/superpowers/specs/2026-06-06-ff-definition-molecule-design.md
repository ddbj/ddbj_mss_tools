# ff_definition の分子トークンを mol_type に応じて可変化 設計

- 日付: 2026-06-06
- 対象: `src/common/source_builder.py`（呼び出し側: `src/egapx2mss/ann_writer.py`, `src/mss_builder/ann_writer.py`）
- 影響ツール: egapx2mss / mss_builder / batch_wgs_builder

## 背景・課題

`ff_definition`（DDBJ Flat File の DEFINITION 行になる文字列）は現在すべて `"<prefix> DNA, <記述子>"` の形で **"DNA" がハードコード**されている。RNA ウイルスや RNA 由来の登録では分子種が DNA でないため、source feature の `mol_type` qualifier に応じて "DNA" 部分を変えたい。

"DNA" のハードコードは3関数・計14テンプレートに散在:
- 経路A `ff_definition()`（`source_builder.py:78-111`、sequence role ベース実値、7パターン）
- 経路B-1 `create_source_feature()`（`:186-233`、category+seq_type メタ記法、5パターン）
- 経路B-2 `_create_source_with_meta()`（`:236-`、SOURCE_IDENTIFIER メタ記法、2パターン）

## 要件

1. `ff_definition` の `"<prefix> DNA"` の "DNA" 部分を、source の `mol_type` qualifier の値から決定する。
2. 全3経路に適用する（一貫性）。
3. mol_type が既定の `genomic DNA` や DNA 系のときは従来通り `DNA`（既存出力は不変＝後方互換）。
4. mol_type が空・未定義・どのパターンにも該当しない場合はデフォルト `DNA`。

## 分子トークン判定（新規ヘルパー）

`src/common/source_builder.py` に純粋関数を追加:

```python
def _molecule_token(mol_type: str | None) -> str:
    """Decide the molecule token used in ff_definition ("<prefix> <token> DNA,").

    Rules (in order):
      1. empty / None              -> "DNA"
      2. contains tRNA/rRNA/mRNA   -> that token (case-SENSITIVE; INSDC fixed spelling)
      3. lowercased contains "dna" -> "DNA"
      4. lowercased contains "rna" -> "RNA"
      5. otherwise (e.g. "protein")-> "DNA" (default)
    """
    if not mol_type:
        return "DNA"
    for token in ("tRNA", "rRNA", "mRNA"):
        if token in mol_type:
            return token
    low = mol_type.lower()
    if "dna" in low:
        return "DNA"
    if "rna" in low:
        return "RNA"
    return "DNA"
```

### 判定表（確定）

| mol_type | 出力トークン |
|---|---|
| `genomic DNA`（既定） / `other DNA` / `unassigned DNA` | `DNA` |
| `mRNA` | `mRNA` |
| `tRNA` | `tRNA` |
| `rRNA` | `rRNA` |
| `genomic RNA` / `transcribed RNA` / `viral cRNA` | `RNA` |
| `protein` / その他非該当 | `DNA`（デフォルト） |
| 空 / None / 未定義 | `DNA`（デフォルト） |

判定の設計上の注記:
- tRNA/rRNA/mRNA を **DNA/RNA より先**に case-sensitive で判定する。これにより `mRNA` の中の `RNA` に引っ張られず `mRNA` を返す。INSDC では DNA系の値（`genomic DNA` 等）が tRNA/rRNA/mRNA 部分文字列を含むことは無いため、実務上「DNA 最優先」と同じ結果になる。
- `viral cRNA` は `tRNA/rRNA/mRNA` を含まないため、小文字化 `dna` 判定（含まない）→ `rna` 判定（含む）で `RNA` になる。

## 適用箇所

各テンプレートの `DNA`（文末記述子の `complete genome` 等ではなく、`<prefix> DNA` の分子部分のみ）を `_molecule_token(...)` の戻り値に置換する。

### 経路A: `ff_definition()` — シグネチャ変更あり

現在:
```python
def ff_definition(entry, seq_id, organism, infraspecific_name_modifier, is_wgs: bool = False) -> str:
```

変更後（`is_wgs` の **前** に必須引数 `mol_type` を挿入。デフォルトなし）:
```python
def ff_definition(entry, seq_id, organism, infraspecific_name_modifier, mol_type, is_wgs: bool = False) -> str:
```

本体: 先頭で `mol = _molecule_token(mol_type)` を作り、全 return 文の `... DNA, ...` を `... {mol}, ...` に置換。7パターン:

| type/status | 変更後テンプレート |
|---|---|
| unplaced + is_wgs | `{prefix} {mol}, {seq_id}` |
| unplaced + 非wgs | `{prefix} {mol}, unplaced sequence {seq_id}` |
| chromosome + complete | `{prefix} {mol}, {chr_part}, complete sequence` |
| chromosome + その他 | `{prefix} {mol}, {chr_part}, unlocalized sequence {seq_id}` |
| organelle + complete | `{prefix} {mol}, {organelle_name}, complete sequence` |
| organelle + その他 | `{prefix} {mol}, {organelle_name}, partial sequence` |
| fallback | `{prefix} {mol}, {seq_id}` |

### 経路B-1: `create_source_feature()`

`mol_type = source_dict.get("mol_type", "genomic DNA")` が既にある（`:194`）。直後に `mol = _molecule_token(mol_type)` を追加し、5テンプレートの `DNA` を `{mol}` に:

| 条件 | 変更後 |
|---|---|
| datatype=WGS | `@@[organism]@@ @@[{modifier}]@@ {mol}, @@[submitter_seqid]@@` |
| complete | `@@[organism]@@ @@[{modifier}]@@ {mol}, complete genome` |
| nearly complete | `@@[organism]@@ @@[{modifier}]@@ {mol}, nearly complete genome` |
| plasmid | `@@[organism]@@ @@[{modifier}]@@ plasmid @@[plasmid]@@ {mol}, complete sequence` |
| その他 | `@@[organism]@@ @@[{modifier}]@@ {mol}, @@[submitter_seqid]@@` |

### 経路B-2: `_create_source_with_meta()`

`mol_type = source_dict.get("mol_type", "genomic DNA")` が既にある（`:252`）。`mol = _molecule_token(mol_type)` を追加し、2テンプレートの `DNA` を `{mol}` に:

| 条件 | 変更後 |
|---|---|
| modifier キーあり | `@@[organism]@@ @@[{source_modifier_key}]@@ {mol}, @@[submitter_seqid]@@` |
| なし | `@@[organism]@@ {mol}, @@[submitter_seqid]@@` |

## 呼び出し側の修正（2箇所）

経路A の `ff_definition()` 呼び出しに `mol_type` を渡す。両箇所とも `base_source`（= common.SOURCE のコピー）から取得できる。

### `src/egapx2mss/ann_writer.py:181`
```python
# 変更前
source_quals["ff_definition"] = ff_definition(
    role_entry, entry_id, organism, infraspecific_name_modifier, is_wgs
)
# 変更後（mol_type を is_wgs の前に）
source_quals["ff_definition"] = ff_definition(
    role_entry, entry_id, organism, infraspecific_name_modifier,
    base_source.get("mol_type", ""), is_wgs
)
```
※ この関数内の `ff_definition` 呼び出しは実際には2箇所ある可能性に注意（is_wgs=True/False 経路）。grep で全呼び出しを確認し、すべてに `mol_type` 引数を追加する。

### `src/mss_builder/ann_writer.py:177`
```python
# 変更前
source_quals["ff_definition"] = ff_definition(
    role_entry, entry_id, organism, infraspecific_name_modifier, is_wgs=False
)
# 変更後
source_quals["ff_definition"] = ff_definition(
    role_entry, entry_id, organism, infraspecific_name_modifier,
    base_source.get("mol_type", ""), is_wgs=False
)
```

mol_type が無い場合は空文字 `""` を渡し、`_molecule_token` が `DNA` を返す（後方互換）。

## 後方互換

- 既定 `genomic DNA` や DNA 系 mol_type → 従来通り `DNA`。既存の全 example 出力は不変。
- mol_type 未設定（`base_source` に mol_type が無い）→ `""` → `DNA`。
- 変化するのは RNA 系 mol_type を明示指定した場合のみ。

## スコープ外

- `plasmid` 行の分子部分も置換対象に含める（一貫性。経路B-1の1パターン）。
- 記述子部分（`complete genome` / `complete sequence` / `unplaced sequence` 等）は変更しない。
- mol_type の値そのものの検証（不正値の拒否）はしない。非該当値は黙って `DNA`。

## テスト方針

`tests/test_ff_definition_molecule.py`（新規、pytest）:
- **`_molecule_token` 単体**: `genomic DNA`→DNA、`mRNA`→mRNA、`tRNA`→tRNA、`rRNA`→rRNA、`genomic RNA`→RNA、`transcribed RNA`→RNA、`viral cRNA`→RNA、`protein`→DNA、`""`→DNA、`None`→DNA。
- **経路A `ff_definition()`**: 同じ entry に対し mol_type を変え、出力の分子トークンが期待通り変わること。既定で `DNA` になる後方互換も確認（例: unplaced+wgs で `genomic DNA` → `... DNA, seq1`、`genomic RNA` → `... RNA, seq1`）。
- **経路B-1 `create_source_feature()`**: source_dict に `mol_type: "genomic RNA"` を入れ、ff_definition 行に `RNA` が入ることを確認。`mol_type` 省略時（既定 genomic DNA）は `DNA`。
- **経路B-2 `create_source_feature(..., use_meta_expression=True)`**: 同様に meta 経路で `RNA` 反映を確認。
- 実行: `PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_ff_definition_molecule.py -q`
