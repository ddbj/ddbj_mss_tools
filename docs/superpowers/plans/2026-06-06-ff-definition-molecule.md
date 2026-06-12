# ff_definition の分子トークンを mol_type に応じて可変化 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `ff_definition` の `"<prefix> DNA, ..."` の "DNA" 部分を、source の `mol_type` qualifier に応じて DNA / RNA / tRNA / rRNA / mRNA に切り替える（全3生成経路）。

**Architecture:** `src/common/source_builder.py` に純粋関数 `_molecule_token(mol_type)` を追加し、3つの ff_definition 生成関数（`ff_definition` 実値経路、`create_source_feature` メタ経路、`_create_source_with_meta` メタ経路）の各テンプレートの "DNA" を `_molecule_token(...)` の戻り値に置換。`ff_definition` には `mol_type` 必須引数を `is_wgs` の前に追加し、呼び出し側2箇所（egapx2mss / mss_builder）で `base_source.get("mol_type", "")` を渡す。

**Tech Stack:** Python 3.10+, pytest

---

## Context & Constraints

- 設計ドキュメント: `docs/superpowers/specs/2026-06-06-ff-definition-molecule-design.md`
- 影響ツール: egapx2mss / mss_builder / batch_wgs_builder（後2者は共有の `create_source_feature` 経由）
- **テスト実行**: `PY=/Users/tanizawa/miniforge3/envs/dr_tools/bin/python3`（pytest 導入済み）。`PYTHONPATH=src $PY -m pytest <args>`
- **判定ルール（確定）**: 空/None→`DNA`、`tRNA`/`rRNA`/`mRNA` を含む（case-sensitive）→その値、小文字化して `dna` を含む→`DNA`、`rna` を含む→`RNA`、それ以外→`DNA`。
- `ff_definition()` の呼び出しは正確に2箇所のみ（grep 済み）: `src/egapx2mss/ann_writer.py:181`, `src/mss_builder/ann_writer.py:177`。

### 現行コードの該当箇所（変更前）

`ff_definition()` 本体（`source_builder.py:88-111`）— prefix 後の各 return の "DNA":
```python
    prefix = f"{organism} {infraspecific_name_modifier}".strip() if infraspecific_name_modifier else organism
    if entry is None or entry.type == "unplaced":
        if is_wgs:
            return f"{prefix} DNA, {seq_id}"
        else:
            return f"{prefix} DNA, unplaced sequence {seq_id}"
    if entry.type == "chromosome":
        chr_part = f"chromosome {entry.seq_name}".strip() if entry.seq_name else "chromosome"
        if entry.status == "complete":
            return f"{prefix} DNA, {chr_part}, complete sequence"
        else:
            return f"{prefix} DNA, {chr_part}, unlocalized sequence {seq_id}"
    if entry.type == "organelle":
        organelle_name = entry.seq_name
        if entry.status == "complete":
            return f"{prefix} DNA, {organelle_name}, complete sequence"
        else:
            return f"{prefix} DNA, {organelle_name}, partial sequence"
    # fallback
    return f"{prefix} DNA, {seq_id}"
```

`create_source_feature()`（`:194-211`）— `mol_type` 取得済み、5テンプレートの "DNA":
```python
    mol_type = source_dict.get("mol_type", "genomic DNA")
    if rules.datatype == "WGS":
        submitter_seqid = "@@[entry]@@"
        ff_def = f"@@[organism]@@ @@[{modifier}]@@ DNA, @@[submitter_seqid]@@"
    else:
        if seq_type in ["c", "complete"]:
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ DNA, complete genome"
        elif seq_type in ["n", "nearly complete", "nearly-complete"]:
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ DNA, nearly complete genome"
        elif seq_type in ["p", "plasmid"]:
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ plasmid @@[plasmid]@@ DNA, complete sequence"
            plasmid = True
        else:
            submitter_seqid = "@@[entry]@@"
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ DNA, @@[submitter_seqid]@@"
```

`_create_source_with_meta()`（`:252-257`）— `mol_type` 取得済み、2テンプレート:
```python
    mol_type = source_dict.get("mol_type", "genomic DNA")
    if source_modifier_key:
        ff_def = f"@@[organism]@@ @@[{source_modifier_key}]@@ DNA, @@[submitter_seqid]@@"
    else:
        ff_def = "@@[organism]@@ DNA, @@[submitter_seqid]@@"
```

### File Structure

| ファイル | 役割 | 操作 |
|---|---|---|
| `src/common/source_builder.py` | `_molecule_token` 追加、3関数の DNA 置換、`ff_definition` に mol_type 引数 | 変更 |
| `src/egapx2mss/ann_writer.py` | `ff_definition()` 呼び出しに mol_type 引数追加 | 変更 |
| `src/mss_builder/ann_writer.py` | 同上 | 変更 |
| `tests/test_ff_definition_molecule.py` | `_molecule_token` 単体 + 3経路結合テスト | 新規 |

---

## Task 1: `_molecule_token` ヘルパー（TDD）

**Files:**
- Modify: `src/common/source_builder.py`
- Test: `tests/test_ff_definition_molecule.py`

- [ ] **Step 1: 失敗するテストを作成**

`tests/test_ff_definition_molecule.py` を新規作成:

```python
"""Tests for mol_type-driven ff_definition molecule token."""

import pytest

from common.source_builder import _molecule_token


@pytest.mark.parametrize("mol_type,expected", [
    ("genomic DNA", "DNA"),
    ("other DNA", "DNA"),
    ("unassigned DNA", "DNA"),
    ("mRNA", "mRNA"),
    ("tRNA", "tRNA"),
    ("rRNA", "rRNA"),
    ("genomic RNA", "RNA"),
    ("transcribed RNA", "RNA"),
    ("viral cRNA", "RNA"),
    ("protein", "DNA"),
    ("", "DNA"),
    (None, "DNA"),
])
def test_molecule_token(mol_type, expected):
    assert _molecule_token(mol_type) == expected
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools && PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_ff_definition_molecule.py -q`
Expected: FAIL（`ImportError: cannot import name '_molecule_token'`）。

- [ ] **Step 3: `_molecule_token` を実装**

`src/common/source_builder.py` の `def ff_definition(` 定義の **直前**（現78行付近）に追加:

```python
def _molecule_token(mol_type: str | None) -> str:
    """Decide the molecule token used in ff_definition ("<prefix> <token>, ...").

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

- [ ] **Step 4: テストが通ることを確認**

Run: 同上
Expected: 12 passed。

- [ ] **Step 5: コミット**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/common/source_builder.py tests/test_ff_definition_molecule.py
git commit -m "feat: add _molecule_token helper for ff_definition

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: 経路A `ff_definition()` に mol_type を反映（TDD）

**Files:**
- Modify: `src/common/source_builder.py`（`ff_definition` 本体・シグネチャ）
- Test: `tests/test_ff_definition_molecule.py`（追記）

- [ ] **Step 1: 失敗するテストを追記**

`tests/test_ff_definition_molecule.py` の末尾に追加。`SequenceRoleEntry` を使い、mol_type 違いで分子トークンが変わること・既定で DNA になることを検証:

```python
from common.source_builder import ff_definition, SequenceRoleEntry


def test_ff_definition_unplaced_wgs_dna():
    # entry=None -> unplaced; is_wgs=True
    out = ff_definition(None, "seq1", "Homo sapiens", "", "genomic DNA", is_wgs=True)
    assert out == "Homo sapiens DNA, seq1"


def test_ff_definition_unplaced_wgs_rna():
    out = ff_definition(None, "seq1", "Homo sapiens", "", "genomic RNA", is_wgs=True)
    assert out == "Homo sapiens RNA, seq1"


def test_ff_definition_unplaced_wgs_mrna():
    out = ff_definition(None, "seq1", "Homo sapiens", "", "mRNA", is_wgs=True)
    assert out == "Homo sapiens mRNA, seq1"


def test_ff_definition_chromosome_complete_rna():
    e = SequenceRoleEntry("seq1", "chromosome", "1", "complete", False)
    out = ff_definition(e, "seq1", "Homo sapiens", "strainX", "genomic RNA", is_wgs=False)
    assert out == "Homo sapiens strainX RNA, chromosome 1, complete sequence"


def test_ff_definition_empty_mol_type_defaults_dna():
    out = ff_definition(None, "seq1", "Homo sapiens", "", "", is_wgs=True)
    assert out == "Homo sapiens DNA, seq1"
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools && PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_ff_definition_molecule.py -q`
Expected: 新規テストが FAIL（`ff_definition` が mol_type 引数を受け取らず TypeError、または "DNA" 固定で RNA 期待が一致しない）。

- [ ] **Step 3: `ff_definition` を変更**

シグネチャ（現78-79行）を変更し、`mol_type` を `is_wgs` の前に必須挿入:

```python
def ff_definition(entry: Optional[SequenceRoleEntry], seq_id: str, organism: str,
                  infraspecific_name_modifier: str, mol_type: str, is_wgs: bool = False) -> str:
```

本体（現88-111行）を以下に置換（prefix 直後に `mol` を作り、全 return の `DNA` を `{mol}` に）:

```python
    prefix = f"{organism} {infraspecific_name_modifier}".strip() if infraspecific_name_modifier else organism
    mol = _molecule_token(mol_type)

    if entry is None or entry.type == "unplaced":
        if is_wgs:
            return f"{prefix} {mol}, {seq_id}"
        else:
            return f"{prefix} {mol}, unplaced sequence {seq_id}"

    if entry.type == "chromosome":
        chr_part = f"chromosome {entry.seq_name}".strip() if entry.seq_name else "chromosome"
        if entry.status == "complete":
            return f"{prefix} {mol}, {chr_part}, complete sequence"
        else:
            return f"{prefix} {mol}, {chr_part}, unlocalized sequence {seq_id}"

    if entry.type == "organelle":
        organelle_name = entry.seq_name
        if entry.status == "complete":
            return f"{prefix} {mol}, {organelle_name}, complete sequence"
        else:
            return f"{prefix} {mol}, {organelle_name}, partial sequence"

    # fallback
    return f"{prefix} {mol}, {seq_id}"
```

- [ ] **Step 4: テストが通ることを確認**

Run: 同上
Expected: 17 passed（12 + 5）。

- [ ] **Step 5: コミット**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/common/source_builder.py tests/test_ff_definition_molecule.py
git commit -m "feat: ff_definition uses mol_type for molecule token (path A)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: 呼び出し側2箇所に mol_type を渡す

**Files:**
- Modify: `src/egapx2mss/ann_writer.py:177-179`
- Modify: `src/mss_builder/ann_writer.py:177-179`

この時点で `ff_definition` は mol_type 必須になったため、呼び出し側を直さないと egapx2mss / mss_builder が TypeError になる。両方修正する。

- [ ] **Step 1: egapx2mss の呼び出しを修正**

`src/egapx2mss/ann_writer.py` の現在のコード:
```python
        source_quals["ff_definition"] = ff_definition(
            role_entry, entry_id, organism, infraspecific_name_modifier, is_wgs
        )
```
を以下に変更（`mol_type` を `is_wgs` の前に）:
```python
        source_quals["ff_definition"] = ff_definition(
            role_entry, entry_id, organism, infraspecific_name_modifier,
            base_source.get("mol_type", ""), is_wgs
        )
```

- [ ] **Step 2: mss_builder の呼び出しを修正**

`src/mss_builder/ann_writer.py` の現在のコード:
```python
            source_quals["ff_definition"] = ff_definition(
                role_entry, entry_id, organism, infraspecific_name_modifier, is_wgs=False
            )
```
を以下に変更:
```python
            source_quals["ff_definition"] = ff_definition(
                role_entry, entry_id, organism, infraspecific_name_modifier,
                base_source.get("mol_type", ""), is_wgs=False
            )
```

- [ ] **Step 3: 両ファイルで `base_source` が存在することを確認**

`base_source` は両ファイルで `ff_definition` 呼び出しより前に定義済み（egapx2mss は `:131`、mss_builder も同様に common.SOURCE から構築）。grep で確認:

Run: `cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools && grep -n "base_source" src/egapx2mss/ann_writer.py src/mss_builder/ann_writer.py`
Expected: 両ファイルで `base_source` の定義行が呼び出し行より前にある。もし mss_builder で変数名が異なる場合は、その箇所の source 辞書（common.SOURCE 由来の dict）から `.get("mol_type", "")` を取るよう合わせる。

- [ ] **Step 4: import が壊れていないか・両ツールが起動するか確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -c "import egapx2mss.ann_writer, mss_builder.ann_writer; print('import OK')"
```
Expected: `import OK`

- [ ] **Step 5: 既存テスト全体が通ることを確認（回帰チェック）**

Run: `cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools && PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest -q`
Expected: 全 pass（17 + 既存 18 = 35 passed）。

- [ ] **Step 6: コミット**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/egapx2mss/ann_writer.py src/mss_builder/ann_writer.py
git commit -m "feat: pass mol_type to ff_definition from egapx2mss and mss_builder

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: 経路B-1/B-2（メタ記法）に mol_type を反映（TDD）

**Files:**
- Modify: `src/common/source_builder.py`（`create_source_feature`, `_create_source_with_meta`）
- Test: `tests/test_ff_definition_molecule.py`（追記）

- [ ] **Step 1: 失敗するテストを追記**

`tests/test_ff_definition_molecule.py` の末尾に追加。メタ記法経路で ff_definition 行の分子トークンが mol_type に応じて変わることを検証。`create_source_feature` の戻り行から ff_definition 値を取り出す:

```python
from common.source_builder import create_source_feature


def _ff_def_value(rows):
    for r in rows:
        if r[3] == "ff_definition":
            return r[4]
    return None


def test_create_source_meta_b2_rna():
    # use_meta_expression=True -> _create_source_with_meta (path B-2)
    src = {"organism": "X", "mol_type": "genomic RNA"}
    rows = create_source_feature("WGS", None, None, None, src,
                                 source_modifier_key="strain",
                                 use_meta_expression=True)
    assert _ff_def_value(rows) == "@@[organism]@@ @@[strain]@@ RNA, @@[submitter_seqid]@@"


def test_create_source_meta_b2_default_dna():
    src = {"organism": "X"}  # no mol_type -> default genomic DNA -> DNA
    rows = create_source_feature("WGS", None, None, None, src,
                                 source_modifier_key="strain",
                                 use_meta_expression=True)
    assert _ff_def_value(rows) == "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"


def test_create_source_b1_complete_rna():
    # path B-1: non-WGS category (GNM), seq_type=complete
    src = {"organism": "X", "mol_type": "genomic RNA"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    assert _ff_def_value(rows) == "@@[organism]@@ @@[strain]@@ RNA, complete genome"


def test_create_source_b1_default_dna():
    src = {"organism": "X"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    assert _ff_def_value(rows) == "@@[organism]@@ @@[strain]@@ DNA, complete genome"
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools && PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_ff_definition_molecule.py -q`
Expected: 4新規テストのうち RNA 系（`_rna`）が FAIL（現状 "DNA" 固定）。`_default_dna` は元々 pass する。

- [ ] **Step 3: `create_source_feature` を変更（経路B-1）**

`src/common/source_builder.py` の `mol_type = source_dict.get("mol_type", "genomic DNA")` の直後に `mol = _molecule_token(mol_type)` を追加し、5テンプレートの `DNA` を `{mol}` に置換:

```python
    mol_type = source_dict.get("mol_type", "genomic DNA")
    mol = _molecule_token(mol_type)

    # WGS-family: source goes in COMMON block, submitter_seqid always set
    if rules.datatype == "WGS":
        submitter_seqid = "@@[entry]@@"
        ff_def = f"@@[organism]@@ @@[{modifier}]@@ {mol}, @@[submitter_seqid]@@"
    else:
        # Per-entry source (GNM, MAG, etc.)
        if seq_type in ["c", "complete"]:
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ {mol}, complete genome"
        elif seq_type in ["n", "nearly complete", "nearly-complete"]:
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ {mol}, nearly complete genome"
        elif seq_type in ["p", "plasmid"]:
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ plasmid @@[plasmid]@@ {mol}, complete sequence"
            plasmid = True
        else:
            submitter_seqid = "@@[entry]@@"
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ {mol}, @@[submitter_seqid]@@"
```

(注: `mol_type` 変数はこの後 `ret.append(["", "source", "1..E", "mol_type", mol_type])` で source 行にそのまま使われているので、`mol_type` 変数自体は残すこと。`mol` は新規追加。)

- [ ] **Step 4: `_create_source_with_meta` を変更（経路B-2）**

`mol_type = source_dict.get("mol_type", "genomic DNA")` の直後に `mol = _molecule_token(mol_type)` を追加し、2テンプレートの `DNA` を `{mol}` に:

```python
    mol_type = source_dict.get("mol_type", "genomic DNA")
    mol = _molecule_token(mol_type)

    if source_modifier_key:
        ff_def = f"@@[organism]@@ @@[{source_modifier_key}]@@ {mol}, @@[submitter_seqid]@@"
    else:
        ff_def = f"@@[organism]@@ {mol}, @@[submitter_seqid]@@"
```

- [ ] **Step 5: テストが通ることを確認**

Run: `cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools && PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_ff_definition_molecule.py -q`
Expected: 21 passed（17 + 4）。

- [ ] **Step 6: 全テスト回帰確認**

Run: `cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools && PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest -q`
Expected: 全 pass（21 + 既存 18 = 39 passed）。

- [ ] **Step 7: コミット**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/common/source_builder.py tests/test_ff_definition_molecule.py
git commit -m "feat: meta-notation ff_definition uses mol_type (paths B-1, B-2)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: 既存 example での後方互換の目視確認

**Files:** （検証のみ）

- [ ] **Step 1: mss_builder を既存 example で実行し ff_definition が DNA のままか確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PY=/Users/tanizawa/miniforge3/envs/dr_tools/bin/python3
PYTHONPATH=src $PY -m mss_builder.cli examples/mss_builder/GruXX_test.fa \
  -o /tmp/molcheck -p out --common examples/mss_builder/example.gnm.common.json \
  --sequence_roles examples/mss_builder/example.sequence_roles.tsv > /tmp/molcheck.log 2>&1
echo "RC=$?"
grep -n "ff_definition" /tmp/molcheck/out.ann
```
Expected: 実行成功（RC=0）。example の common JSON は mol_type 未指定 or `genomic DNA` のため、ff_definition 行は従来通り `... DNA, ...`（RNA 等にならない）。もし example ファイルが存在せずエラーなら、`examples/mss_builder/` 内の他の .fa + common JSON + sequence_roles の組で代替（`ls examples/mss_builder/`）。

- [ ] **Step 2: RNA を指定したケースの動作確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PY=/Users/tanizawa/miniforge3/envs/dr_tools/bin/python3
PYTHONPATH=src $PY -c "
import json, copy
base = json.load(open('examples/mss_builder/example.gnm.common.json'))
base.setdefault('SOURCE', {})['mol_type'] = 'genomic RNA'
json.dump(base, open('/tmp/rna_common.json','w'))
print('wrote /tmp/rna_common.json with mol_type=genomic RNA')
"
PYTHONPATH=src $PY -m mss_builder.cli examples/mss_builder/GruXX_test.fa \
  -o /tmp/molcheck_rna -p out --common /tmp/rna_common.json \
  --sequence_roles examples/mss_builder/example.sequence_roles.tsv > /tmp/molcheck_rna.log 2>&1
echo "RC=$?"
grep -n "ff_definition" /tmp/molcheck_rna/out.ann
```
Expected: ff_definition 行が `... RNA, ...` になっている（mol_type=genomic RNA が反映）。

（このタスクはコミット不要。動作確認のみ。example が無い場合は Step を読み替えて手元のファイルで確認。）

---

## Verification Strategy

- **単体（Task 1）**: `_molecule_token` の判定表12ケースを parametrize で網羅。
- **経路A（Task 2）**: `ff_definition()` で mol_type 違い（DNA/RNA/mRNA）と既定DNA・chromosome経路を確認。
- **呼び出し側（Task 3）**: import OK と全テスト回帰で、必須引数化が既存呼び出しを壊していないことを確認。
- **経路B-1/B-2（Task 4）**: `create_source_feature` の通常経路とメタ経路で RNA 反映・既定 DNA を確認。
- **結合・後方互換（Task 5）**: 実 example で mss_builder を走らせ、mol_type なし→DNA、genomic RNA→RNA を目視。
- **全体**: 最後に `PYTHONPATH=src $PY -m pytest -q` で全 39 件 pass。

## Risks

1. **必須引数化で呼び出し側が壊れる**: Task 3 で2箇所を同時修正し、Task 3 Step 4/5 の import + 全テストで担保。`ff_definition` の呼び出しは grep 済みで2箇所のみ。
2. **`base_source` 変数名の差異**: mss_builder 側で source 辞書の変数名が `base_source` でない可能性。Task 3 Step 3 の grep で確認し、実際の変数に合わせる（common.SOURCE 由来の dict から mol_type を取れればよい）。
3. **mol_type 変数の二重用途**: `create_source_feature` / `_create_source_with_meta` では `mol_type` 変数が source 行（`mol_type` qualifier の値）にも使われる。`mol`（トークン）と `mol_type`（生の値）を取り違えないこと。Task 4 の注記参照。
4. **既存 example の mol_type**: example の common JSON が将来 RNA を持つと出力が変わるが、現状は DNA 系のため後方互換は保たれる（Task 5 で確認）。
