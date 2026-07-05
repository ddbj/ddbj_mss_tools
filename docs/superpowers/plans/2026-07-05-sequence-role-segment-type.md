# sequence role `type=segment` 対応 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `--sequence_roles` TSV の `type` に `segment` を追加し、submission 内の segment 数と `status` に応じて source の `/segment` qualifier と ff_definition を生成する。

**Architecture:** 判定ロジックは共有モジュール `src/common/source_builder.py` の純粋関数 `ff_definition()` / `source_qualifier()` に `segment_count` 引数を足して実装する。`src/mss_builder/ann_writer.py` と `src/egapx2mss/ann_writer.py` は既存の `chromosome_count` と同様に `segment_count` を集計して両関数へ渡す。既存 `chromosome`（count 依存）と `plasmid`（status 依存）の組み合わせ挙動に相当する。

**Tech Stack:** Python 3.10+, pytest, pydantic v2, BioPython。

## Global Constraints

- Python 3.10 以上（`match` 文・`X | Y` 型ヒント使用可）。
- 新引数 `segment_count` は必ずデフォルト値 `0` を持たせ、既存呼び出し・既存テストを壊さないこと。
- `type` にホワイトリスト検証は存在しない。検証コードは追加しない。
- テストは PYTHONPATH に `src` を通して実行する（例: `PYTHONPATH=src pytest ...`）。既存 CI/ローカルの慣習に従う。
- 全コミットの末尾に以下のトレーラを付ける:
  ```
  Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01KsvgqPHeLSD8bXwvoSf7BD
  ```
- 作業ブランチは `feat/sequence-role-segment-type`（既に spec コミット済み）。

---

### Task 1: `ff_definition` の segment 分岐

**Files:**
- Modify: `src/common/source_builder.py`（`ff_definition` 関数、現状 `def ff_definition(...)` 128 行付近〜171 行付近）
- Test: `tests/test_ff_definition_molecule.py`

**Interfaces:**
- Consumes: 既存 `SequenceRoleEntry(seq_id, type_, seq_name, status, is_circular)`、`_molecule_token()`。
- Produces: `ff_definition(entry, seq_id, organism, infraspecific_name_modifier, mol_type, is_wgs=False, chromosome_count=0, segment_count=0) -> str`（`segment_count` を新たに追加）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_ff_definition_molecule.py` の末尾（plasmid セクションの後）に追記:

```python
# ── segment: count + status 依存 ─────────────────────────────────────────
def test_ff_segment_single_complete_genome():
    e = SequenceRoleEntry("seg1", "segment", "", "complete", False)
    out = ff_definition(e, "seg1", "Influenza A virus", "", "viral cRNA",
                        is_wgs=False, segment_count=1)
    assert out == "Influenza A virus RNA, complete genome"


def test_ff_segment_single_partial_genome():
    e = SequenceRoleEntry("seg1", "segment", "", "partial", False)
    out = ff_definition(e, "seg1", "Influenza A virus", "", "viral cRNA",
                        is_wgs=False, segment_count=1)
    assert out == "Influenza A virus RNA, partial genome"


def test_ff_segment_multi_complete_sequence():
    e = SequenceRoleEntry("seg4", "segment", "4", "complete", False)
    out = ff_definition(e, "seg4", "Influenza A virus", "isolate X", "viral cRNA",
                        is_wgs=False, segment_count=8)
    assert out == "Influenza A virus isolate X RNA, segment 4, complete sequence"


def test_ff_segment_multi_partial_sequence():
    e = SequenceRoleEntry("seg4", "segment", "4", "partial", False)
    out = ff_definition(e, "seg4", "Influenza A virus", "", "viral cRNA",
                        is_wgs=False, segment_count=8)
    assert out == "Influenza A virus RNA, segment 4, partial sequence"


def test_ff_segment_multi_empty_name_fallback():
    e = SequenceRoleEntry("seg1", "segment", "", "complete", False)
    out = ff_definition(e, "seg1", "Influenza A virus", "", "viral cRNA",
                        is_wgs=False, segment_count=2)
    assert out == "Influenza A virus RNA, segment, complete sequence"
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `PYTHONPATH=src pytest tests/test_ff_definition_molecule.py -k segment -v`
Expected: FAIL（`segment_count` を渡すと `TypeError: unexpected keyword argument`、または現状は fallback に落ちて `"Influenza A virus RNA, seg1"` になり assert 不一致）

- [ ] **Step 3: 最小実装**

`src/common/source_builder.py` の `ff_definition` シグネチャに `segment_count` を追加:

```python
def ff_definition(entry: Optional[SequenceRoleEntry], seq_id: str, organism: str,
                  infraspecific_name_modifier: str, mol_type: str, is_wgs: bool = False,
                  chromosome_count: int = 0, segment_count: int = 0) -> str:
```

`plasmid` 分岐の直後、`# fallback` コメントの直前に segment 分岐を挿入:

```python
    if entry.type == "segment":
        if segment_count <= 1:
            if entry.status == "complete":
                return f"{prefix} {mol}, complete genome"
            return f"{prefix} {mol}, partial genome"
        seg_part = f"segment {entry.seq_name}".strip() if entry.seq_name else "segment"
        if entry.status == "complete":
            return f"{prefix} {mol}, {seg_part}, complete sequence"
        return f"{prefix} {mol}, {seg_part}, partial sequence"
```

同関数の docstring の `chromosome_count` 説明の後に 1 行追記:

```python
    *segment_count* is the number of segment-type entries in the whole submission;
    a single segment uses 'complete/partial genome' (no 'segment' word), otherwise
    'segment {seq_name}, complete/partial sequence'.
```

- [ ] **Step 4: テストが通ることを確認**

Run: `PYTHONPATH=src pytest tests/test_ff_definition_molecule.py -k segment -v`
Expected: PASS（5 件）

- [ ] **Step 5: 回帰確認**

Run: `PYTHONPATH=src pytest tests/test_ff_definition_molecule.py -v`
Expected: PASS（既存 + 新規すべて）

- [ ] **Step 6: コミット**

```bash
git add src/common/source_builder.py tests/test_ff_definition_molecule.py
git commit -m "feat(source): ff_definition supports type=segment"
```

---

### Task 2: `source_qualifier` の segment 分岐

**Files:**
- Modify: `src/common/source_builder.py`（`source_qualifier` 関数 58〜78 行付近、および `SequenceRoleEntry` の type コメント 20 行）
- Test: `tests/test_ff_definition_molecule.py`

**Interfaces:**
- Produces: `source_qualifier(entry, seq_id, is_wgs=False, segment_count=0) -> dict[str, str]`（`segment_count` を新たに追加）。segment・単一 → `{}`、segment・複数かつ `seq_name` 非空 → `{"segment": seq_name}`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_ff_definition_molecule.py` の末尾（source_qualifier セクション付近）に追記:

```python
def test_source_qualifier_segment_single_omitted():
    e = SequenceRoleEntry("seg1", "segment", "", "complete", False)
    assert source_qualifier(e, "seg1", segment_count=1) == {}


def test_source_qualifier_segment_single_name_ignored():
    e = SequenceRoleEntry("seg1", "segment", "4", "complete", False)
    assert source_qualifier(e, "seg1", segment_count=1) == {}


def test_source_qualifier_segment_multi():
    e = SequenceRoleEntry("seg4", "segment", "4", "complete", False)
    assert source_qualifier(e, "seg4", segment_count=8) == {"segment": "4"}


def test_source_qualifier_segment_multi_empty_name_omitted():
    e = SequenceRoleEntry("seg4", "segment", "", "complete", False)
    assert source_qualifier(e, "seg4", segment_count=8) == {}
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `PYTHONPATH=src pytest tests/test_ff_definition_molecule.py -k "source_qualifier_segment" -v`
Expected: FAIL（`TypeError: unexpected keyword argument 'segment_count'`）

- [ ] **Step 3: 最小実装**

`source_qualifier` シグネチャに `segment_count` を追加:

```python
def source_qualifier(entry: Optional[SequenceRoleEntry], seq_id: str,
                     is_wgs: bool = False, segment_count: int = 0) -> dict[str, str]:
```

`plasmid` 分岐の直後、最後の `return {}` の直前に segment 分岐を挿入:

```python
    if entry.type == "segment":
        if segment_count >= 2 and entry.seq_name:
            return {"segment": entry.seq_name}
        return {}
```

同関数の docstring のルール列挙に 1 行追記:

```python
    - segment: single -> no qualifier; multiple -> segment = seq_name (omitted when empty)
```

`SequenceRoleEntry` の type コメント（20 行付近）を更新:

```python
        self.type = type_            # "chromosome" | "organelle" | "plasmid" | "segment" | "unplaced"
```

- [ ] **Step 4: テストが通ることを確認**

Run: `PYTHONPATH=src pytest tests/test_ff_definition_molecule.py -k "source_qualifier_segment" -v`
Expected: PASS（4 件）

- [ ] **Step 5: コミット**

```bash
git add src/common/source_builder.py tests/test_ff_definition_molecule.py
git commit -m "feat(source): source_qualifier emits /segment for multi-segment"
```

---

### Task 3: 両 ann_writer に `segment_count` を配線 + 統合テスト

**Files:**
- Modify: `src/mss_builder/ann_writer.py`（`chromosome_count` 集計 125〜129 行付近、呼び出し 182〜187 行付近）
- Modify: `src/egapx2mss/ann_writer.py`（`chromosome_count` 集計 155〜159 行付近、呼び出し 186〜191 行付近）
- Test: `tests/test_mss_segment.py`（新規作成）

**Interfaces:**
- Consumes: Task 1/2 の `ff_definition(..., segment_count=...)` と `source_qualifier(..., segment_count=...)`、既存 `write_mss_ann(fsa_path, ann_path, common=None, sequence_roles=None, submission_category=None)`、`CommonModel`、`SequenceRoleEntry`。

- [ ] **Step 1: 失敗する統合テストを書く**

`tests/test_mss_segment.py` を新規作成:

```python
"""End-to-end: type=segment in sequence roles -> mss_builder .ann output."""

from common.models import CommonModel
from common.source_builder import SequenceRoleEntry
from mss_builder.ann_writer import write_mss_ann


def _common():
    return CommonModel.model_validate({
        "DBLINK": {"project": "P", "sample": "S"},
        "SOURCE": {"organism": "Influenza A virus", "mol_type": "viral cRNA"},
    })


def test_multi_segment_emits_segment_qualifier(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">seg4\nACGTACGTAC\n>seg6\nTTTTGGGGCC\n")
    ann = tmp_path / "out.ann"
    roles = {
        "seg4": SequenceRoleEntry("seg4", "segment", "4", "complete", False),
        "seg6": SequenceRoleEntry("seg6", "segment", "6", "complete", False),
    }
    write_mss_ann(str(fasta), str(ann), common=_common(), sequence_roles=roles)
    text = ann.read_text()
    assert "\t\t\tsegment\t4\n" in text
    assert "Influenza A virus RNA, segment 4, complete sequence" in text


def test_single_segment_omits_qualifier_and_uses_genome(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">seg1\nACGTACGTAC\n")
    ann = tmp_path / "out.ann"
    roles = {"seg1": SequenceRoleEntry("seg1", "segment", "", "complete", False)}
    write_mss_ann(str(fasta), str(ann), common=_common(), sequence_roles=roles)
    text = ann.read_text()
    assert "Influenza A virus RNA, complete genome" in text
    assert "\t\t\tsegment\t" not in text
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `PYTHONPATH=src pytest tests/test_mss_segment.py -v`
Expected: FAIL（`multi` 側で `/segment` 行が出ず、ff_definition が `... RNA, seg4`（fallback）になる。ann_writer が `segment_count` をまだ渡していないため）

- [ ] **Step 3: mss_builder ann_writer を実装**

`src/mss_builder/ann_writer.py` の `chromosome_count` 集計ブロック（125〜129 行付近）の直後に追加:

```python
    segment_count = 0
    if sequence_roles:
        segment_count = sum(
            1 for e in sequence_roles.values() if e.type == "segment"
        )
```

呼び出し（182〜187 行付近）を変更:

```python
            source_quals.update(source_qualifier(role_entry, entry_id, is_wgs=False, segment_count=segment_count))
            source_quals["ff_definition"] = ff_definition(
                role_entry, entry_id, organism, infraspecific_name_modifier,
                base_source.get("mol_type", ""), is_wgs=False,
                chromosome_count=chromosome_count, segment_count=segment_count,
            )
```

- [ ] **Step 4: egapx2mss ann_writer を実装**

`src/egapx2mss/ann_writer.py` の `chromosome_count` 集計ブロック（155〜159 行付近）の直後に追加:

```python
    segment_count = 0
    if sequence_roles:
        segment_count = sum(
            1 for e in sequence_roles.values() if e.type == "segment"
        )
```

呼び出し（186〜191 行付近）を変更:

```python
        source_quals.update(source_qualifier(role_entry, entry_id, is_wgs, segment_count=segment_count))
        source_quals["ff_definition"] = ff_definition(
            role_entry, entry_id, organism, infraspecific_name_modifier,
            base_source.get("mol_type", ""), is_wgs,
            chromosome_count=chromosome_count, segment_count=segment_count,
        )
```

- [ ] **Step 5: 統合テストが通ることを確認**

Run: `PYTHONPATH=src pytest tests/test_mss_segment.py -v`
Expected: PASS（2 件）

- [ ] **Step 6: 全体回帰確認**

Run: `PYTHONPATH=src pytest -q`
Expected: PASS（`slow` マーク以外の全テスト。既存の egapx2mss/mss_builder テストが壊れていないこと）

- [ ] **Step 7: コミット**

```bash
git add src/mss_builder/ann_writer.py src/egapx2mss/ann_writer.py tests/test_mss_segment.py
git commit -m "feat(mss): wire segment_count into mss_builder and egapx2mss"
```

---

### Task 4: ドキュメント更新

**Files:**
- Modify: `CLAUDE.md`（sequence role の type 一覧、ff_definition 表）
- Modify: `README.md`（sequence role 説明）

**Interfaces:** なし（ドキュメントのみ）。

- [ ] **Step 1: CLAUDE.md を更新**

「sequence role (`--sequence_roles` TSV) と ff_definition」節の `type` の説明文（`chromosome` / `organelle` / `plasmid` / `unplaced` のいずれか）に `segment` を追加し、`| type | status | ff_definition |` 表の `unplaced` 行の前に以下 2 行を追加:

```markdown
| segment（submission 全体で1件のみ） | complete/partial | `{prefix} {mol}, complete genome` / `{prefix} {mol}, partial genome`（`/segment` は付与しない） |
| segment（複数） | complete/partial | `{prefix} {mol}, segment {seq_name}, complete sequence` / `... partial sequence`（source に `/segment` を付与） |
```

- [ ] **Step 2: README.md を更新**

README.md の sequence role / `type` を説明している箇所（`chromosome` / `organelle` / `plasmid` / `unplaced` を列挙している行。日本語版・英語版の両方があれば両方）に `segment` を追加し、単一=`complete genome`（`/segment` なし）・複数=`segment {seq_name}, complete sequence`（`/segment` あり）である旨を 1〜2 文で追記する。

- [ ] **Step 3: 記述の齟齬がないか確認**

Run: `grep -n "segment" CLAUDE.md README.md`
Expected: 追加した行が表示され、単一/複数の挙動が spec と一致していること。

- [ ] **Step 4: コミット**

```bash
git add CLAUDE.md README.md
git commit -m "docs: document type=segment in sequence roles"
```

---

## Self-Review

**1. Spec coverage:**
- source `/segment`（単一なし・複数あり）→ Task 2 + Task 3。
- ff_definition 4 パターン → Task 1 + Task 3（統合）。
- segment_count 集計を両 ann_writer に配線 → Task 3。
- topology/circular は既存挙動のまま（変更なし）→ 変更不要と spec に明記済み、対応タスク不要。
- is_wgs で segment を placed 扱い → 既存 `_is_unplaced` が `unplaced` 以外を placed 扱いするため変更不要（Task 3 の統合テストが非 WGS 経路を通ることで担保）。
- ドキュメント（CLAUDE.md / README.md）→ Task 4。
- 全カバー、gap なし。

**2. Placeholder scan:** "TBD"/"適切に"/"等を処理" 等の曖昧語なし。README のみ既存文面を grep して合わせる指示だが、追記内容（単一/複数の挙動）は明示済み。

**3. Type consistency:**
- `segment_count` の名称・デフォルト `0` は Task 1/2/3 で一貫。
- `ff_definition(..., chromosome_count=0, segment_count=0)` と `source_qualifier(..., is_wgs=False, segment_count=0)` のシグネチャが Task と統合テスト・両 ann_writer 呼び出しで一致。
- `SequenceRoleEntry(seq_id, type_, seq_name, status, is_circular)` の位置引数順がテスト全体で一致。
