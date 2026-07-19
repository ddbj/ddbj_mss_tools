# ff_definition メタ記法化 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `ff_definition()`（sequence role の per-entry パス）を具体値展開から DDBJ MSS の `@@[...]@@` メタ記法テンプレート生成に変更する。

**Architecture:** 共有 `src/common/source_builder.py` の `ff_definition()` のシグネチャとロジックを変更（`organism`/`modifier` 値を廃止し `source_identifier` 名を受け取る、seq_name をメタ参照化、chromosome を count 分岐、count≥2 chromosome/segment と plasmid で seq_name 空を `ValueError`）。3 呼び出し側（mss_builder/egapx2mss/gff2mss）を新シグネチャに合わせ、不要になった値算出を削除。テスト2ファイルとドキュメント2ファイルを更新。

**Tech Stack:** Python 3.10+, pytest, pydantic v2, BioPython。

## Global Constraints

- Python 3.10 以上。
- テスト実行: `PYTHONPATH=src pytest ...`。`pytest` が PATH に無ければ `PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest ...` を使う（このenvにpytestあり）。gff2mss 関連テストは `ddbj_gff` 未導入により collection error になる既知事象（23件）で、本変更とは無関係。回帰は `-m "not slow" --continue-on-collection-errors` で `78 passed`（＋既存23 errors）を基準とする。
- `source_qualifier()` は変更しない。`create_source_feature()`（WGS/COMMON メタ経路）も変更しない。
- 新 `ff_definition` シグネチャ: `ff_definition(entry, source_identifier, mol_type, is_wgs=False, chromosome_count=0, segment_count=0) -> str`。
- prefix: `source_identifier` 非空→`@@[organism]@@ @@[{source_identifier}]@@`、空/None→`@@[organism]@@`。`{mol}` は `_molecule_token(mol_type)` の具体値。
- seq_name 必須（空なら `ValueError`）: plasmid（count 問わず）／chromosome count≥2／segment count≥2。単一 chromosome/segment は seq_name 空を許容。
- コミット末尾トレーラ:
  ```
  Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01KsvgqPHeLSD8bXwvoSf7BD
  ```
- 作業ブランチ `feat/ff-definition-meta-notation`（spec コミット済み）。

---

### Task 1: ff_definition のメタ記法化（関数＋3 caller＋テスト）

**Files:**
- Modify: `src/common/source_builder.py`（`ff_definition` 関数 133–190 行）
- Modify: `src/mss_builder/ann_writer.py`（146・151 行の未使用変数削除、189–193 行の呼び出し）
- Modify: `src/egapx2mss/ann_writer.py`（140・144 行の未使用変数削除、193–197 行の呼び出し）
- Modify: `src/gff2mss/assemble.py`（47・49 行の未使用変数削除、79 行の呼び出し）
- Test: `tests/test_ff_definition_molecule.py`（ff_definition を呼ぶテスト群を差し替え）
- Test: `tests/test_mss_segment.py`（統合アサーション 26・36 行）

**Interfaces:**
- Consumes: 既存 `SequenceRoleEntry(seq_id, type_, seq_name, status, is_circular)`、`_molecule_token`、`_organelle_code`。各 writer が保持する `source_id_key`/`src_id_key`（qualifier 名、None 可）。
- Produces: `ff_definition(entry, source_identifier, mol_type, is_wgs=False, chromosome_count=0, segment_count=0) -> str`（旧 `seq_id`/`organism`/`infraspecific_name_modifier` 引数を廃止）。

- [ ] **Step 1: ff_definition の単体テストを差し替え（失敗する新テスト）**

`tests/test_ff_definition_molecule.py` から、`ff_definition(` を呼ぶ**旧テスト関数をすべて削除**する（下記の名前。`_molecule_token`/`_organelle_code`/`create_source_*`/`source_qualifier_*` のテストは残す）:
`test_ff_definition_unplaced_wgs_dna`, `test_ff_definition_unplaced_wgs_rna`, `test_ff_definition_unplaced_wgs_mrna`, `test_ff_definition_chromosome_complete_rna`, `test_ff_definition_empty_mol_type_defaults_dna`, `test_ff_chromosome_single_complete_genome`, `test_ff_chromosome_multi_complete_sequence`, `test_ff_chromosome_default_count_complete_sequence`, `test_ff_chromosome_partial_no_localization_suffix`, `test_ff_organelle_mitochondrion_complete`, `test_ff_organelle_partial_genome`, `test_ff_organelle_plastid_chloroplast`, `test_ff_organelle_passthrough_name`, `test_ff_organelle_rna_token_position`, `test_ff_plasmid_complete`, `test_ff_plasmid_partial`, `test_ff_segment_single_complete_genome`, `test_ff_segment_single_partial_genome`, `test_ff_segment_multi_complete_sequence`, `test_ff_segment_multi_partial`, `test_ff_segment_multi_empty_name_fallback`。

そのうえで、ファイル末尾に次の新テストを追加する（`ff_definition` と `SequenceRoleEntry` と `pytest` は既に import 済み）:

```python
# ── ff_definition: meta-notation (new signature) ───────────────────────
def _entry(type_, seq_name="", status="complete"):
    return SequenceRoleEntry("sid", type_, seq_name, status, False)


@pytest.mark.parametrize("entry,source_identifier,mol_type,is_wgs,chrom,seg,expected", [
    # unplaced / None
    (None, None, "genomic DNA", True, 0, 0,
     "@@[organism]@@ DNA, @@[submitter_seqid]@@"),
    (None, None, "genomic DNA", False, 0, 0,
     "@@[organism]@@ DNA, unplaced sequence @@[entry]@@"),
    (None, "cultivar", "genomic RNA", False, 0, 0,
     "@@[organism]@@ @@[cultivar]@@ RNA, unplaced sequence @@[entry]@@"),
    (_entry("unplaced"), "strain", "genomic DNA", True, 0, 0,
     "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"),
    # chromosome — single (count<=1): no number
    (_entry("chromosome", "1", "complete"), "strain", "genomic DNA", False, 1, 0,
     "@@[organism]@@ @@[strain]@@ DNA, chromosome, complete genome"),
    (_entry("chromosome", "1", "partial"), "strain", "genomic DNA", False, 1, 0,
     "@@[organism]@@ @@[strain]@@ DNA, chromosome"),
    (_entry("chromosome", "", "complete"), None, "genomic DNA", False, 1, 0,
     "@@[organism]@@ DNA, chromosome, complete genome"),
    # chromosome — multiple (count>=2): @@[chromosome]@@
    (_entry("chromosome", "1", "complete"), "strain", "genomic DNA", False, 2, 0,
     "@@[organism]@@ @@[strain]@@ DNA, chromosome @@[chromosome]@@, complete sequence"),
    (_entry("chromosome", "1", "partial"), "strain", "genomic DNA", False, 2, 0,
     "@@[organism]@@ @@[strain]@@ DNA, chromosome @@[chromosome]@@"),
    # organelle — prefix meta, adjective concrete
    (_entry("organelle", "mitochondrion", "complete"), "", "genomic DNA", False, 0, 0,
     "@@[organism]@@ mitochondrial DNA, complete genome"),
    (_entry("organelle", "mitochondrion", "partial"), None, "genomic DNA", False, 0, 0,
     "@@[organism]@@ mitochondrial DNA, partial genome"),
    (_entry("organelle", "plastid:chloroplast", "complete"), "isolate", "genomic DNA", False, 0, 0,
     "@@[organism]@@ @@[isolate]@@ chloroplast DNA, complete genome"),
    (_entry("organelle", "mitochondrion", "complete"), None, "genomic RNA", False, 0, 0,
     "@@[organism]@@ mitochondrial RNA, complete genome"),
    # plasmid
    (_entry("plasmid", "pLG1", "complete"), "strain", "genomic DNA", False, 0, 0,
     "@@[organism]@@ @@[strain]@@ plasmid @@[plasmid]@@ DNA, complete sequence"),
    (_entry("plasmid", "pLG1", "partial"), "strain", "genomic DNA", False, 0, 0,
     "@@[organism]@@ @@[strain]@@ plasmid @@[plasmid]@@ DNA, partial sequence"),
    # segment — single (count<=1): no 'segment' word
    (_entry("segment", "", "complete"), "strain", "viral cRNA", False, 0, 1,
     "@@[organism]@@ @@[strain]@@ RNA, complete genome"),
    (_entry("segment", "", "partial"), "strain", "viral cRNA", False, 0, 1,
     "@@[organism]@@ @@[strain]@@ RNA, partial genome"),
    # segment — multiple (count>=2): @@[segment]@@
    (_entry("segment", "4", "complete"), "strain", "viral cRNA", False, 0, 8,
     "@@[organism]@@ @@[strain]@@ RNA, segment @@[segment]@@, complete sequence"),
    (_entry("segment", "4", "partial"), "strain", "viral cRNA", False, 0, 8,
     "@@[organism]@@ @@[strain]@@ RNA, segment @@[segment]@@"),
    # fallback (unknown type)
    (_entry("weird", "", "complete"), None, "genomic DNA", False, 0, 0,
     "@@[organism]@@ DNA, @@[entry]@@"),
    # mol token via mRNA
    (None, None, "mRNA", True, 0, 0,
     "@@[organism]@@ mRNA, @@[submitter_seqid]@@"),
    (None, None, "", True, 0, 0,
     "@@[organism]@@ DNA, @@[submitter_seqid]@@"),
])
def test_ff_definition_meta(entry, source_identifier, mol_type, is_wgs, chrom, seg, expected):
    out = ff_definition(entry, source_identifier, mol_type, is_wgs=is_wgs,
                        chromosome_count=chrom, segment_count=seg)
    assert out == expected


@pytest.mark.parametrize("entry,chrom,seg", [
    (_entry("plasmid", "", "complete"), 0, 0),
    (_entry("plasmid", "", "partial"), 0, 0),
    (_entry("chromosome", "", "complete"), 2, 0),
    (_entry("chromosome", "", "partial"), 2, 0),
    (_entry("segment", "", "complete"), 0, 2),
    (_entry("segment", "", "partial"), 0, 2),
])
def test_ff_definition_empty_seqname_raises(entry, chrom, seg):
    with pytest.raises(ValueError):
        ff_definition(entry, "strain", "genomic DNA", is_wgs=False,
                      chromosome_count=chrom, segment_count=seg)


def test_ff_definition_single_empty_seqname_allowed():
    # single chromosome / single segment with empty seq_name must NOT raise
    assert ff_definition(_entry("chromosome", "", "complete"), None, "genomic DNA",
                         chromosome_count=1) == "@@[organism]@@ DNA, chromosome, complete genome"
    assert ff_definition(_entry("segment", "", "complete"), None, "genomic DNA",
                         segment_count=1) == "@@[organism]@@ DNA, complete genome"
```

- [ ] **Step 2: 新テストが失敗することを確認（RED）**

Run: `PYTHONPATH=src pytest tests/test_ff_definition_molecule.py -k "meta or empty_seqname or single_empty" -q`
Expected: FAIL（旧 `ff_definition` は位置引数 `seq_id, organism, ...` を取るため、新シグネチャ呼び出しで `TypeError`、または期待メタ文字列と不一致）

- [ ] **Step 3: ff_definition を新シグネチャ・メタ記法に書き換え**

`src/common/source_builder.py` の `ff_definition`（133–190 行）を丸ごと以下に置換:

```python
def ff_definition(entry: Optional[SequenceRoleEntry], source_identifier: Optional[str],
                  mol_type: str, is_wgs: bool = False,
                  chromosome_count: int = 0, segment_count: int = 0) -> str:
    """
    Build the ff_definition qualifier value as a DDBJ MSS @@[...]@@ meta-notation
    template. Values are substituted by MSS at submission time from the source
    feature's own qualifiers (/organism, the SOURCE_IDENTIFIER qualifier,
    /chromosome, /plasmid, /segment, /submitter_seqid) or MSS-provided @@[entry]@@.

    *source_identifier* is the NAME of the SOURCE_IDENTIFIER qualifier
    (e.g. "cultivar", "strain", "isolate"); None/empty omits the modifier ref.
    *is_wgs* is True when all entries in the submission are unplaced (WGS mode).
    *chromosome_count* / *segment_count* are the counts of chromosome- / segment-type
    entries in the whole submission (single vs multiple changes the wording).

    Raises ValueError when a required seq_name is empty: plasmid (always), and
    chromosome / segment when their count >= 2.
    """
    if source_identifier:
        prefix = f"@@[organism]@@ @@[{source_identifier}]@@"
    else:
        prefix = "@@[organism]@@"
    mol = _molecule_token(mol_type)

    if entry is None or entry.type == "unplaced":
        if is_wgs:
            return f"{prefix} {mol}, @@[submitter_seqid]@@"
        return f"{prefix} {mol}, unplaced sequence @@[entry]@@"

    if entry.type == "chromosome":
        if chromosome_count >= 2:
            if not entry.seq_name:
                raise ValueError("chromosome entry requires a non-empty seq_name when count >= 2")
            if entry.status == "complete":
                return f"{prefix} {mol}, chromosome @@[chromosome]@@, complete sequence"
            return f"{prefix} {mol}, chromosome @@[chromosome]@@"
        # single chromosome (count <= 1): no number
        if entry.status == "complete":
            return f"{prefix} {mol}, chromosome, complete genome"
        return f"{prefix} {mol}, chromosome"

    if entry.type == "organelle":
        converted = _organelle_code(entry.seq_name)
        if entry.status == "complete":
            return f"{prefix} {converted} {mol}, complete genome"
        return f"{prefix} {converted} {mol}, partial genome"

    if entry.type == "plasmid":
        if not entry.seq_name:
            raise ValueError("plasmid entry requires a non-empty seq_name")
        if entry.status == "complete":
            return f"{prefix} plasmid @@[plasmid]@@ {mol}, complete sequence"
        return f"{prefix} plasmid @@[plasmid]@@ {mol}, partial sequence"

    if entry.type == "segment":
        if segment_count >= 2:
            if not entry.seq_name:
                raise ValueError("segment entry requires a non-empty seq_name when count >= 2")
            if entry.status == "complete":
                return f"{prefix} {mol}, segment @@[segment]@@, complete sequence"
            return f"{prefix} {mol}, segment @@[segment]@@"
        # single segment (count <= 1): no 'segment' word
        if entry.status == "complete":
            return f"{prefix} {mol}, complete genome"
        return f"{prefix} {mol}, partial genome"

    # fallback (unknown type)
    return f"{prefix} {mol}, @@[entry]@@"
```

- [ ] **Step 4: 単体テストが通ることを確認（GREEN）**

Run: `PYTHONPATH=src pytest tests/test_ff_definition_molecule.py -q`
Expected: PASS（この時点で caller は未修正のため、他ファイルは壊れていてよい。このファイル単体で全 PASS を確認）

- [ ] **Step 5: mss_builder の呼び出しを新シグネチャに更新＋未使用変数削除**

`src/mss_builder/ann_writer.py`:
1. 146 行 `    organism = base_source.get("organism", "")` を**削除**。
2. 151 行 `    infraspecific_name_modifier = base_source.get(source_id_key, "") if source_id_key else ""` を**削除**（`source_id_key` の算出 147–150 行は残す）。
3. 189–193 行の呼び出しを置換:

```python
            source_quals["ff_definition"] = ff_definition(
                role_entry, source_id_key, base_source.get("mol_type", ""),
                is_wgs=False, chromosome_count=chromosome_count, segment_count=segment_count,
            )
```

- [ ] **Step 6: egapx2mss の呼び出しを更新＋未使用変数削除**

`src/egapx2mss/ann_writer.py`:
1. 140 行 `    organism = base_source.get("organism", "")` を**削除**。
2. 144 行 `    infraspecific_name_modifier = base_source.get(source_id_key, "") if source_id_key else ""` を**削除**（`source_id_key` の算出 141–143 行は残す）。139 行のコメント `# infraspecific_name_modifier: ...` も削除。
3. 193–197 行の呼び出しを置換:

```python
        source_quals["ff_definition"] = ff_definition(
            role_entry, source_id_key, base_source.get("mol_type", ""),
            is_wgs, chromosome_count=chromosome_count, segment_count=segment_count,
        )
```

- [ ] **Step 7: gff2mss の呼び出しを更新＋未使用変数削除**

`src/gff2mss/assemble.py`:
1. 47 行 `    organism = base_source.get("organism", "")` を**削除**。
2. 49 行 `    infra = base_source.get(src_id_key, "") if src_id_key else ""` を**削除**（48 行 `src_id_key = common.SOURCE_IDENTIFIER` は残す）。
3. 79 行の呼び出しを置換（gff2mss は従来どおり chromosome_count を渡さない＝既存仕様）:

```python
        src["ff_definition"] = ff_definition(role, src_id_key, mol_type, is_wgs, segment_count=segment_count)
```

- [ ] **Step 8: 統合テスト（test_mss_segment.py）のアサーションをメタ記法に更新**

`tests/test_mss_segment.py` の 26 行:
```python
    assert "@@[organism]@@ RNA, segment @@[segment]@@, complete sequence" in text
```
36 行:
```python
    assert "@@[organism]@@ RNA, complete genome" in text
```
（common に SOURCE_IDENTIFIER が無いため prefix は `@@[organism]@@` のみ。mol_type=viral cRNA→RNA。`/segment` 行の検証 25 行・37 行は変更しない。）

- [ ] **Step 9: 全体回帰を確認**

Run: `PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest -q -m "not slow" --continue-on-collection-errors`
Expected: `78 passed`（＋既存 23 errors は `ddbj_gff` 未導入によるもので無関係）。egapx2mss/mss_builder 系テストが壊れていないこと。

- [ ] **Step 10: コミット**

```bash
git add src/common/source_builder.py src/mss_builder/ann_writer.py src/egapx2mss/ann_writer.py src/gff2mss/assemble.py tests/test_ff_definition_molecule.py tests/test_mss_segment.py
git commit -m "feat(source): ff_definition emits @@[...]@@ meta-notation"
```

---

### Task 2: ドキュメント更新（CLAUDE.md / README.md）

**Files:**
- Modify: `CLAUDE.md`（「sequence role と ff_definition」節の表・説明）
- Modify: `README.md`（sequence role 節の ff_definition 記述、日・英）

**Interfaces:** なし（ドキュメントのみ）。

- [ ] **Step 1: CLAUDE.md の ff_definition 表をメタ記法版に差し替え**

`CLAUDE.md` の該当節（`| type | status | ff_definition |` 表と前後の説明）を、spec `docs/superpowers/specs/2026-07-19-ff-definition-meta-notation.md` の「決定テーブル」および prefix/mol 説明に合わせて書き換える。具体的には:
- prefix の説明を「`@@[organism]@@ @@[{source_identifier}]@@`（source_identifier 空なら `@@[organism]@@` のみ）」に。
- 表を #1〜#14（chromosome は count==1 と count≥2 を分離、`@@[chromosome]@@`/`@@[plasmid]@@`/`@@[segment]@@` を使用、unplaced は is_wgs で `@@[submitter_seqid]@@` / `unplaced sequence @@[entry]@@`）に更新。
- 「count≥2 chromosome/segment と plasmid は seq_name 必須（空は ValueError）」の注記を追加。
- organelle の形容詞形は変換値のまま（メタ化しない）旨を明記。

- [ ] **Step 2: README.md の sequence role 節（日・英）を更新**

`README.md` の sequence role / ff_definition を説明している箇所（日本語・英語の両方）を、上記メタ記法・count 分岐・seq_name 必須ルールに合わせて更新する。`grep -n "ff_definition\|chromosome\|segment\|unplaced" README.md` で該当箇所を特定して修正する。

- [ ] **Step 3: 記述の整合を確認**

Run: `grep -nE "@@\[|ff_definition|seq_name" CLAUDE.md README.md`
Expected: メタ記法の記述が spec と一致し、旧具体値表記（`{organism} {識別子}` や `unlocalized sequence` 等）が残っていないこと。

- [ ] **Step 4: コミット**

```bash
git add CLAUDE.md README.md
git commit -m "docs: document ff_definition meta-notation"
```

---

## Self-Review

**1. Spec coverage:**
- 新シグネチャ（source_identifier 名、seq_id/organism/modifier 削除）→ Task 1 Step 3・5–7。
- prefix メタ（空 source_id で @@[organism]@@ のみ）→ Step 3、テスト（Step 1 の該当ケース）。
- 決定テーブル #1〜#14（chromosome count 分岐、unplaced is_wgs 分岐、fallback @@[entry]@@）→ Step 3、Step 1 パラメタライズ。
- seq_name 必須の ValueError（plasmid / chromosome≥2 / segment≥2）→ Step 3、Step 1 `test_ff_definition_empty_seqname_raises`。単一許容 → `test_ff_definition_single_empty_seqname_allowed`。
- organelle 形容詞据え置き → Step 3、テストの organelle ケース。
- 3 caller 更新＋未使用変数削除 → Step 5–7。
- 統合テスト更新 → Step 8。回帰 → Step 9。
- ドキュメント → Task 2。
- source_qualifier/create_source_feature 不変 → いずれの Step でも変更しない（明記）。
- 全カバー、gap なし。

**2. Placeholder scan:** Task 1 は全コードを提示。Task 2 のみ「grep して該当箇所を更新」という指示だが、更新内容（prefix/表/ValueError 注記/organelle）を spec 参照で明示しており、対象は既存の限定的な節。曖昧語なし。

**3. Type consistency:**
- 新シグネチャ `ff_definition(entry, source_identifier, mol_type, is_wgs=False, chromosome_count=0, segment_count=0)` が関数定義（Step 3）・3 caller（Step 5–7）・テスト（Step 1）で一致。
- `source_id_key`（mss_builder/egapx2mss）/`src_id_key`（gff2mss）を第2引数に渡す点が各 caller で一致。
- `SequenceRoleEntry(seq_id, type_, seq_name, status, is_circular)` の位置引数順がテストヘルパ `_entry` と一致。
