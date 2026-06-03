# source flag 型 qualifier の真偽判定 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** source フィーチャーの flag 型 qualifier（`environmental_sample` 等）を、TSV/JSON のどちらからでも真偽値で「付与する／しない」指定できるようにする（`false`/`no`/空欄 → 非付与、それ以外 → 値なし flag 行で付与）。

**Architecture:** 判定ロジックを `src/common/source_builder.py` の既存ヘルパー `_source_qualifier_rows` に集約する。flag 型 qualifier のホワイトリスト `_FLAG_QUALIFIERS` と真偽解釈関数 `_flag_is_set` を追加し、`_source_qualifier_rows` が flag 型なら真偽で 1 行/0 行に分岐、非 flag 型は従来通り（値そのまま、list は複数行）。per-entry 経路と COMMON/meta 経路は既にこのヘルパー経由なので、ヘルパー変更だけで両経路に反映される。入力型は `models.py` の `SOURCE` に `bool` を追加。

**Tech Stack:** Python 3.10+, pydantic v2, pytest

---

## Context & Constraints

- 設計ドキュメント: `docs/superpowers/specs/2026-06-03-source-flag-qualifier-design.md`
- 影響ツール: mss_builder / egapx2mss / batch_wgs_builder（いずれも `create_source_feature` を共有）
- **テスト実行環境:** システム `python3`(3.9) には Bio/pydantic が無い。依存入りは conda 環境 `/Users/tanizawa/miniforge3/envs/dr_tools/bin/python3`（ただし pytest 未導入）。本計画では **Task 1 で pytest を導入**し、以降 `PYTHON=/Users/tanizawa/miniforge3/envs/dr_tools/bin/python3` を使い `PYTHONPATH=src $PYTHON -m pytest ...` で実行する。
- 既存ヘルパー（実装済み・現状）:
  ```python
  # src/common/source_builder.py:116-124
  def _source_qualifier_rows(key: str, value: "str | list[str]") -> list[Row]:
      if isinstance(value, list):
          return [["", "", "", key, str(v)] for v in value]
      return [["", "", "", key, str(value)]]
  ```
  per-entry ループ（`create_source_feature` 内）と COMMON/meta ループ（`_create_source_with_meta` 内）が両方ともこれを呼ぶ。
- 現 `models.py`: `SOURCE: Optional[dict[str, str | list[str]]] = None`
- `Row = list[str]`（`source_builder.py:9`）

### File Structure

| ファイル | 役割 | 操作 |
|---|---|---|
| `src/common/source_builder.py` | `_FLAG_QUALIFIERS`・`_flag_is_set` 追加、`_source_qualifier_rows` を flag 対応に | 変更 |
| `src/common/models.py` | `SOURCE` 値型に `bool` 追加 | 変更 |
| `tests/test_source_flag_qualifier.py` | flag 判定の単体・結合テスト | 新規 |
| `CLAUDE.md` | flag qualifier の真偽記法を追記 | 変更 |

---

## Task 1: テスト環境の確認と pytest 導入

**Files:** （環境準備のみ）

- [ ] **Step 1: conda 環境に pytest を導入**

Run:
```bash
/Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pip install pytest
```
Expected: `Successfully installed pytest-...`（既に入っていれば `Requirement already satisfied`）。

- [ ] **Step 2: 既存テストが収集できることを確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest -q
```
Expected: collected 0 items（テストはまだ無い）。エラーなく終了すること。`pyproject.toml` に `[tool.pytest.ini_options] testpaths = ["tests"]` が既にある。

（コミット不要。環境準備のみ。）

---

## Task 2: `_flag_is_set` 真偽判定関数（TDD）

**Files:**
- Modify: `src/common/source_builder.py`
- Test: `tests/test_source_flag_qualifier.py`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_source_flag_qualifier.py` を新規作成:

```python
"""Tests for source flag-qualifier truthiness handling."""

from common.source_builder import _flag_is_set


def test_flag_is_set_false_strings():
    for v in ["false", "False", "FALSE", "no", "No", "NO"]:
        assert _flag_is_set(v) is False, f"{v!r} should be off"


def test_flag_is_set_false_bool():
    assert _flag_is_set(False) is False


def test_flag_is_set_true_bool():
    assert _flag_is_set(True) is True


def test_flag_is_set_truthy_strings():
    for v in ["yes", "Yes", "true", "True", "1", "x", "ATCC"]:
        assert _flag_is_set(v) is True, f"{v!r} should be on"


def test_flag_is_set_empty_string_is_on():
    # backward compatibility: JSON "environmental_sample": "" means on
    assert _flag_is_set("") is True


def test_flag_is_set_whitespace_false():
    assert _flag_is_set("  no  ") is False
```

- [ ] **Step 2: テストが失敗することを確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_source_flag_qualifier.py -q
```
Expected: FAIL（`ImportError: cannot import name '_flag_is_set'`）。

- [ ] **Step 3: 最小実装**

`src/common/source_builder.py` の `_source_qualifier_rows` 定義（116行付近）の**直前**に追加:

```python
# INSDC source qualifiers written without a value (the bare key).
_FLAG_QUALIFIERS = frozenset({
    "environmental_sample",
    "transgenic",
    "germline",
    "rearranged",
    "proviral",
    "macronuclear",
    "metagenomic",
    "focus",
})

# Values interpreted as "off" for a flag qualifier (case-insensitive).
_FALSE_VALUES = frozenset({"false", "no"})


def _flag_is_set(value) -> bool:
    """Interpret a flag-qualifier value as on/off.

    Off: boolean False, or the strings "false"/"no" (case-insensitive,
    surrounding whitespace ignored).
    On:  everything else — boolean True, "" (empty string, kept on for
    backward compatibility with ``"environmental_sample": ""``), "true",
    "yes", "1", and any other non-empty string.
    """
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() not in _FALSE_VALUES
```

- [ ] **Step 4: テストが通ることを確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_source_flag_qualifier.py -q
```
Expected: 6 passed。

- [ ] **Step 5: コミット**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/common/source_builder.py tests/test_source_flag_qualifier.py
git commit -m "feat: add _flag_is_set for source flag qualifiers

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: `_source_qualifier_rows` を flag 対応に（TDD）

**Files:**
- Modify: `src/common/source_builder.py:116-124`（`_source_qualifier_rows` 本体）
- Test: `tests/test_source_flag_qualifier.py`（追記）

- [ ] **Step 1: 失敗するテストを追記**

`tests/test_source_flag_qualifier.py` の末尾に追加:

```python
from common.source_builder import _source_qualifier_rows


def test_rows_flag_on_emits_valueless_row():
    assert _source_qualifier_rows("environmental_sample", "yes") == [
        ["", "", "", "environmental_sample", ""]
    ]


def test_rows_flag_empty_string_on():
    assert _source_qualifier_rows("environmental_sample", "") == [
        ["", "", "", "environmental_sample", ""]
    ]


def test_rows_flag_off_emits_nothing():
    assert _source_qualifier_rows("environmental_sample", "no") == []
    assert _source_qualifier_rows("transgenic", False) == []


def test_rows_nonflag_keeps_value():
    # a non-flag qualifier whose value happens to be "No" is NOT treated as a flag
    assert _source_qualifier_rows("strain", "No") == [
        ["", "", "", "strain", "No"]
    ]


def test_rows_nonflag_list_multiple():
    assert _source_qualifier_rows("culture_collection", ["ATCC:1", "NBRC:2"]) == [
        ["", "", "", "culture_collection", "ATCC:1"],
        ["", "", "", "culture_collection", "NBRC:2"],
    ]
```

- [ ] **Step 2: テストが失敗することを確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_source_flag_qualifier.py -q
```
Expected: 新規5件が FAIL（`environmental_sample` が値付き行を返す／`no` で空にならない）。

- [ ] **Step 3: `_source_qualifier_rows` を置き換え**

`src/common/source_builder.py` の現 `_source_qualifier_rows`（116-124行）を以下に置き換え:

```python
def _source_qualifier_rows(key: str, value) -> list[Row]:
    """Emit source qualifier row(s) for one (key, value).

    Flag qualifiers (in ``_FLAG_QUALIFIERS``) emit a single valueless row
    when the value is truthy (see :func:`_flag_is_set`) and are omitted
    entirely when falsy. Other qualifiers emit one row per value — a list
    yields one row per element (e.g. multiple ``culture_collection``).
    """
    if key in _FLAG_QUALIFIERS:
        return [["", "", "", key, ""]] if _flag_is_set(value) else []
    if isinstance(value, list):
        return [["", "", "", key, str(v)] for v in value]
    return [["", "", "", key, str(value)]]
```

- [ ] **Step 4: テストが通ることを確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_source_flag_qualifier.py -q
```
Expected: 11 passed（Task2 の 6 + 今回 5）。

- [ ] **Step 5: コミット**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/common/source_builder.py tests/test_source_flag_qualifier.py
git commit -m "feat: emit flag source qualifiers by truthiness

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: `SOURCE` 入力型に `bool` を追加（TDD）

**Files:**
- Modify: `src/common/models.py`（`SOURCE` フィールド定義）
- Test: `tests/test_source_flag_qualifier.py`（追記）

- [ ] **Step 1: 失敗するテストを追記**

`tests/test_source_flag_qualifier.py` の末尾に追加:

```python
from common.models import CommonModel


def test_source_accepts_bool_value():
    data = {
        "DBLINK": {"project": "P", "sample": "S"},
        "SOURCE": {"organism": "E. coli", "environmental_sample": False},
    }
    m = CommonModel.model_validate(data)
    assert m.SOURCE["environmental_sample"] is False


def test_source_accepts_list_value_still():
    data = {
        "DBLINK": {"project": "P", "sample": "S"},
        "SOURCE": {"organism": "E. coli", "culture_collection": ["ATCC:1", "NBRC:2"]},
    }
    m = CommonModel.model_validate(data)
    assert m.SOURCE["culture_collection"] == ["ATCC:1", "NBRC:2"]
```

- [ ] **Step 2: テストが失敗することを確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_source_flag_qualifier.py -k bool -q
```
Expected: `test_source_accepts_bool_value` が FAIL（pydantic ValidationError: bool is not a valid string）。

- [ ] **Step 3: `SOURCE` の型を広げる**

`src/common/models.py` の `SOURCE` 行を変更:

```python
# 変更前（現行）
    SOURCE: Optional[dict[str, str | list[str]]] = None
# 変更後
    SOURCE: Optional[dict[str, str | bool | list[str]]] = None
```

- [ ] **Step 4: テストが通ることを確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_source_flag_qualifier.py -q
```
Expected: 13 passed。

- [ ] **Step 5: コミット**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/common/models.py tests/test_source_flag_qualifier.py
git commit -m "feat: accept bool values for SOURCE qualifiers

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: 結合テスト（`create_source_feature` 経由）

**Files:**
- Test: `tests/test_source_flag_qualifier.py`（追記）

- [ ] **Step 1: 結合テストを追記**

`create_source_feature` は per-entry / meta 両経路を持つ。`environmental_sample` はカテゴリ自動付与の影響を避けるため、auto_source_qualifiers を持たない素のカテゴリ（`GNM`）で検証する。`tests/test_source_flag_qualifier.py` 末尾に追加:

```python
from common.source_builder import create_source_feature


def _quals(rows):
    """(key, value) pairs from 5-column source rows, skipping the feature line."""
    return [(r[3], r[4]) for r in rows if r[3]]


def test_create_source_flag_off_omitted_per_entry():
    src = {"organism": "E. coli", "environmental_sample": "no"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    keys = [k for k, _ in _quals(rows)]
    assert "environmental_sample" not in keys
    assert "organism" in keys


def test_create_source_flag_on_valueless_per_entry():
    src = {"organism": "E. coli", "environmental_sample": "yes"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    assert ("environmental_sample", "") in _quals(rows)


def test_create_source_nonflag_no_value_kept():
    # strain="No" must remain a normal valued qualifier, not be dropped
    src = {"organism": "E. coli", "strain": "No"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    assert ("strain", "No") in _quals(rows)


def test_create_source_flag_on_meta_path():
    src = {"organism": "E. coli", "transgenic": True}
    rows = create_source_feature("WGS", None, None, None, src,
                                 source_modifier_key="strain",
                                 use_meta_expression=True)
    assert ("transgenic", "") in _quals(rows)


def test_create_source_flag_off_meta_path():
    src = {"organism": "E. coli", "transgenic": False}
    rows = create_source_feature("WGS", None, None, None, src,
                                 source_modifier_key="strain",
                                 use_meta_expression=True)
    keys = [k for k, _ in _quals(rows)]
    assert "transgenic" not in keys
```

- [ ] **Step 2: テストが通ることを確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest tests/test_source_flag_qualifier.py -q
```
Expected: 18 passed。

これは Task 2-4 の実装が両経路に効いていることの確認。新規実装は無いが、もし FAIL した場合は実装側（`_source_qualifier_rows`）の問題なので、Task 3 の Step 3 を見直す。

- [ ] **Step 3: コミット**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add tests/test_source_flag_qualifier.py
git commit -m "test: integration tests for source flag qualifiers via create_source_feature

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: ドキュメント追記

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: CLAUDE.md に flag qualifier の記法を追記**

`CLAUDE.md` の「## common_example.json の形式」節の末尾（`末尾カンマ (JSON5スタイル) は許容される。` の直後）に以下を追記:

```markdown

### source の flag 型 qualifier（値なし qualifier）

`environmental_sample` など値を持たない qualifier（INSDC flag 型）は、真偽値で付与を制御できる。

- 付与しない: `false` / `no`（大文字小文字問わず）、JSON boolean `false`、または **TSV の空欄**
- 付与する（値なし行 `/qualifier` として出力）: 上記以外。**推奨は `yes`**。JSON では `true`（boolean）や空文字 `""`（後方互換）も付与扱い。

対象 flag: `environmental_sample`, `transgenic`, `germline`, `rearranged`, `proviral`, `macronuclear`, `metagenomic`, `focus`。
flag 型でない通常 qualifier（`strain` 等）はこの判定の対象外で、値はそのまま出力される（例: `strain` の値 `No` は文字列 `No` のまま）。

注意: ENV/MAG/MAG-WGS など `environmental_sample` を自動付与するカテゴリでは、`false` を指定しても自動付与が優先される。
```

- [ ] **Step 2: 追記内容を目視確認**

Run:
```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git --no-pager diff CLAUDE.md
```
Expected: common_example.json 節の後ろに上記ブロックのみが追加されている。

- [ ] **Step 3: コミット**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add CLAUDE.md
git commit -m "docs: document source flag qualifier truthiness

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Verification Strategy

- **単体（Task 2-4）:** `_flag_is_set` の真偽表、`_source_qualifier_rows` の flag on/off・非 flag・list、`SOURCE` の bool/list 受理を pytest で確認。
- **結合（Task 5）:** `create_source_feature` の per-entry 経路・COMMON/meta 経路の両方で flag on→値なし行・off→行なし、非 flag（`strain=No`）が誤判定されないことを確認。
- **全体:** 最後に `PYTHONPATH=src /Users/tanizawa/miniforge3/envs/dr_tools/bin/python3 -m pytest -q` で全 18 件 pass を確認。
- **回帰:** 既存の biosample/culture_collection 対応（list 値）はテスト `test_rows_nonflag_list_multiple` / `test_source_accepts_list_value_still` でカバー。

## Risks

1. **flag リストの過不足**: ホワイトリストに無い flag 型 qualifier は従来通り値ありで出る。今回は INSDC 標準8件で固定。将来の追加は `_FLAG_QUALIFIERS` への追記で対応（スコープ外）。
2. **カテゴリ自動付与との優先関係**: `auto_source_qualifiers` カテゴリでは `false` 指定が無効。これは仕様（spec のスコープ外）で、Task 6 のドキュメントで明示。`environmental_sample` の結合テストは自動付与を持たない `GNM` カテゴリで行い干渉を避ける。
3. **pytest 未導入**: Task 1 で導入。CI 等が別環境の場合は要再確認だが、本リポジトリにCI設定は無い。
