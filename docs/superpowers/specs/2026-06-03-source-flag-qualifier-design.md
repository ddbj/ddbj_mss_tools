# source flag 型 qualifier の真偽判定 設計

- 日付: 2026-06-03
- 対象: `src/common/source_builder.py`（および入力型の `src/common/models.py`）
- 影響ツール: mss_builder / egapx2mss / batch_wgs_builder（いずれも `create_source_feature` を共有）

## 背景・課題

`environmental_sample` のような「値を持たない（flag 型）」source qualifier の扱いに、入力経路ごとの不統一がある。

- **JSON 経路**（mss_builder / egapx2mss の `--common`）: `"environmental_sample": ""` のようにキーを書けば値なし行として出力される。「キーがあれば付与」という規約で、値の中身は実質見ていない。
- **TSV 経路**（batch_wgs_builder の sample_list）: `core.py:61-62` が空セルを `continue` でスキップするため、空欄は source_dict に入らず非付与。逆に何か書くと flag 型でも値付き（例 `environmental_sample\tyes`）の不正な行が出る。

つまり TSV から「値なし flag を、意図して付ける／外す」を表現する手段がなく、JSON とも整合していない。

## 要件

1. flag 型 qualifier を TSV / JSON のどちらからでも、真偽値で「付与する／しない」を指定できる。
2. 付与時は必ず**値なし行** `["", "", "", key, ""]` で出力する（`key\tyes` のような不正値を出さない）。
3. 既存の JSON `"environmental_sample": ""`（空文字＝付与）を壊さない（後方互換重視）。
4. flag 型でない通常 qualifier（`strain`, `culture_collection` 等）には影響しない。`strain` の値が `No` でも文字列 `No` としてそのまま出力する（誤判定回避）。
5. 判定ロジックは1箇所に集約し、per-entry 経路と COMMON/meta 経路で食い違わない。

## 設計

### 1. flag 型 qualifier のホワイトリスト

`source_builder.py` にモジュール定数として定義する。INSDC で source フィーチャーに値なしで記載しうる qualifier を中心に、将来の feature 拡張も見据えて広めに持つ。

```python
# INSDC qualifiers that take no value (written as the bare key).
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
```

このリストに**ある** qualifier のみ真偽判定の対象。リストに無い qualifier は従来通り（`_source_qualifier_rows` で値そのまま、list 値は複数行）。

### 2. 真偽判定関数

flag 型 qualifier の値を真偽に解釈する単一関数を `source_builder.py` に置く。

```python
_FALSE_VALUES = frozenset({"false", "no"})

def _flag_is_set(value) -> bool:
    """Interpret a flag-qualifier value as on/off.

    Off: boolean False, or the strings "false"/"no" (case-insensitive).
    On:  everything else, including "" (empty string), "true", "yes", "1",
         boolean True, and any other non-empty string.
    """
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() not in _FALSE_VALUES
```

判定表:

| 入力値 | 結果 |
|---|---|
| 文字列 `false` / `no`（大文字小文字無視） | off（付与しない） |
| boolean `False` | off |
| boolean `True` | on（付与） |
| `true` / `yes` / `1` / その他任意の非空文字列 | on |
| `""`（空文字） | **on（付与）— 後方互換** |

### 3. qualifier 行生成の分岐

`_source_qualifier_rows(key, value)` を、flag 型と通常型で分岐するよう拡張する。

```python
def _source_qualifier_rows(key: str, value) -> list[Row]:
    """Emit source qualifier row(s) for one (key, value).

    Flag qualifiers (in _FLAG_QUALIFIERS) are emitted as a single valueless
    row when the value is truthy, and omitted entirely when falsy.
    Other qualifiers emit one row per value (a list yields multiple rows).
    """
    if key in _FLAG_QUALIFIERS:
        return [["", "", "", key, ""]] if _flag_is_set(value) else []
    if isinstance(value, list):
        return [["", "", "", key, str(v)] for v in value]
    return [["", "", "", key, str(value)]]
```

per-entry 経路（`create_source_feature` 内のループ、現 189-191 行）と COMMON/meta 経路（`_create_source_with_meta` 内のループ、現 232-234 行）は既に `_source_qualifier_rows` を呼んでいるため、**この関数の変更だけで両経路に反映される**（呼び出し側の変更不要）。

### 4. 入力型の拡張（boolean を許可）

JSON で `true`/`false`（boolean リテラル）を書けるよう SOURCE の値型を広げる。

```python
# models.py（現: dict[str, str | list[str]]）
SOURCE: Optional[dict[str, str | bool | list[str]]] = None
```

TSV 経路は値が常に文字列（`core.py:61` で `str(...).strip()`）なので boolean 対応は不要だが、`_flag_is_set` は文字列・boolean 両方を受けるため両経路で同じ関数が使える。

### 5. TSV 空欄の扱い（変更なし）

`core.py:61-62` の空セルスキップはそのまま。TSV 空欄＝未記入＝非付与。`false`/`no` を明示的に書けば（flag 型なら）非付与にできる。

## 入力経路まとめ

| 経路 / 入力 | 結果 |
|---|---|
| JSON `"environmental_sample": ""` | 付与（後方互換）|
| JSON `"environmental_sample": "yes"` | 付与（値は捨て、値なし行）|
| JSON `"environmental_sample": "no"` | 付与しない（新）|
| JSON `"environmental_sample": true` | 付与（新）|
| JSON `"environmental_sample": false` | 付与しない（新）|
| TSV 空欄 | 付与しない（core.py スキップ）|
| TSV `yes` | 付与 |
| TSV `No` / `false` | 付与しない（新）|
| 通常 qualifier（`strain` 等）の値 `No` | 文字列 `No` として出力（flag 判定対象外）|

## スコープ外

- **カテゴリ自動付与との関係**: ENV/MAG/MAG-WGS 等は `auto_source_qualifiers` で `environmental_sample` を無条件付与する（`source_builder.py:155,185 / 217,230`）。この経路はユーザー入力より優先され、`false` 指定があっても自動付与される。本変更では auto 機構は変えず、ドキュメントに「自動付与カテゴリでは false 指定は効かない」と注記するに留める。
- flag リストの外部設定化（JSON 化）はしない（YAGNI。必要になれば別途）。

## ドキュメント

- `CLAUDE.md` の common JSON 説明、および batch_wgs_builder の TSV 説明に、flag 型 qualifier の真偽記法を追記。
- 推奨記法: 付与は `yes`、非付与は `no`。
- 注記: 「`false`/`no`（大文字小文字無視）と空欄のみ非付与、それ以外は付与」「flag 型でない qualifier には影響しない」「自動付与カテゴリでは false は効かない」。

## テスト方針

`tests/` に `source_builder` の単体テストを追加（pytest）:
- `_flag_is_set`: `false`/`no`/`False`/`FALSE`/`No` → False、`""`/`yes`/`true`/`True`/`1`/任意文字列 → True。
- `_source_qualifier_rows`: flag 型 on→1行（値空）、off→0行、通常型→値そのまま、list→複数行。
- `create_source_feature` 結合: `environmental_sample: "no"` で行が出ない／`"yes"` で値なし行が出る。`strain: "No"` は `No` 値のまま出る（誤判定しない）。
