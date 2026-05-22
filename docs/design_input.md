# MQW 設計入力ファイル（design-input JSON）

MQW の縦構造パラメータを JSON ファイルで指定し、各スクリプトに `--design-input` オプションで渡せます。

## 概要

従来は CLI 引数（`--as-well 0.567 --well-strain -0.006 ...`）で個々のパラメータを指定していましたが、
**組成をモル分率で直接指定** し、ひずみは自動算出させることができます。

## JSON スキーマ

### InGaAsP（In₁₋ₓGaₓAsᵧP₁₋ᵧ）

```json
{
  "family": "ingaasp",
  "wells": 5,
  "well_nm": 7.0,
  "barrier_nm": 10.0,
  "well": {
    "y_As": 0.567,
    "x_Ga": 0.1755
  },
  "barrier": {
    "y_As": 0.30,
    "x_Ga": 0.1797
  }
}
```

### AlGaInAs（Al_x Ga_y In₁₋ₓ₋ᵧ As）

```json
{
  "family": "algainas",
  "wells": 5,
  "well_nm": 7.0,
  "barrier_nm": 10.0,
  "well": {
    "x_Al": 0.14,
    "y_Ga": 0.34
  },
  "barrier": {
    "x_Al": 0.30,
    "y_Ga": 0.24
  }
}
```

## フィールド一覧

### トップレベル

| キー | 型 | 必須 | 説明 |
|---|---|---|---|
| `family` | string | いいえ | `"ingaasp"` または `"algainas"`。省略時 `"ingaasp"` |
| `wells` | int | いいえ | 量子井戸数。省略時 5 |
| `well_nm` | float | いいえ | 井戸幅 [nm]。省略時 7.0 |
| `barrier_nm` | float | いいえ | 障壁幅 [nm]。省略時 10.0 |
| `qc` | float | いいえ | 伝導帯オフセット比 ΔEc/ΔEg。省略時はファミリーデフォルト |
| `eg_offset_well_eV` | float | いいえ | 井戸バンドギャップ補正 [eV]。省略時 0.0 |
| `eg_offset_barrier_eV` | float | いいえ | 障壁バンドギャップ補正 [eV]。省略時 0.0 |
| `well` | object | いいえ | 井戸層の組成（下記参照） |
| `barrier` | object | いいえ | 障壁層の組成（下記参照） |

### InGaAsP の well / barrier ブロック

| キー | 型 | 説明 |
|---|---|---|
| `y_As` | float | As モル分率 $y$。省略時はファミリーデフォルト |
| `x_Ga` | float | Ga モル分率 $x$。省略時はひずみターゲットから逆算 |
| `strain` | float | ひずみターゲット $\varepsilon_\parallel$。`x_Ga` と排他（`x_Ga` 指定時は無視） |

### AlGaInAs の well / barrier ブロック

| キー | 型 | 説明 |
|---|---|---|
| `x_Al` | float | Al モル分率 $x$。省略時はファミリーデフォルト |
| `y_Ga` | float | Ga モル分率 $y$。省略時はひずみターゲットから逆算 |
| `strain` | float | ひずみターゲット $\varepsilon_\parallel$。`y_Ga` と排他（`y_Ga` 指定時は無視） |

## 組成指定モード

| モード | 指定する項目 | ひずみ | 典型用途 |
|---|---|---|---|
| **全モル分率指定** | `y_As` + `x_Ga`（または `x_Al` + `y_Ga`） | 自動算出 | エピ成長レシピの検証 |
| **片方 + ひずみ** | `y_As` + `strain`（または `x_Al` + `strain`） | 指定値 | 従来どおりの設計探索 |
| **片方のみ** | `y_As` のみ（または `x_Al` のみ） | デフォルト値 | 最小入力 |

## 使用例

### BasicMQWDesign.py

```bash
uv run python -B src/BasicMQWDesign.py \
  --design-input my_mqw.json \
  --json out/design.json \
  --lsf out/design.lsf
```

### MQWGainDesign.py

```bash
uv run python -B src/MQWGainDesign.py \
  --design-input my_mqw.json \
  --calibration calibrations/ingaasp_oband_example.json
```

### MQWGainSweep.py

キャリア密度または線幅の sweep と組み合わせ可能です。

```bash
uv run python -B src/MQWGainSweep.py \
  --design-input my_mqw.json \
  --sweep carrier-density \
  --start 1e18 --stop 3e18 --points 5
```

> **注意**: `--design-input` 使用時の sweep は `carrier-density` と `broadening-eV` のみ対応です。
> 井戸幅やひずみを掃引する場合は、従来の CLI 引数方式を使ってください。

## `--design-input` と `--design-json` の違い

| | `--design-input` | `--design-json` |
|---|---|---|
| ファイル内容 | 設計**入力パラメータ**（組成・層厚） | 計算済み **DesignDict** 全体 |
| 処理 | `design_default()` を呼んで構造を計算 | 計算済み結果をそのまま使用 |
| ファイルサイズ | 小（10〜15 行） | 大（数百行） |
| 用途 | 構造候補の管理・入力 | 同一構造の再計算高速化 |

## サンプルファイル

`tests/fixtures/` にサンプルが用意されています。

| ファイル | 内容 |
|---|---|
| `tests/fixtures/ingaasp_oband_input.json` | InGaAsP O-band（~1.31 µm） |
| `tests/fixtures/ingaasp_cband_input.json` | InGaAsP C-band（~1.55 µm） |
| `tests/fixtures/suzuki2018_cband_input.json` | Suzuki, JJAP 57, 094101 (2018) Table I の C-band MQW |
