# III-V MQW Design

InP 基板上の InGaAsP / AlGaInAs 多重量子井戸（MQW）活性層を一次設計するための Python スクリプトです。
O-band、特に 1.31 um 近傍の SOA 活性層を想定し、組成、ひずみ、量子閉じ込め準位、遷移波長、ひずみバランスを概算します。
主な使用対象は InGaAsP MQW で、既定設定も InGaAsP になっています。

物理モデルの詳細は [docs/explanation.md](docs/explanation.md) にまとめています。

## できること

- InP 基板上で指定ひずみになる InGaAsP / AlGaInAs 組成の算出
- (001) 擬似格子整合層の面内ひずみ・面直ひずみの計算
- 変形ポテンシャルによる HH / LH ひずみ込みバンドギャップの概算
- 一次元有効質量有限井戸による e1 / hh1 / lh1 閉じ込め準位の計算
- e1-hh1 / e1-lh1 遷移波長の概算
- MQW 全体の平均ひずみ、ひずみ厚み積、簡易臨界膜厚の見積もり
- Lumerical MQW 計算へ渡すための `.lsf` 入力断片の生成

## 注意

このプログラムは一次設計用の簡易ツールです。
校正済みの 6x6 / 8x8 k.p シミュレータ、自己無撞着 Poisson-Schrödinger 計算、実測 PL / XRD によるエピ校正を置き換えるものではありません。

出力される遷移波長は「ひずみと量子閉じ込めを含めたバンド端遷移の概算値」です。
最終的な利得ピーク、PL ピーク、TE/TM 利得差は Lumerical MQW、nextnano などで再検証してください。

## 必要環境

- Python 3.12 以上
- uv

依存パッケージは `pyproject.toml` に記載されています。

## セットアップ

```bash
uv sync
```

## 実行例

既定では InGaAsP 系の O-band 向け候補を計算します。

```bash
uv run python -B src/BasicMQWDesign.py
```

実行すると、既定では次のファイルが `out/` フォルダ内に生成されます。

- `out/ingaasp_design_result.json`
- `out/ingaasp_lumerical_input.lsf`

AlGaInAs 系を計算する場合:

```bash
uv run python -B src/BasicMQWDesign.py --family algainas
```

出力ファイル名を指定する場合:

```bash
uv run python -B src/BasicMQWDesign.py \
  --family ingaasp \
  --json out/ingaasp_custom_design.json \
  --lsf out/ingaasp_custom_lumerical_input.lsf
```

井戸数、井戸幅、障壁幅を変更する場合:

```bash
uv run python -B src/BasicMQWDesign.py \
  --wells 7 \
  --well-nm 6.5 \
  --barrier-nm 9.0
```

ひずみや組成を指定する場合:

```bash
uv run python -B src/BasicMQWDesign.py \
  --family ingaasp \
  --well-strain -0.006 \
  --barrier-strain 0.003 \
  --as-well 0.567 \
  --as-barrier 0.30
```

AlGaInAs で Al 分率を指定する場合:

```bash
uv run python -B src/BasicMQWDesign.py \
  --family algainas \
  --al-well 0.14 \
  --al-barrier 0.30
```


## 臨界膜圧(膜応力)の計算スクリプト

別ファイルとして `src/CriticalFilmStress.py` を用意しています。
InP 基板上の単層について、次を算出します。

- 面内ひずみ `eps=(a_sub-a_layer)/a_layer`
- 二軸弾性近似の膜応力 [GPa]
- Matthews-Blakeslee の簡易臨界膜厚 [nm]

実行例 (InGaAsP):

```bash
uv run python -B src/CriticalFilmStress.py \
  --family ingaasp \
  --strain -0.006 \
  --as-frac 0.567 \
  --thickness-nm 7.0 \
  --json out/critical_film_stress.json
```

実行例 (AlGaInAs):

```bash
uv run python -B src/CriticalFilmStress.py \
  --family algainas \
  --strain -0.007 \
  --al-frac 0.14
```

## 主な引数

| 引数 | 意味 | 既定値 |
| --- | --- | --- |
| `--family` | 材料系。`algainas` または `ingaasp` | `ingaasp` |
| `--wells` | 井戸数 | `5` |
| `--well-nm` | 井戸幅 [nm] | `7.0` |
| `--barrier-nm` | 障壁幅 [nm] | `10.0` |
| `--qc` | 伝導帯オフセット比 | 材料系ごとの既定値 |
| `--well-strain` | 井戸層の目標面内ひずみ | 材料系ごとの既定値 |
| `--barrier-strain` | 障壁層の目標面内ひずみ | 材料系ごとの既定値 |
| `--al-well` | AlGaInAs 井戸層の Al 分率 | `0.14` |
| `--al-barrier` | AlGaInAs 障壁層の Al 分率 | `0.30` |
| `--as-well` | InGaAsP 井戸層の As 分率 | `0.567` |
| `--as-barrier` | InGaAsP 障壁層の As 分率 | `0.30` |
| `--json` | JSON 出力先 | `out/ingaasp_design_result.json` |
| `--lsf` | Lumerical script 出力先 | `out/ingaasp_lumerical_input.lsf` |

ひずみの符号は次の定義です。

```text
strain = (a_substrate - a_layer) / a_layer
```

このため、負の値は圧縮ひずみ、正の値は引張ひずみです。

## 出力される主な情報

標準出力には、設計候補の要約が表示されます。

```text
family             : 材料系
wells              : 井戸数
well/barrier       : 井戸幅 / 障壁幅
well material      : 井戸層組成
barrier material   : 障壁層組成
well strain        : 井戸層ひずみ
barrier strain     : 障壁層ひずみ
avg stack strain   : MQW 全体の平均ひずみ
DeltaEc/DeltaEv_hh : 伝導帯 / HH 価電子帯側の井戸深さ
e1/hh1/lh1         : 閉じ込め準位
e1-hh1 wavelength  : e1-hh1 遷移波長
e1-lh1 wavelength  : e1-lh1 遷移波長
```

JSON には、材料パラメータ、ひずみ、遷移エネルギー、遷移波長、ひずみバランス、簡易臨界膜厚などが保存されます。

LSF には、Lumerical MQW 計算に渡すための材料構造、層厚、ひずみ、初期キャリア密度掃引などが出力されます。

出力先フォルダが存在しない場合は自動作成されます。

## 既定設計の目安

現時点の既定値では、おおよそ次のような e1-hh1 遷移波長になります。

| 材料系 | 井戸 / 障壁 | e1-hh1 遷移波長の目安 |
| --- | --- | --- |
| `algainas` | 7.0 nm / 10.0 nm | 約 1310 nm |
| `ingaasp` | 7.0 nm / 10.0 nm | 約 1310 nm |

数値は材料パラメータ、バンドオフセット比、井戸幅、ひずみ設定に強く依存します。

## 詳細説明

計算している物理量、近似、式、限界については次を参照してください。

- [docs/explanation.md](docs/explanation.md)
