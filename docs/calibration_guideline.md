# k.p 材料ゲイン計算の校正ガイドライン

この文書は、`src/MQWGainDesign.py` の簡易 k.p 材料ゲイン計算を、Lumerical MQW、PL、吸収、modal gain などの基準データに合わせるための校正手順をまとめる。

本ツールのゲイン計算はスクリーニング用であり、絶対ゲイン値は未校正のままでは設計判断の最終値として扱わない。
まず遷移波長、次にスペクトル幅と TE/TM 傾向、最後に絶対ゲインを合わせる。

## 1. 校正対象

### 1.1 遷移波長

最初に合わせる量は e1-hh1 遷移波長または gain / absorption / PL のピーク波長である。

主に効くパラメータ:

- 井戸層・障壁層のバンドギャップ `Eg`
- 伝導帯オフセット比 `qc`
- ひずみ変形ポテンシャル `ac`, `av`, `b`
- 有効質量と Luttinger parameter
- 実効井戸幅

推奨基準データ:

- Lumerical MQW の e1-hh1 遷移または gain peak
- nextnano などの k.p 計算結果
- 低注入 PL peak
- 吸収端

注意:

- PL peak は Stokes shift、温度、注入密度の影響を受ける。
- PL だけで band offset や有効質量を一意に決めることはできない。

### 1.2 TE/TM 比

次に TE/TM の相対傾向を合わせる。
圧縮ひずみ井戸では HH が上がり、TE が優勢になることが期待される。

主に効くパラメータ:

- HH/LH splitting
- ひずみ量
- 価電子帯変形ポテンシャル `b`
- `gamma1`, `gamma2`, `gamma3`
- valence band mixing の扱い
- 井戸幅

比較対象:

- Lumerical MQW の TE/TM gain
- TE/TM absorption
- 偏波分解 PL

見るべき指標:

- `peak_TE / peak_TM`
- 同一波長での `gain_TE - gain_TM`
- キャリア密度 sweep に対する TE/TM 傾向

### 1.3 スペクトル幅

ピーク高さを合わせる前に、スペクトル幅を合わせる。

主に効くパラメータ:

- `--broadening-eV`
- 線幅関数 `--line-shape`
- キャリア密度
- 温度
- 不均一広がり

初期値の目安:

- `--broadening-eV 0.020` から `0.050`
- まずは `lorentzian`
- PL の不均一広がりが強い場合は `gaussian` も比較する

### 1.4 絶対ゲイン

最後に絶対ゲイン値を合わせる。

主に効くパラメータ:

- `--gain-scale-cm`
- `--broadening-eV`
- 遷移行列要素
- active volume の定義
- carrier density の定義
- confinement factor

実測 modal gain と比較する場合は、材料ゲインと modal gain の関係を分けて扱う。

```text
g_modal = Gamma * g_material - alpha_i
```

ここで `Gamma` は光閉じ込め係数、`alpha_i` は内部損失である。
`Gamma` と `alpha_i` が不明なまま `gain_scale_cm` を調整すると、材料ゲインではなくデバイス損失まで吸収した経験係数になる。

## 2. 推奨校正順序

### Step 1: 入力構造を完全に揃える

基準データとツール側で、以下を一致させる。

- 材料系
- 井戸数
- 井戸幅
- 障壁幅
- 井戸組成
- 障壁組成
- 井戸ひずみ
- 障壁ひずみ
- 温度
- キャリア密度
- 線幅
- TE/TM の定義

Lumerical MQW と比較する場合は、同じ構造を `BasicMQWDesign.py` の `.lsf` 出力と照合する。

### Step 2: e1-hh1 遷移波長を合わせる

まず `BasicMQWDesign.py` と `MQWGainDesign.py` の低 `k_t` 遷移を確認する。

確認項目:

- e1-hh1 遷移波長
- e1-lh1 遷移波長
- 井戸幅を広げたときの長波長化
- 圧縮ひずみを強めたときの HH/LH 分裂

ずれの原因候補:

- 波長が全体に短い、または長い: `Eg` モデルを疑う
- 井戸幅依存性が合わない: 有効質量または band offset を疑う
- e1-hh1 と e1-lh1 の間隔が合わない: ひずみ、`b`、Luttinger parameter を疑う

この段階では `--gain-scale-cm` は触らない。

### Step 3: `qc` と band offset を確認する

既存実装では、伝導帯井戸深さを `qc` で近似している。
遷移波長と井戸幅依存性が合わない場合は、`qc` を sweep する。

例:

```bash
uv run python -B src/MQWGainDesign.py --qc 0.35
uv run python -B src/MQWGainDesign.py --qc 0.40
uv run python -B src/MQWGainDesign.py --qc 0.45
```

判断:

- `qc` を上げると電子閉じ込めが強くなり、電子準位が変化する。
- 1 点のピーク波長だけで `qc` を決めない。
- 複数の井戸幅または複数の組成で同じ `qc` が使えるか確認する。

### Step 4: ひずみ依存性と TE/TM 傾向を合わせる

井戸ひずみを数点振り、TE/TM の変化を見る。

例:

```bash
uv run python -B src/MQWGainDesign.py --well-strain -0.004
uv run python -B src/MQWGainDesign.py --well-strain -0.006
uv run python -B src/MQWGainDesign.py --well-strain -0.008
```

確認項目:

- 圧縮ひずみで TE が強くなるか
- e1-hh1 と e1-lh1 の分離が基準データと近いか
- TM gain が過大または過小でないか

TE/TM が合わない場合の優先確認:

1. ひずみ定義と符号
2. `b_eV`
3. `gamma1`, `gamma2`, `gamma3`
4. valence mixing モデル

### Step 5: スペクトル幅を合わせる

ピーク位置が合った後で、`--broadening-eV` を調整する。

例:

```bash
uv run python -B src/MQWGainDesign.py --broadening-eV 0.020
uv run python -B src/MQWGainDesign.py --broadening-eV 0.030
uv run python -B src/MQWGainDesign.py --broadening-eV 0.050
```

判断:

- ピーク幅が狭すぎる場合は `broadening-eV` を増やす。
- ピーク幅が広すぎる場合は `broadening-eV` を減らす。
- PL に合わせる場合は不均一広がりが入りやすい。
- gain spectrum に合わせる場合は測定条件の温度・注入密度を確認する。

### Step 6: carrier density sweep を確認する

絶対ゲイン校正前に、キャリア密度に対する傾向を確認する。

例:

```bash
uv run python -B src/MQWGainSweep.py \
  --calibration calibrations/ingaasp_oband_example.json \
  --sweep carrier-density \
  --values 1e18,1.5e18,2e18,3e18 \
  --out-json out/gain_sweep_density.json \
  --out-csv out/gain_sweep_density.csv \
  --plot out/gain_sweep_density_peak.png \
  --spectra-csv out/gain_sweep_density_spectra.csv \
  --spectra-plot out/gain_sweep_density_spectra.png
```

確認項目:

- キャリア密度を上げると gain peak が増えるか
- 透明キャリア密度の目安が基準データと極端にずれていないか
- ピーク波長のシフトが不自然でないか
- `--spectra-plot` の波長軸上で、スペクトル全体の横移動とピーク形状が自然か

### Step 7: 最後に `gain_scale_cm` を合わせる

ピーク波長、TE/TM 傾向、スペクトル幅、carrier density 依存性が概ね合った後で、`--gain-scale-cm` を調整する。

例:

```bash
uv run python -B src/MQWGainDesign.py --gain-scale-cm 1800
uv run python -B src/MQWGainDesign.py --gain-scale-cm 2400
uv run python -B src/MQWGainDesign.py --gain-scale-cm 3000
```

注意:

- 1 点の peak gain だけで合わせない。
- 複数の carrier density で同じ `gain_scale_cm` が使えるか確認する。
- 実測 modal gain と合わせる場合は `Gamma` と `alpha_i` を別途扱う。

`gain_scale_cm` の候補をまとめて比較する場合は、`MQWGainSweep.py` ではなく単発計算を複数回行う。
これは `gain_scale_cm` がピーク高さをほぼ線形に変える校正係数であり、波長依存性の sweep よりも基準データとの残差評価として扱う方がよいためである。

### Step 8: 簡易 fit で校正 JSON を生成する

ピーク波長、スペクトル幅、ピークゲインの目標値がある場合は、`FitCalibration.py` で校正 JSON を生成できる。
初期実装では、次の順序で最小限のパラメータだけを合わせる。

```text
Eg_offset_well_eV -> broadening_eV -> gain_scale_cm
```

例:

```bash
uv run python -B src/FitCalibration.py \
  --calibration-in calibrations/ingaasp_oband_example.json \
  --target-peak-wavelength-nm 1310 \
  --target-te-peak-gain-cm 1200 \
  --target-fwhm-meV 35 \
  --out calibrations/fitted/ingaasp_fit.json
```

生成された JSON には `reference` と `fit_result` が保存される。
ただし、これは単一点ターゲットに対する簡易 fit であり、複数の carrier density、井戸幅、または Lumerical / 実測スペクトル全体で再確認してから使う。

## 3. 基準データ別の使い方

### 3.1 Lumerical MQW

最初の校正基準として最も扱いやすい。

合わせる項目:

- 層構造
- 組成
- ひずみ
- 温度
- carrier density
- 線幅
- TE/TM 定義

比較順:

1. e1-hh1 / e1-lh1 遷移
2. TE/TM gain spectrum のピーク波長
3. TE/TM 比
4. carrier density sweep
5. 絶対 gain scale

### 3.2 PL

主に遷移波長の校正に使う。

注意:

- PL peak は gain peak と一致しない場合がある。
- 励起強度依存性を確認する。
- 温度を揃える。
- Stokes shift を考慮する。

PL だけで絶対ゲインを校正しない。

### 3.3 吸収スペクトル

遷移エネルギーと TE/TM 選択則の確認に有用である。

見る項目:

- absorption edge
- e1-hh1 / e1-lh1 の分離
- TE/TM absorption 比

### 3.4 実測 modal gain

最終的な絶対値校正に使う。

必要な追加情報:

- optical confinement factor `Gamma`
- internal loss `alpha_i`
- active volume
- 注入電流から carrier density への変換
- 温度
- waveguide mode の TE/TM 定義

材料ゲインへ戻して比較する場合:

```text
g_material = (g_modal + alpha_i) / Gamma
```

## 4. 最小限の校正パラメータ

最初から多くの自由度を持たせると、物理的に意味のある校正になりにくい。
初期段階では次に限定する。

```text
Eg_offset_well_eV
Eg_offset_barrier_eV
qc
broadening_eV
gain_scale_cm
```

校正 JSON では次のように指定する。

```json
{
  "name": "ingaasp_oband_example",
  "description": "Example only; not calibrated to measured data.",
  "band": {
    "qc": 0.4,
    "Eg_offset_well_eV": 0.0,
    "Eg_offset_barrier_eV": 0.0
  },
  "gain": {
    "broadening_eV": 0.03,
    "line_shape": "lorentzian",
    "gain_scale_cm": 2400.0
  }
}
```

実行例:

```bash
uv run python -B src/MQWGainDesign.py \
  --calibration calibrations/ingaasp_oband_example.json \
  --carrier-density-cm3 2e18 \
  --out-json out/gain_calibrated.json \
  --out-csv out/gain_calibrated.csv \
  --plot out/gain_calibrated.png
```

CLI 引数で同じ値を明示した場合は、CLI が校正ファイルより優先される。

```bash
uv run python -B src/MQWGainDesign.py \
  --calibration calibrations/ingaasp_oband_example.json \
  --gain-scale-cm 3000
```

この場合、`gain_scale_cm` は校正ファイル中の値ではなく `3000` が使われる。
適用された値と CLI override の有無は、出力 JSON の `calibration` に保存される。

`MQWGainSweep.py` でも同じ校正ファイルを使用できる。

```bash
uv run python -B src/MQWGainSweep.py \
  --calibration calibrations/ingaasp_oband_example.json \
  --sweep qc \
  --values 0.35,0.40,0.45 \
  --out-json out/gain_sweep_qc.json \
  --out-csv out/gain_sweep_qc.csv \
  --plot out/gain_sweep_qc_peak.png \
  --spectra-csv out/gain_sweep_qc_spectra.csv \
  --spectra-plot out/gain_sweep_qc_spectra.png
```

`--sweep qc` と `--sweep broadening-eV` では、sweep 値が校正ファイルや CLI の固定値より優先される。
sweep JSON には、校正ファイル由来の base 設定と、sweep 対象がどの校正フィールドを上書きしているかが保存される。

必要になった場合だけ追加する候補:

```text
ac_scale
av_scale
b_scale
gamma_scale
effective_well_width_nm
effective_barrier_width_nm
```

## 5. ずれ方ごとの診断

### 5.1 peak wavelength が全体にずれる

優先して確認するもの:

1. `Eg` モデル
2. 井戸組成
3. 実効井戸幅
4. `qc`

### 5.2 井戸幅依存性が合わない

優先して確認するもの:

1. 電子有効質量
2. HH/LH 有効質量
3. band offset
4. 有限井戸境界条件

### 5.3 TE/TM 比が合わない

優先して確認するもの:

1. ひずみ符号
2. `b_eV`
3. Luttinger parameter
4. valence band mixing

### 5.4 peak 幅が合わない

優先して確認するもの:

1. `broadening_eV`
2. `line_shape`
3. 温度
4. carrier density

### 5.5 peak 高さだけが合わない

最後に確認するもの:

1. `gain_scale_cm`
2. optical confinement factor
3. internal loss
4. active volume
5. carrier density 換算

## 6. 推奨する校正ログ

校正ごとに、最低限以下を記録する。

```text
date
reference_data
material_family
wells
well_nm
barrier_nm
well_composition
barrier_composition
well_strain
barrier_strain
temperature_K
carrier_density_cm3
qc
broadening_eV
line_shape
gain_scale_cm
peak_TE_gain_cm-1
peak_TE_wavelength_nm
peak_TM_gain_cm-1
peak_TM_wavelength_nm
notes
```

将来的には、このログを CSV または JSON として保存し、校正セットを再利用できるようにする。

## 7. 基本方針

校正の順序は次を守る。

```text
遷移波長 -> 井戸幅依存性 -> ひずみ依存性 / TE-TM -> スペクトル幅 -> 絶対ゲイン
```

`gain_scale_cm` は最後に触る。
最初から `gain_scale_cm` で合わせると、バンドギャップ、band offset、ひずみ、線幅の誤差をすべて吸収してしまい、他条件への外挿性が落ちる。
