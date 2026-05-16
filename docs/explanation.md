# BasicMQWDesign.py / MQWGainDesign.py で実施している物理計算の詳細

この文書は、`src/BasicMQWDesign.py` が InP 基板上の InGaAsP または AlGaInAs 多重量子井戸
（MQW: multiple quantum well）活性層を一次設計するために、物理的に何を計算しているかを説明する。
また、`src/MQWGainDesign.py` で追加された簡易 k.p 材料ゲイン計算の扱いも説明する。
実行方法やコマンドライン引数の説明は README に分離し、ここでは計算内容、仮定、近似、出力値の物理的意味に集中する。

本プログラムは、O-band、特に 1.31 um 近傍の InP 系 SOA 活性層を想定した簡易設計ツールである。
厳密な 6x6 / 8x8 k.p 計算、自己無撞着 Poisson-Schrödinger 計算、利得スペクトル計算そのものを置き換えるものではない。
目的は、次のような量を短時間で見積もり、Lumerical MQW ソルバなどへ渡す初期構造を作ることである。

- InP 基板上で所望のひずみを持つ井戸層・障壁層の合金組成
- 擬似格子整合成長した (001) 層の面内ひずみ、面直ひずみ、静水圧ひずみ
- ひずみによる伝導帯端・価電子帯端のシフト
- heavy hole (HH) と light hole (LH) のひずみ分裂
- 有限量子井戸における電子、HH、LH の一次元閉じ込め準位
- e1-hh1 および e1-lh1 遷移エネルギーと遷移波長
- MQW 全体の平均ひずみ、ひずみ厚み積、簡易臨界膜厚
- Lumerical 用 MQW 入力スクリプトに必要な層厚、材料、ひずみ情報
- 簡易 k.p サブバンド分散と TE/TM 材料ゲインスペクトル

## 0. k.p ゲイン計算の位置づけ

`MQWGainDesign.py` は、一次設計で得た井戸層・障壁層・ひずみ・バンド端を z 方向グリッドへ展開し、材料ゲインの初期スクリーニングを行う。
初期実装では、伝導帯はスカラー有効質量 Hamiltonian、価電子帯は軸近似 4x4 Luttinger-Kohn Hamiltonian の 1 つの HH/LH Kramers ブロックを 2x2 有限差分行列として扱う。

この実装で含めているものは次である。

- z 方向の有限差分閉じ込め
- 面内波数 `k_t` sweep
- HH/LH の面内 k 依存混合
- ひずみによる HH/LH バンド端分裂
- キャリア密度からの電子・正孔擬フェルミ準位の数値解
- Lorentzian または Gaussian 線幅による TE/TM ゲインスペクトル

一方、次はまだ簡略化している。

- 6x6 / 8x8 k.p の split-off band / conduction-valence coupling
- 自己無撞着 Poisson-Schrödinger
- nonparabolicity
- 絶対運動量行列要素の材料依存校正
- 実測 PL / gain / absorption に基づく絶対ゲイン校正

したがって、`MQWGainDesign.py` の絶対ゲイン値は `--gain-scale-cm` で校正する前提のスクリーニング値である。
ピーク波長、TE/TM の相対傾向、井戸幅・ひずみ・キャリア密度 sweep の比較を主目的に使う。

## 1. 全体の計算フロー

プログラムは大きく次の順序で計算する。

1. 井戸層と障壁層の材料系を選ぶ。
   対応している材料系は `InGaAsP` と `AlGaInAs` である。
2. 入力された As 組成、Al 組成、目標ひずみなどから、InP 基板上で指定ひずみになる合金組成を決める。
3. 二元化合物の材料パラメータから、多元混晶の格子定数、バンドギャップ、有効質量、Luttinger パラメータ、変形ポテンシャル、弾性定数などを補間する。
4. InP 基板上に擬似格子整合していると仮定し、面内ひずみ `eps_parallel` と面直ひずみ `eps_zz` を計算する。
5. 変形ポテンシャル理論により、ひずみによる伝導帯端と HH/LH 価電子帯端のエネルギーシフトを計算する。
6. 井戸層と障壁層のひずみ込みバンドギャップ差から、伝導帯井戸深さと価電子帯井戸深さを近似的に割り振る。
7. 一次元有効質量 Schrödinger 方程式を有限差分で解き、電子、HH、LH の束縛準位を求める。
8. ひずみ込み井戸バンドギャップに閉じ込めエネルギーを加えて、e1-hh1 と e1-lh1 の遷移エネルギーを見積もる。
9. 遷移エネルギーを波長へ変換し、MQW スタック全体のひずみバランスと簡易臨界膜厚も併せて出力する。

## 2. 材料パラメータ

プログラム内では、まず二元化合物の代表的な 300 K 近傍の材料パラメータを `BIN` に保持している。
登録されている二元材料は次の 5 種類である。

- InP
- InAs
- GaAs
- GaP
- AlAs

各材料は `Material` データ構造として扱われ、主に次の物理量を持つ。

- `a_A`: 格子定数 [Angstrom]
- `Eg_eV`: 無ひずみ Γ 点バンドギャップ [eV]
- `me`: 電子有効質量、自由電子質量 `m0` で規格化
- `gamma1`, `gamma2`, `gamma3`: Luttinger パラメータ
- `ac_eV`: 伝導帯の静水圧変形ポテンシャル [eV]
- `av_eV`: 価電子帯の静水圧変形ポテンシャル [eV]
- `b_eV`: 価電子帯の一軸変形ポテンシャル [eV]
- `C11_GPa`, `C12_GPa`: 弾性定数 [GPa]
- `eps_static`: 静的誘電率

このデータベースは、研究室ごとの校正値ではなく、一般的な III-V 材料パラメータ表に基づくスターター値である。
したがって、PL、XRD、Hall 測定、既知の成長実績に基づく補正を行う前の一次見積もりとして扱う必要がある。

## 3. 混晶組成の表記

### 3.1 InGaAsP

InGaAsP は次の組成で表す。

```text
In_{1-x}Ga_x As_y P_{1-y}
```

プログラム内の変数名では、`x_Ga = x`、`y_As = y` である。
すなわち `x_Ga` は III 族サイトにおける Ga 分率、`y_As` は V 族サイトにおける As 分率である。

InGaAsP では、ユーザーが `y_As` を指定し、`x_Ga` を指定しない場合、目標ひずみから必要な `x_Ga` を逆算する。
このとき、格子定数は Vegard 則で線形補間される。

```text
a = (1 - x)y a_InAs + xy a_GaAs + (1 - x)(1 - y) a_InP + x(1 - y) a_GaP
```

目標ひずみを `strain_target` とすると、目標格子定数 `a_target` は次で与えられる。

```text
a_target = a_InP / (1 + strain_target)
```

プログラムのひずみ定義は後述の通り、

```text
eps_parallel = (a_substrate - a_layer) / a_layer
```

であるため、`strain_target` に負の値を入れると、InP より格子定数が大きい圧縮ひずみ井戸が得られる。
逆に正の値を入れると、InP より格子定数が小さい引張ひずみ層が得られる。

上式の `a = a_target` を `x` について解くことで、指定した As 分率と目標ひずみに対応する Ga 分率を求めている。

### 3.2 AlGaInAs

AlGaInAs は次の組成で表す。

```text
Al_x Ga_y In_{1-x-y} As
```

プログラム内では `x_Al = x`、`y_Ga = y`、`z_In = 1 - x - y` である。
AlGaInAs は V 族が As に固定された III 族混晶として扱う。

`x_Al` が指定され、`y_Ga` が指定されない場合、InP 基板上で目標ひずみを満たすように `y_Ga` を求める。
格子定数は次の線形補間である。

```text
a = x a_AlAs + y a_GaAs + (1 - x - y) a_InAs
```

この式を `y` について解くことで、指定した Al 分率と目標ひずみに対応する Ga 分率を得る。
得られた In 分率 `z_In` が負になる場合は、物理的に不正な組成としてエラーにする。

## 4. 混晶物性の補間

### 4.1 線形補間される物性

格子定数、電子有効質量、Luttinger パラメータ、変形ポテンシャル、弾性定数、誘電率などは、基本的に二元材料の重み付き平均で補間している。

これは Vegard 則または単純な組成線形補間であり、局所的な無秩序、非線形性、温度依存性、成長条件依存性は含まない。
一次設計では、井戸深さやひずみ分裂のおおまかな見積もりには有用だが、精密な利得ピークや偏波利得の予測には不十分である。

### 4.2 InGaAsP のバンドギャップ

InGaAsP の Γ 点バンドギャップは、単純な二元材料の線形補間では実験的な格子整合 InGaAsP の値から大きくずれやすい。
そのためプログラムでは、InP 格子整合 InGaAsP に対してよく用いられる経験式を基準にしている。

```text
Eg_lm(y) = 1.35 - 0.72 y + 0.12 y^2  [eV]
```

ここで `y` は As 分率である。
この式は、InP に格子整合する InGaAsP のバンドギャップを As 分率の関数として近似するものである。

一方、本プログラムでは圧縮ひずみ井戸や引張ひずみ障壁も扱うため、格子整合条件から外れた `x_Ga` も使う。
その場合は、同じ `y_As` で InP 格子整合となる Ga 分率 `x_Ga_lm` を求め、次の補正を加える。

```text
Eg = Eg_lm(y) + [Eg_linear(y, x_Ga) - Eg_linear(y, x_Ga_lm)]
```

ここで `Eg_linear` は二元化合物のバンドギャップを重み付き平均した線形補間値である。
つまり、「As 分率で決まる格子整合 InGaAsP の実験的な傾向」を主成分として使い、そこから Ga 分率をずらして格子不整合を作った分だけ線形補正している。

この扱いは、完全な bowing モデルではない。
しかし、単純線形補間だけで InGaAsP のバンドギャップを決めるより、InP 格子整合系の既知の波長スケールと整合しやすい。

### 4.3 AlGaInAs のバンドギャップ

AlGaInAs では、まず AlAs、GaAs、InAs の Γ 点バンドギャップを線形補間し、その後に pairwise bowing を引いている。

```text
Eg_bowed = Eg_linear
           - b_InGaAs z_In y_Ga
           - b_InAlAs z_In x_Al
           - b_AlGaAs x_Al y_Ga
```

プログラムで使っている bowing 係数は次である。

```text
b_InGaAs = 0.477 eV
b_InAlAs = 0.70  eV
b_AlGaAs = 0.127 eV
```

これは Γ 点直接遷移を意識した簡易モデルであり、X 谷や L 谷、間接遷移、非放物性、温度依存性は明示的には扱わない。

## 5. ひずみの定義と符号規約

本プログラムでは、InP 基板上に擬似格子整合している薄膜を仮定する。
すなわち、成長後の層の面内格子定数は基板に強制的に合わせられる。

面内ひずみは次で定義される。

```text
eps_parallel = (a_substrate - a_layer) / a_layer
```

ここで `a_substrate` は InP の格子定数、`a_layer` は自由状態での層の格子定数である。

この符号規約では、層の自然格子定数が InP より大きい場合、

```text
a_layer > a_substrate
eps_parallel < 0
```

となる。
この層は基板に合わせるため面内方向に縮められるので、圧縮ひずみである。
したがって、本プログラムでは負の `eps_parallel` が圧縮ひずみを意味する。

逆に、層の自然格子定数が InP より小さい場合、

```text
a_layer < a_substrate
eps_parallel > 0
```

となる。
この層は面内方向に引き伸ばされるため、引張ひずみである。

この符号は Lumerical MQW のひずみ符号規約に合わせている。

## 6. (001) 擬似格子整合層の面直ひずみ

InP (001) 基板上の薄膜を想定し、面内方向には等方的に同じひずみがかかると仮定する。

```text
eps_xx = eps_yy = eps_parallel
```

薄膜表面の面直方向には外力がない近似を置くと、立方晶の弾性関係から面直ひずみは次で与えられる。

```text
eps_zz = -2 (C12 / C11) eps_parallel
```

圧縮ひずみ層では `eps_parallel < 0` なので、一般に `eps_zz > 0` となる。
これは、面内で縮められた層が面直方向に伸びることを表している。

引張ひずみ層では `eps_parallel > 0` なので、一般に `eps_zz < 0` となる。
これは、面内で引き伸ばされた層が面直方向に縮むことを表している。

静水圧ひずみ成分はトレースであり、次のように計算する。

```text
hydrostatic_strain = eps_xx + eps_yy + eps_zz
                   = 2 eps_parallel + eps_zz
```

この静水圧成分は主にバンドギャップ全体の増減に効き、後述の一軸成分は HH/LH 分裂に効く。

## 7. ひずみによるバンド端シフト

ひずみが入ると、伝導帯端と価電子帯端が変形ポテンシャルを通じてシフトする。
本プログラムでは、(001) 二軸ひずみの一次近似として次を計算している。

伝導帯端のシフトは、静水圧変形ポテンシャル `ac` を用いて、

```text
dEc = ac * (2 eps_parallel + eps_zz)
```

とする。

価電子帯の静水圧シフトは、

```text
dEv_hydro = av * (2 eps_parallel + eps_zz)
```

である。

価電子帯ではさらに、二軸ひずみによって HH と LH が分裂する。
プログラムでは一軸変形ポテンシャル `b` を用いて、HH と LH の価電子帯端シフトを次のように置いている。

```text
dEv_hh = dEv_hydro - b (eps_zz - eps_parallel)
dEv_lh = dEv_hydro + b (eps_zz - eps_parallel)
```

この結果、ひずみ込みの HH 遷移に対応するバンドギャップと LH 遷移に対応するバンドギャップは、

```text
Eg_hh = Eg + dEc - dEv_hh
Eg_lh = Eg + dEc - dEv_lh
```

として計算される。

`Eg_hh` は、伝導帯端から HH 価電子帯端までのエネルギー差である。
`Eg_lh` は、伝導帯端から LH 価電子帯端までのエネルギー差である。

HH/LH の分裂量は、

```text
hh_lh_split = dEv_hh - dEv_lh
```

として出力される。

圧縮ひずみ井戸では、一般に HH バンドが LH バンドより上に来やすく、TE 偏波利得に有利な設計になる。
引張ひずみでは逆に LH 的な成分が相対的に上がり、TM 偏波や偏波無依存設計に関係する。
ただし、このプログラムは HH/LH の厳密な混合を解いていないため、偏波利得そのものは計算していない。

## 8. 量子井戸のバンドオフセット近似

井戸層と障壁層のひずみ込みバンドギャップを計算した後、プログラムは障壁と井戸のバンドギャップ差を量子井戸深さとして扱う。

HH 遷移については、

```text
dEg_hh = max(0, Eg_hh_barrier - Eg_hh_well)
```

LH 遷移については、

```text
dEg_lh = max(0, Eg_lh_barrier - Eg_lh_well)
```

を計算する。
`max(0, ...)` としているのは、障壁の方が井戸より狭ギャップになった場合に、井戸として束縛できない状況を負の井戸深さとして扱わないためである。

伝導帯オフセット比 `qc` を使って、HH 遷移に対するギャップ差を伝導帯側と価電子帯側に分配する。

```text
dEc    = qc       * dEg_hh
dEv_hh = (1 - qc) * dEg_hh
```

LH については、価電子帯側の井戸深さを次のように置く。

```text
dEv_lh = (1 - qc) * dEg_lh
```

このモデルは非常に単純なバンドオフセット近似である。
実際には、伝導帯端と価電子帯端の絶対位置、model-solid theory、電子親和力、界面双極子、組成依存 valence band offset などを考慮する方が望ましい。
本プログラムでは、初期設計のために「バンドギャップ差の一定割合が伝導帯に行く」という近似を採用している。

## 9. 有限量子井戸の一次元有効質量計算

量子閉じ込め準位は、成長方向 `z` の一次元有限井戸問題として解く。
井戸幅は `well_nm`、左右の障壁厚はそれぞれ `barrier_nm` とし、計算領域全体は次である。

```text
total_nm = well_nm + 2 barrier_nm
```

井戸中心を `z = 0` とし、

```text
|z| <= well_nm / 2
```

の領域を井戸、外側を障壁として扱う。

ポテンシャルは井戸内を 0 eV、障壁内を井戸深さ `V0` とする矩形有限井戸である。

```text
V(z) = 0   in well
V(z) = V0  in barrier
```

ここで `V0` は電子なら `dEc`、HH なら `dEv_hh`、LH なら `dEv_lh` である。

解いている方程式は、有効質量近似の一次元 Schrödinger 方程式である。

```text
[- d/dz { (hbar^2 / 2 m*(z)) d/dz } + V(z)] psi(z) = E psi(z)
```

有効質量 `m*(z)` は、井戸内では井戸材料の有効質量、障壁内では障壁材料の有効質量に切り替わる。
この空間依存有効質量を有限差分で離散化し、疎行列固有値問題として最小固有値から求めている。

離散化に使っている定数は、

```text
hbar^2 / (2 m0) = 0.0380998212 eV nm^2
```

である。
したがって、長さを nm、エネルギーを eV、有効質量を `m0` 規格化で扱える。

境界条件は、計算領域の両端で波動関数が 0 になる Dirichlet 境界条件である。
左右に有限厚の障壁を置いているため、井戸中の束縛状態は障壁側へ指数関数的に染み出し、その外側の領域端で 0 になる近似で計算される。
障壁厚が十分大きければ、領域端の境界条件が井戸準位へ与える影響は小さい。

求めた固有値のうち、障壁高さより低いものだけを束縛準位として採用する。

```text
E < 0.999 V0
```

電子では最も低い準位が `e1`、HH では `hh1`、LH では `lh1` として使われる。

## 10. HH と LH の有効質量

電子の成長方向有効質量は、材料パラメータ `me` をそのまま使う。

正孔については、Luttinger パラメータから (001) 成長方向の有効質量を近似している。

HH の成長方向有効質量は、

```text
mhh_z = 1 / (gamma1 - 2 gamma2)
```

LH の成長方向有効質量は、

```text
mlh_z = 1 / (gamma1 + 2 gamma2)
```

である。

この式は、バンド混合を無視した単純な軸方向有効質量である。
実際の量子井戸では HH/LH の混合、ひずみによる分裂、面内波数依存、非放物性が効くため、偏波利得や高キャリア密度でのスペクトル形状には k.p 計算が必要になる。
ただし、基底準位の閉じ込めエネルギーとおおまかな遷移波長を見積もる目的には有用である。

## 11. 遷移エネルギーと波長

井戸層のひずみ込みバンドギャップに、電子と正孔の閉じ込めエネルギーを加えて遷移エネルギーを見積もる。

HH 遷移は、

```text
E_transition_hh = Eg_hh_well + Ee1 + Ehh1
```

LH 遷移は、

```text
E_transition_lh = Eg_lh_well + Ee1 + Elh1
```

である。

ここで `Ee1` は電子の基底準位閉じ込めエネルギー、`Ehh1` は HH の基底準位閉じ込めエネルギー、`Elh1` は LH の基底準位閉じ込めエネルギーである。
閉じ込めエネルギーは井戸端から測った正の値なので、量子井戸を薄くすると一般に遷移エネルギーは高くなり、波長は短くなる。

エネルギーから波長への変換には、

```text
hc = 1.239841984 eV um
lambda_um = hc / E_transition_eV
```

を使う。

この遷移波長は、吸収端または利得ピークの厳密値ではない。
実際の発光・利得ピークは、キャリア注入によるバンドフィリング、バンドギャップ renormalization、励起子効果、温度、内部電場、井戸幅揺らぎ、界面粗さ、散乱幅などで変わる。
本プログラムの値は、MQW 構造の初期設計における「量子閉じ込め込みのバンド端遷移波長」として読むべきである。

## 12. MQW スタックとひずみバランス

プログラムは、井戸層と障壁層から次のような対称スタックを構成する。

```text
barrier / (well / barrier) repeated N times
```

井戸数を `N` とすると、層数は `2N + 1` である。
先頭と末尾に障壁があり、その間に井戸と障壁が交互に入る。

各層のひずみ `eps_i` と厚み `t_i` から、全厚は、

```text
total_nm = sum_i t_i
```

平均ひずみは、

```text
average_strain = sum_i (eps_i t_i) / sum_i t_i
```

である。

また、符号付きひずみ厚み積は、

```text
signed_strain_thickness_nm = sum_i (eps_i t_i)
```

絶対値ひずみ厚み積は、

```text
absolute_strain_thickness_nm = sum_i (|eps_i| t_i)
```

として出力される。

圧縮ひずみ井戸と引張ひずみ障壁を組み合わせると、平均ひずみを 0 付近へ近づけることができる。
これは strain compensation と呼ばれる考え方で、MQW 全体としての緩和や転位発生を抑えながら、個々の井戸には大きめの圧縮ひずみを入れて TE 利得を高める設計に使われる。

ただし、平均ひずみが 0 であっても、個々の層のひずみ、井戸数、総膜厚、成長温度、界面品質によって緩和条件は変わる。
平均ひずみだけで結晶品質が保証されるわけではない。

## 13. Matthews-Blakeslee 臨界膜厚の簡易見積もり

プログラムは、単層の格子不整合に対して Matthews-Blakeslee 型の簡易臨界膜厚を見積もる。

入力として使うミスマッチは、

```text
f = |eps_parallel|
```

である。

InP 基板の格子定数から Burgers ベクトルを近似的に、

```text
b = a_substrate / sqrt(2)
```

としている。
単位は nm に直して使う。

60 度転位を仮定し、Poisson 比 `nu = 0.35`、転位角 `alpha = 60 deg` のもとで、次のような係数を作る。

```text
pref = b / (4 pi f) * (1 - nu cos^2 alpha) / (1 + nu)
```

そのうえで、対数項を含む臨界膜厚方程式を反復で解く。
実装では数値的に極端な値になりすぎないよう、初期値や下限を置いたスクリーニング用の式になっている。

この臨界膜厚は、あくまで単層の目安である。
MQW では、圧縮層と引張層の積層、界面の転位核生成、成長中の緩和、熱履歴などが絡むため、実際の許容膜厚はこの値だけでは決まらない。
それでも、井戸層や障壁層単体のミスマッチが過大でないかを初期段階で見る指標として役に立つ。

## 14. 既定設計値の物理的意図

### 14.1 AlGaInAs 既定設計

AlGaInAs 系では、既定で次のような考え方を置いている。

- 井戸層は Al 分率を比較的低くし、狭ギャップにする。
- 井戸層には圧縮ひずみを入れる。
- 障壁層は Al 分率を高くし、井戸より広ギャップにする。
- 障壁層には引張ひずみを入れ、井戸の圧縮ひずみを補償する。
- 伝導帯オフセット比 `qc` は 0.65 とする。

AlGaInAs は伝導帯オフセットを比較的大きく取りやすく、電子閉じ込めに有利な系として扱っている。
圧縮ひずみ井戸により HH 遷移を基底側に置き、O-band 近傍の TE 利得を狙う一次設計である。

### 14.2 InGaAsP 既定設計

InGaAsP 系では、既定で次のような考え方を置いている。

- 井戸層は As 分率を高め、狭ギャップにする。
- 井戸層には圧縮ひずみを入れる。
- 障壁層は As 分率を低めにし、広ギャップにする。
- 障壁層には引張ひずみを入れ、井戸の圧縮ひずみを補償する。
- 伝導帯オフセット比 `qc` は 0.40 とする。

InGaAsP では、量子閉じ込めによって遷移エネルギーが高くなり、波長が短くなる。
そのため、井戸層の無ひずみ・無閉じ込めバンドギャップは、最終的に狙う遷移波長より長波長側になるように選ぶ必要がある。
プログラムの既定値は、閉じ込めエネルギーとひずみシフトを含めた e1-hh1 遷移が O-band 近傍になるように置かれている。

## 15. Lumerical 用出力に含まれる物理量

プログラムは、計算した設計候補から Lumerical script command 用の MQW 入力断片も生成する。
この出力は、Lumerical 側の MQW ゲイン計算へ渡すためのスターター構造である。

含まれる主な物理量は次である。

- 井戸材料と障壁材料の組成
- 各層の厚み
- 各層の面内ひずみ
- 井戸数に対応した barrier / well / barrier ... の積層順序
- Lorentzian linewidth の初期値
- キャリア密度掃引の初期値
- 価電子帯モデル指定の初期値

ここで重要なのは、Python 側で計算している遷移波長は簡易有効質量モデルによる見積もりであり、Lumerical 側の k.p MQW 計算とはモデルの詳細が異なるという点である。
Python 側の結果は、Lumerical に渡す初期構造を選ぶためのスクリーニング値として扱う。

## 16. このプログラムが含まない物理

本プログラムは一次設計用なので、次の物理は明示的には含まない。

- 6x6 / 8x8 k.p による価電子帯混合
- 面内波数 `k_parallel` 依存のサブバンド分散
- 非放物性
- スピン軌道分裂帯との混合
- 自己無撞着 Poisson-Schrödinger 計算
- キャリア注入によるバンドフィリング
- バンドギャップ renormalization
- 励起子効果
- 多体効果
- 光学利得、吸収、自然放出スペクトルの直接計算
- TE/TM 行列要素の厳密計算
- 界面粗さ、井戸幅揺らぎ、組成揺らぎ
- 温度依存性
- X 谷、L 谷、Γ-X クロスオーバーの詳細評価
- 絶対バンドアラインメントに基づく厳密なバンドオフセット
- piezoelectric field や built-in field
- ドーピング、自由キャリア遮蔽、内部損失

特に、偏波利得や発振/増幅スペクトルを議論するには、Python の簡易値だけで判断せず、k.p ベースの MQW ソルバで再計算する必要がある。

## 17. 出力値を読むときの注意

`well_strain` と `barrier_strain` は、各層単体が InP 基板に擬似格子整合したときの面内ひずみである。
負なら圧縮、正なら引張である。

`dEc_eV`、`dEv_hh_eV`、`dEv_lh_eV` は、単純なバンドギャップ差分配モデルによる井戸深さである。
実際の conduction band offset / valence band offset と完全に一致するとは限らない。

`e1_eV`、`hh1_eV`、`lh1_eV` は、それぞれ井戸端から測った閉じ込めエネルギーである。
バンドギャップそのものではない。

`E_transition_hh_eV` と `lambda_transition_hh_um` は、ひずみ込み HH バンドギャップに e1 と hh1 の閉じ込めエネルギーを加えた遷移である。
通常、圧縮ひずみ井戸ではこの e1-hh1 遷移が TE 利得設計で重要になる。

`E_transition_lh_eV` と `lambda_transition_lh_um` は、同様に LH に対する見積もりである。
HH/LH の分離が十分大きいか、LH 遷移がどの程度短波長側にあるかを見る目安になる。

`average_strain` は MQW 全厚で重み付けした平均ひずみである。
0 に近いほど strain-balanced に近いが、結晶緩和を完全に予測する量ではない。

`critical_thickness_*_nm_est` は Matthews-Blakeslee 型の単層臨界膜厚の簡易見積もりである。
MQW 全体の成長可能性を保証するものではなく、過大ひずみを検出するためのスクリーニング指標である。

## 18. 実設計へ進める際の推奨確認

このプログラムで候補構造を得た後、実際のエピ設計やデバイス設計へ進める場合は、少なくとも次を確認するのが望ましい。

- 使用する MOCVD/MBE 条件に対して校正済みの格子定数・バンドギャップ・組成関係に置き換える。
- XRD で得られる組成・ひずみと、プログラムの組成・ひずみの対応を確認する。
- PL または吸収測定により、井戸幅と組成に対する遷移波長の補正量を決める。
- Lumerical MQW、nextnano などで k.p サブバンドと TE/TM 行列要素を計算する。
- キャリア密度依存の利得ピーク、透明キャリア密度、偏波利得差を評価する。
- 障壁高さが電子・正孔の熱漏れに十分かを、動作温度で確認する。
- 平均ひずみだけでなく、総 MQW 厚、個別層臨界膜厚、成長実績を含めて緩和リスクを判断する。
- SOA の場合は、活性層だけでなく SCH 層、クラッド層、導波路閉じ込め係数、内部損失も含めて最終設計する。

## 19. まとめ

本プログラムは、InP 系 MQW 活性層に対して、組成、ひずみ、ひずみ込みバンドギャップ、有限井戸閉じ込め準位、遷移波長、ひずみバランスを一貫して見積もる。

物理モデルとしては、次の近似を組み合わせている。

- 混晶物性の Vegard 則・簡易 bowing 補間
- InP 基板上の (001) 擬似格子整合二軸ひずみ
- 変形ポテンシャルによる伝導帯・価電子帯シフト
- HH/LH の単純なひずみ分裂
- 一次元有効質量有限井戸
- バンドギャップ差を `qc` で分配する簡易バンドオフセット
- 厚み重み付き平均によるひずみバランス評価
- Matthews-Blakeslee 型の単層臨界膜厚スクリーニング

したがって、出力される遷移波長は「ひずみと量子閉じ込めを含めた一次近似のバンド端遷移」であり、最終的な利得ピークや実測 PL 波長そのものではない。
この値を起点に、校正済み材料パラメータと k.p MQW 計算へ進める、という位置づけで使うのが適切である。

## 20. 使用変数と説明の対応

### 20.1 物理定数

| name | explanation |
| --- | --- |
| `C_LIGHT` | 真空中の光速 [m/s]。 |
| `M0` | 自由電子質量 [kg]。一部の物理定数確認用として保持している。 |
| `HC_EV_UM` | `h c` を `eV um` 単位で表した値。遷移エネルギー [eV] から波長 [um] へ変換するために使う。 |
| `H_PRANCK` | Planck 定数 [eV s]。 |
| `HBAR` | 換算 Planck 定数 [eV s]。 |
| `HBAR2_OVER_2M0_EV_NM2` | `hbar^2 / (2 m0)` を `eV nm^2` 単位で表した値。一次元有効質量 Schrödinger 方程式の有限差分 Hamiltonian に使う。 |

### 20.2 材料パラメータ

| name | explanation |
| --- | --- |
| `Material` | 材料組成と物性値をまとめるデータ構造。二元材料、InGaAsP、AlGaInAs のいずれもこの形で扱う。 |
| `name` | 材料名または組成式を表す文字列。例: `In0.8245Ga0.1755As0.5670P0.4330`。 |
| `family` | 材料ファミリ。`binary`、`InGaAsP`、`AlGaInAs` など。 |
| `x_Ga` | InGaAsP の Ga 分率。組成 `In_{1-x}Ga_x As_y P_{1-y}` の `x`。 |
| `y_As` | InGaAsP の As 分率。組成 `In_{1-x}Ga_x As_y P_{1-y}` の `y`。 |
| `x_Al` | AlGaInAs の Al 分率。組成 `Al_x Ga_y In_{1-x-y} As` の `x`。 |
| `y_Ga` | AlGaInAs の Ga 分率。組成 `Al_x Ga_y In_{1-x-y} As` の `y`。 |
| `z_In` | AlGaInAs の In 分率。`z_In = 1 - x_Al - y_Ga`。 |
| `a_A` | 自由状態の格子定数 [Angstrom]。InP 基板との格子不整合とひずみ計算に使う。 |
| `Eg_eV` | 無ひずみ Γ 点バンドギャップ [eV]。ひずみ込みバンドギャップ計算の基準値。 |
| `me` | 電子有効質量。自由電子質量 `m0` で規格化した値。 |
| `gamma1` | Luttinger パラメータ。正孔の有効質量計算に使う。 |
| `gamma2` | Luttinger パラメータ。HH/LH の成長方向有効質量計算に使う。 |
| `gamma3` | Luttinger パラメータ。現行の一次元計算では保持しているが、Hamiltonian には直接入れていない。 |
| `ac_eV` | 伝導帯の静水圧変形ポテンシャル [eV]。ひずみによる伝導帯端シフト `dEc` に使う。 |
| `av_eV` | 価電子帯の静水圧変形ポテンシャル [eV]。価電子帯の静水圧シフトに使う。 |
| `b_eV` | 価電子帯の一軸変形ポテンシャル [eV]。HH/LH ひずみ分裂に使う。 |
| `C11_GPa` | 弾性定数 `C11` [GPa]。面直ひずみ `eps_zz` の計算に使う。 |
| `C12_GPa` | 弾性定数 `C12` [GPa]。面直ひずみ `eps_zz` の計算に使う。 |
| `eps_static` | 静的誘電率。ひずみの epsilon ではない。現行の閉じ込め準位計算には直接使っていない。 |
| `source_note` | 材料パラメータの由来や補間方法を記録するメモ。 |
| `mhh_z` | (001) 成長方向の HH 有効質量。`1 / (gamma1 - 2 gamma2)`。 |
| `mlh_z` | (001) 成長方向の LH 有効質量。`1 / (gamma1 + 2 gamma2)`。 |
| `BIN` | InP、InAs、GaAs、GaP、AlAs の二元材料パラメータ辞書。 |
| `PARAM_KEYS` | 混晶補間する材料パラメータ名のリスト。 |

### 20.3 入力引数・設計条件

| name | explanation |
| --- | --- |
| `family` / `--family` | 材料系。`ingaasp` または `algainas`。既定は `ingaasp`。 |
| `wells` / `--wells` | 井戸数。MQW スタックは `barrier / (well / barrier) x wells` として構成される。 |
| `well_nm` / `--well-nm` | 井戸層厚 [nm]。有限井戸幅として使う。 |
| `barrier_nm` / `--barrier-nm` | 障壁層厚 [nm]。有限井戸計算領域の左右障壁厚、および MQW スタックの障壁厚として使う。 |
| `q_c` / `--qc` | 伝導帯オフセット比。バンドギャップ差のうち伝導帯側へ割り当てる割合。 |
| `well_strain` / `--well-strain` | 井戸層の目標面内ひずみ。負は圧縮、正は引張。 |
| `barrier_strain` / `--barrier-strain` | 障壁層の目標面内ひずみ。負は圧縮、正は引張。 |
| `strain_target` | 組成を逆算するときの目標面内ひずみ。`eps_static` と混同しないよう、ひずみ側はこの名前にしている。 |
| `al_well` / `--al-well` | AlGaInAs 井戸層の Al 分率 `x_Al`。 |
| `al_barrier` / `--al-barrier` | AlGaInAs 障壁層の Al 分率 `x_Al`。 |
| `as_well` / `--as-well` | InGaAsP 井戸層の As 分率 `y_As`。 |
| `as_barrier` / `--as-barrier` | InGaAsP 障壁層の As 分率 `y_As`。 |
| `x_Ga` | InGaAsP の Ga 分率を直接指定する場合の値。未指定なら `strain_target` から逆算する。 |
| `y_Ga` | AlGaInAs の Ga 分率を直接指定する場合の値。未指定なら `strain_target` から逆算する。 |
| `outer_barrier_nm` | MQW スタック最外側障壁の厚み [nm]。未指定なら `barrier_nm` と同じ。 |
| `dz_nm` | 有限差分メッシュ間隔 [nm]。小さいほど精度は上がるが計算量が増える。 |
| `n_eigs` | 有限井戸ソルバで求める固有値数の上限。 |
| `json` / `--json` | JSON 出力先。既定は `out/ingaasp_design_result.json`。 |
| `lsf` / `--lsf` | Lumerical script 出力先。既定は `out/ingaasp_lumerical_input.lsf`。 |

### 20.4 ひずみ・バンド端計算

| name | explanation |
| --- | --- |
| `substrate` | 基板材料。既定では InP。 |
| `a_target` | 目標ひずみを満たすための自由状態格子定数 [Angstrom]。`a_InP / (1 + strain_target)`。 |
| `a_no_ga` | InGaAsP で Ga 分率を 0 としたときの格子定数。`x_Ga` 逆算に使う。 |
| `a_all_ga` | InGaAsP で Ga 分率を 1 としたときの格子定数。`x_Ga` 逆算に使う。 |
| `x_Ga_lm` | 同じ `y_As` で InP 格子整合になる InGaAsP の Ga 分率。 |
| `eg_lm` | InP 格子整合 InGaAsP の経験式バンドギャップ [eV]。 |
| `eg_linear_delta` | 格子整合組成から `x_Ga` をずらしたことによる線形補間バンドギャップ補正 [eV]。 |
| `eg_bowed` | AlGaInAs の pairwise bowing 補正込み Γ 点バンドギャップ [eV]。 |
| `eps_parallel` | 面内ひずみ。`(a_substrate - a_layer) / a_layer`。負は圧縮、正は引張。 |
| `eps_zz` | 面直ひずみ。`-2 (C12 / C11) eps_parallel`。 |
| `hydrostatic_strain` / `hydro` | 静水圧ひずみ成分。`2 eps_parallel + eps_zz`。 |
| `dEc` / `dEc_eV` | ひずみによる伝導帯端シフト、または遷移計算では伝導帯側井戸深さ [eV]。文脈により意味が変わるため注意。 |
| `dEv_hydro` | ひずみによる価電子帯の静水圧シフト [eV]。 |
| `dEv_hh` / `dEv_hh_eV` | HH 価電子帯端シフト、または遷移計算では HH 価電子帯側井戸深さ [eV]。 |
| `dEv_lh` / `dEv_lh_eV` | LH 価電子帯端シフト、または遷移計算では LH 価電子帯側井戸深さ [eV]。 |
| `Eg_hh_eV` | ひずみ込みの e-HH バンドギャップ [eV]。 |
| `Eg_lh_eV` | ひずみ込みの e-LH バンドギャップ [eV]。 |
| `hh_lh_split_eV` | HH と LH の価電子帯端分裂量 [eV]。`dEv_hh - dEv_lh`。 |
| `well_strain` | 井戸層のひずみ計算結果一式。`eps_parallel`、`eps_zz`、`Eg_hh_eV` などを含む辞書。 |
| `barrier_strain` | 障壁層のひずみ計算結果一式。 |

### 20.5 有限井戸・遷移計算

| name | explanation |
| --- | --- |
| `well` | 井戸層材料の `Material`。 |
| `barrier` | 障壁層材料の `Material`。 |
| `total_nm` | 有限井戸計算領域の全幅、または MQW スタック全厚 [nm]。文脈により意味が変わる。 |
| `z` | 有限差分計算における成長方向座標 [nm]。 |
| `interior` | Dirichlet 境界の両端を除いた内部格子点。 |
| `in_well` | 各格子点が井戸内にあるかを表す真偽値配列。 |
| `V` | 有限井戸ポテンシャル [eV]。井戸内 0、障壁内 `barrier_height_eV`。 |
| `m` | 位置依存有効質量。井戸内は井戸材料、障壁内は障壁材料の値。 |
| `invm` | `1 / m`。位置依存有効質量 Hamiltonian の離散化に使う。 |
| `H` | 有限差分で作った一次元有効質量 Hamiltonian 疎行列。 |
| `vals` | Hamiltonian の固有値候補 [eV]。障壁高さより低いものを束縛準位として返す。 |
| `barrier_height_eV` | 有限井戸の障壁高さ [eV]。電子では `dEc_eV`、HH では `dEv_hh_eV`、LH では `dEv_lh_eV`。 |
| `dEg_hh_eV` | 障壁と井戸の HH ひずみ込みバンドギャップ差 [eV]。 |
| `dEg_lh_eV` | 障壁と井戸の LH ひずみ込みバンドギャップ差 [eV]。 |
| `electron_levels_eV` | 電子の束縛準位リスト [eV]。井戸伝導帯端から測った閉じ込めエネルギー。 |
| `hh_levels_eV` | HH の束縛準位リスト [eV]。井戸 HH 価電子帯端から測った閉じ込めエネルギー。 |
| `lh_levels_eV` | LH の束縛準位リスト [eV]。井戸 LH 価電子帯端から測った閉じ込めエネルギー。 |
| `Ee1` / `e1_eV` | 電子基底準位の閉じ込めエネルギー [eV]。 |
| `Ehh1` / `hh1_eV` | HH 基底準位の閉じ込めエネルギー [eV]。 |
| `Elh1` / `lh1_eV` | LH 基底準位の閉じ込めエネルギー [eV]。 |
| `Etr_hh` / `E_transition_hh_eV` | e1-hh1 遷移エネルギー [eV]。`Eg_hh_well + Ee1 + Ehh1`。 |
| `Etr_lh` / `E_transition_lh_eV` | e1-lh1 遷移エネルギー [eV]。`Eg_lh_well + Ee1 + Elh1`。 |
| `lambda_transition_hh_um` | e1-hh1 遷移波長 [um]。`HC_EV_UM / E_transition_hh_eV`。 |
| `lambda_transition_lh_um` | e1-lh1 遷移波長 [um]。`HC_EV_UM / E_transition_lh_eV`。 |

### 20.6 MQW スタック・出力

| name | explanation |
| --- | --- |
| `layers` | MQW スタックを表す `(Material, thickness_nm)` のリスト。 |
| `sum_eps_t` | `sum(eps_i t_i)`。符号付きひずみ厚み積 [nm]。 |
| `sum_abs_eps_t` | `sum(|eps_i| t_i)`。絶対値ひずみ厚み積 [nm]。 |
| `average_strain` | 厚み重み付き平均ひずみ。`sum_eps_t / total_nm`。 |
| `signed_strain_thickness_nm` | 符号付きひずみ厚み積 [nm]。圧縮と引張の相殺を見る。 |
| `absolute_strain_thickness_nm` | 絶対値ひずみ厚み積 [nm]。総ひずみ量の目安。 |
| `abs_mismatch` | Matthews-Blakeslee 臨界膜厚見積もりに使う格子不整合の絶対値。 |
| `critical_thickness_well_nm_est` | 井戸層単体の簡易臨界膜厚見積もり [nm]。 |
| `critical_thickness_barrier_nm_est` | 障壁層単体の簡易臨界膜厚見積もり [nm]。 |
| `design` | 設計結果全体をまとめた辞書。JSON 出力と Lumerical 出力の元データ。 |
| `transition` / `tr` | 遷移計算結果をまとめた辞書。井戸深さ、閉じ込め準位、遷移エネルギー、遷移波長を含む。 |
| `strain_balance` / `sb` | MQW スタックのひずみバランス計算結果をまとめた辞書。 |
| `well_mat` | Lumerical script 内で使う井戸材料構造体名。 |
| `barrier_mat` | Lumerical script 内で使う障壁材料構造体名。 |
| `stack.length` | Lumerical script に出力する各層厚 [m]。 |
| `stack.material` | Lumerical script に出力する各層の材料名配列。 |
| `stack.strain` | Lumerical script に出力する各層の面内ひずみ配列。 |
| `stack.gamma` | Lumerical script の Lorentzian linewidth 初期値 [eV]。 |
| `sim.T` | Lumerical script の温度 [K]。 |
| `sim.kt` | Lumerical script の面内波数サンプリング範囲。 |
| `sim.cden` | Lumerical script の平均キャリア密度掃引 [m^-3]。 |
