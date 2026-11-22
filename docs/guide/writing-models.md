# モデルのフック

`MaterialModel` に書けるフック (メソッドと宣言) の仕様を示す。モデルを書いて検査する手順は [モデルを書く](model-python.md) にあり、このページは個々のフックを扱う。

コードは、組み込みモデルのソース (`src/returnmap/models/`) からの抜粋である。

| フック | 区分 | 役割 | 節 |
|---|---|---|---|
| `param_names`、`param_meta`、`dimension` | 宣言 | パラメータと次元 | 2 |
| `Explicit` / `Implicit` | 宣言 | 状態変数 | 3、4 |
| `yield_function(state)` | 必須 | 降伏関数 | 3 |
| `update_state(new, dlambda, inc)` | `Explicit` の状態がある場合 | 閉形式の更新 | 3 |
| `state_residual(new, dlambda, inc)` | `Implicit` の状態がある場合 | 残差 | 4 |
| `elastic_stiffness(state)` | 任意 | 弾性剛性 | 3 |
| `flow(state)` | 任意 | 流れ方向 | 5 |
| `elastic_check(inc)`、`elastic_update_state(inc)` | 任意 | 弾塑性の判定、弾性ステップでの状態の更新 | 7 |
| `elastic_predictor(state_n, *, strain_inc)`、`elastic_tangent(inc)` | 任意 | 弾性予測子と、その接線 | 8 |
| `user_defined_return_mapping(inc)`、`user_defined_tangent(new, dlambda, inc)` | 任意 | 手書きの解法と接線 | 9 |

フックは位置引数で呼ばれる。引数リストが規定と異なるモデルは、クラス定義時に `TypeError` になる。

## 1. フックへの入力

予測子以外のフックは、1 ステップの入力をまとめた `Increment` (`inc`) を受け取る。構成は UMAT と同じである (冒頭で試行状態を計算し、そこから分岐と補正を行う)。

| 属性 | 内容 |
|---|---|
| `inc.state_n` | ステップ開始時の状態。`inc.state_n["stress"]` は σ_n |
| `inc.strain_inc` | ひずみ増分 Δε (せん断は工学ひずみ) |
| `inc.trial` | 弾性試行状態。`inc.trial["stress"]` は σ_trial |

`new` はステップ終了時の状態、`dlambda` は塑性乗数の増分 Δλ である。状態は読み取り専用であり、作業用の dict は `dict(inc.state_n)` で作る。

## 2. パラメータと次元

| 項目 | 規定 |
|---|---|
| `param_names` | パラメータの唯一の宣言。基底クラスが同名のキーワード引数を受け取り、同名の属性 (`self.E` など) に設定する。順序は UMAT の `PROPS` の順である |
| 既定値 | 持てない。全パラメータが必須である |
| `__init__` | 書かない。例外は下記の「パラメータ以外の設定」のみ |
| 定数 | クラス属性に書く |
| パラメータから導く値 | `@property` に書く |
| `param_meta` | 各パラメータの説明、単位、目安の値、範囲。任意 |

| 次元の決め方 | 動作 |
|---|---|
| クラス本体に `dimension` を書かない | 構築時の `dimension=` で選ぶ (省略時は 3D) |
| クラス本体に `dimension = X` を書く | そのクラスと子クラスは、その次元に固定される。別の `dimension=` を渡すと `TypeError` になる |

閉形式のリターンマッピングや手導出ブロックのように、代数が特定の次元でのみ成り立つクラスは、次元を固定する。物理は次元に依存しない形で親クラスに書き、次元を固定したクラスは `dimension` のみを持つ。

**パラメータ以外の設定** (数値的な正則化の係数など) が必要な場合に限り、`__init__` を定義し、残りを `super().__init__(**params)` に渡す。

```{literalinclude} ../../src/returnmap/models/j2_damage/unilateral.py
:language: python
:pyobject: J2IsotropicDamageUnilateral.__init__
:dedent: 4
```

## 3. Explicit の状態

閉形式で更新できる状態は `Explicit` で宣言し、`update_state` に更新を書く。`Explicit` の状態は、応力更新の未知数に含まれない。未知数は `[σ, Δλ]` である。

```{literalinclude} ../../src/returnmap/models/j2_isotropic/model.py
:language: python
:start-at: "ep = Explicit(shape=SCALAR"
:end-at: 'return [self.ep(inc.state_n["ep"] + dlambda)]'
:dedent: 4
```

| 項目 | 規定 |
|---|---|
| 宣言 | `名前 = Explicit(shape=SCALAR)` (スカラー)、`Explicit(shape=NTENS)` (テンソル) |
| `update_state` の戻り値 | `self.<名前>(新しい値)` のリスト |
| 応力 | 宣言しない。σ は常に未知数であり、残差は既定の関連流れのものが用いられる |

### 弾性剛性

`elastic_stiffness(state)` の既定は、`self.E` と `self.nu` から作る等方剛性である。剛性が状態によって変わるモデルのみ、定義する。

```{literalinclude} ../../src/returnmap/models/j2_damage/damage.py
:language: python
:pyobject: J2IsotropicDamage.elastic_stiffness
:dedent: 4
```

## 4. Implicit の状態

閉形式で解けない発展則を持つ状態は `Implicit` で宣言し、`state_residual` に残差を書く。`Implicit` の状態は、σ、Δλ と並んで応力更新の未知数になる (`[σ, Δλ, implicit の状態]`)。

次は `OWKinematic` の残差である。状態は `stress = Implicit(shape=NTENS)`、`alpha = Implicit(shape=NTENS)`、`ep = Implicit(shape=SCALAR)` と宣言されている。

```{literalinclude} ../../src/returnmap/models/ow_kinematic/model.py
:language: python
:pyobject: OWKinematic.state_residual
:dedent: 4
```

| 項目 | 規定 |
|---|---|
| `state_residual` の戻り値 | `self.<名前>(残差)` のリスト |
| σ の残差 | `stress = Implicit(shape=NTENS)` を宣言した場合は、σ の行の残差も書く。宣言しない場合は、既定の関連流れの残差 (`default_stress_residual`) が用いられる |
| Δλ の行 | 既定は `yield_function = 0`。粘塑性などで差し替える場合は、`self.dlambda(残差)` を返す |
| `Explicit` との混在 | 可能である |
| 行のスケール | 同じ根を持つ残差であれば、行のスケールを変えても接線は変わらない |

`inc.trial` と `inc.strain_inc` は、どの残差の行に、どの係数で用いてもよい。整合接線は、残差を Δε について自動微分して求められるため、`R_σ = σ − σ_trial + …` の形を仮定しない。次は、σ の行をスケールした例 (`J2IsotropicDamage`) である。

```{literalinclude} ../../src/returnmap/models/j2_damage/damage.py
:language: python
:pyobject: J2IsotropicDamage.state_residual
:dedent: 4
```

## 5. 流れ方向

既定の流れ方向は関連流れであり、`∂f/∂σ` を自動微分で求める。塑性ポテンシャルが降伏関数と異なる場合、または手導出の式を用いる場合は、`flow` を定義する。`default_stress_residual` は変更しない (残差の代数は変わらない)。

流れ方向には、せん断成分が 2 倍異なる 2 つの Voigt 表現がある。`flow` は配列ではなく `FlowVector` を返し、使用する側が表現を指定する ([構成と規約](../concepts.md) の 5 節)。

| 場面 | 書き方 |
|---|---|
| 工学ひずみのせん断で求めた方向を返す | `return self.strain_flow(n)` |
| テンソル成分 (教科書の形) で求めた方向を返す | `return self.stress_flow(n)` |
| ひずみとして用いる (`C @ n`、Δε_p) | `self.flow(state).strain_like` |
| 応力として用いる (背応力の発展則) | `self.flow(state).stress_like` |

配列をそのまま返した場合、および `FlowVector` に直接演算した場合は `TypeError` になる。モデルの中に ×2、÷2 の変換を書かない。

次は、有効応力空間の法線を返す例 (`J2IsotropicDamage`) である。

```{literalinclude} ../../src/returnmap/models/j2_damage/damage.py
:language: python
:pyobject: J2IsotropicDamage.flow
:dedent: 4
```

## 6. 滑らかな近似

モデルの式は自動微分され、反復で解かれる。絶対値、ノルム、`max`、条件分岐のような折れや段差がある点では、導関数が定まらず (`sqrt(0)` の微分は無限大になる)、反復が収束しない。`returnmap.utils` は、折れの無い近似を提供する。

| 関数 | 置き換える対象 | 式 |
|---|---|---|
| `smooth_sqrt(x)` | `sqrt(x)` | √(x + ε²) |
| `smooth_abs(x)` | `abs(x)` | √(x² + ε²) |
| `smooth_norm(v)` | ベクトルのノルム | √(v·v + ε²) |
| `smooth_direction(v)` | `v / ‖v‖` | v / `smooth_norm(v)` |
| `smooth_macaulay(x)` | `max(x, 0)` | (x + `smooth_abs(x)`) / 2 |
| `smooth_max(a, b)`、`smooth_min(a, b)` | `max`、`min` | — |
| `smooth_heaviside(x, beta=500)` | `x > 0` による切り替え | ½ (1 + tanh(βx/2)) |

使用例は、4 節の `J2IsotropicDamage.state_residual` (しきい値を超えた分を `smooth_macaulay` で求める) と、`YUKinematic` の硬化係数の切り替え (`smooth_heaviside`) である。

| 項目 | 内容 |
|---|---|
| `eps` | 0 での除算を防ぐ値 (既定 1e-12)。応答には現れない |
| `smooth_heaviside` の幅 | 引数の単位で決まる。既定の `beta=500` では、0 から 1 への遷移 (1 % → 99 %) の幅は、引数の値で ±0.01 である。引数が応力 (MPa) の場合は ±0.01 MPa であり、実質的に段差になる。ひずみのように小さい値の量を渡す場合は、`beta` を調整するか、引数に倍率を掛ける |
| 手書きの経路と UMAT | 同じ関数を、同じ係数で書く。係数が異なると、遷移の帯の中で 2 つの実装が異なる方程式を解く |
| 主値の符号による分解 | これらの関数ではなく、主値分解の関数を用いる (8 節) |

### 切り替えに掛ける量

切り替えに掛ける量が、切り替わる点で 0 にならない形は避ける。

| 形 | 例 | 切り替わる点での挙動 |
|---|---|---|
| 避ける形 | `R_new = R_n + smooth_heaviside(g) * dR` | ゲートが 0 → 1 に動くと、`R` が `dR` だけ不連続に変わる |
| 連続な形 | `r_new = r_n + self.h * smooth_macaulay(g)` | 掛ける量が `g` とともに 0 になる |

避ける形では、解が切り替えの縁に来るステップで、残差がほぼ不連続になる。その結果は次のとおりである。

- 反復が切り替えの両側を往復し、収束しないことがある。
- 収束した場合も、増分をわずかに変えると収束しなくなることがあり、接線の差分検査が成り立たないことがある。
- 平滑化の幅を広げても、段差の大きさは変わらない。

式の定義としてこの形を避けられないモデルもある (Yoshida-Uemori の限界曲面の膨張は、停留の判定がこの形になる)。その場合の扱いは次のとおりである。

- 切り替わるステップでは、解き直しや細かい刻みが必要になることがある。
- 切り替えを定数として微分した手導出は、そのステップでのみ近似になる。許容値を変更せず、`approx_when` で宣言する ([検査の一覧](verification.md) の 3 節)。

## 7. 弾塑性の判定

既定の弾性判定は `yield_function(試行状態) <= 0` である。古典的な塑性では、既定のままでよい。

| フック | 定義する場合 |
|---|---|
| `elastic_check(inc) -> bool` | 弾性域が降伏関数の内部と一致しないモデル。下負荷面モデルが典型である (判定が 2 条件の論理和になり、1 つの滑らかな関数の符号で表せない) |
| `elastic_update_state(inc)` | 弾性ステップでも変化する状態 (応力に追従する正規降伏比など) を持つモデル。弾性と判定された後に適用される |

次は、テスト用の下負荷面モデルの定義である (`tests/integration/simulation/integrator/test_base.py`)。

```{literalinclude} ../../tests/integration/simulation/integrator/test_base.py
:language: python
:pyobject: _SubloadingNatural.elastic_check
:dedent: 4
```

```{literalinclude} ../../tests/integration/simulation/integrator/test_base.py
:language: python
:pyobject: _SubloadingNatural.elastic_update_state
:dedent: 4
```

- 追従する状態を予測子 (8 節) ではなく `elastic_update_state` に書く理由: 試行状態に含めると `yield_function(試行状態)` が恒等的に 0 になり、既定の判定が載荷と除荷を区別できなくなる。
- 弾性ステップでの更新が応力に戻る発展則 (弾性エネルギーで駆動される損傷、クリープ) は扱えない。

## 8. 弾性予測子

既定の弾性予測子は線形であり、`σ_trial = σ_n + C(state_n)·Δε` である。次のモデルは、`elastic_predictor` を定義し、試行状態をモデル側で与える。

- 弾性則が区分線形であるモデル (片側損傷のように、主値の符号で剛性が変わる)。
- σ 以外の状態 (`Implicit` の有効応力など) が弾性的に変化するモデル。

```{literalinclude} ../../src/returnmap/models/j2_damage/unilateral.py
:language: python
:pyobject: J2IsotropicDamageUnilateral.elastic_predictor
:dedent: 4
```

| 項目 | 規定 |
|---|---|
| 戻り値 | `self.<名前>.update(値)` のリスト。`stress` は必須である。`Implicit` の状態を `self.<名前>(値)` と書くと、残差として扱われる |
| 試行状態の用途 | 戻り値が `inc.trial` になる。弾性ステップの解、弾性判定、塑性の反復の初期値、接線のすべてに用いられる |
| 微分可能性 | 自動微分で追跡できる形で書く |
| `elastic_update_state` との違い | 予測子は毎ステップ計算され、判定に用いられる。`elastic_update_state` は、弾性と判定された後の更新である |
| 定義しない場合 | 線形の予測子が用いられ、弾性接線は `elastic_stiffness(state_n)` である |

### 弾性接線

予測子を定義した場合、弾性接線は予測子の自動微分で求められる。UMAT に移植する式を検査する場合は、`elastic_tangent(inc)` に `dσ_trial/dΔε` を書く。

| integrator | 弾性接線 |
|---|---|
| `PythonAnalyticalIntegrator` | `elastic_tangent(inc)` を用いる (未定義の場合は `NotImplementedError`) |
| `PythonNumericalIntegrator` | 常に予測子を自動微分する |

両者を `crosscheck` で比較すると、手導出を検査できる。

### 主値分解

主値に依存する量には、`returnmap.utils` の主値分解の関数を用いる。`anp.linalg.eigh` の微分は、主値が重なる状態 (単軸、平面応力、純せん断) で NaN になるためである。

| 関数 | 内容 |
|---|---|
| `positive_part(sigma, dimension)`、`negative_part(...)` | 主値の正の部分、負の部分で作るテンソル σ₊ = Σ⟨σ_i⟩ n_i⊗n_i。接線も厳密である |
| `positive_part_sq(sigma, dimension)`、`negative_part_sq(...)` | 主値の正の部分、負の部分の 2 乗和 (Σ⟨σ_i⟩₊² など) |
| `macaulay_sq(x)`、`neg_macaulay_sq(x)` | スカラーに対する ⟨x⟩₊²、⟨x⟩₋² |

主値 0 の点は平滑化せず、規約 H(0) = c で扱う。既定は `DEFAULT_HEAVISIDE_AT_ZERO = 1/2` であり、暫定の値である (UMAT 側の規約に合わせて確定する)。上記の関数は `c=` で変更できる。参照実装は `src/returnmap/models/j2_damage/unilateral.py`、設計の理由は `utils/spectral.py` の docstring にある。

## 9. 手書きの解法

`user_defined_return_mapping(inc)` を定義すると、`PythonIntegrator` はそれを用いる。`None` を返した場合は、自動微分の解法が用いられる。整合接線の閉形式は `user_defined_tangent(new, dlambda, inc)` に書く。検査の手順は [手書きの解法を照合する](model-hand-written.md) にある。

次は `J2Isotropic3D` の閉形式の解である。

```{literalinclude} ../../src/returnmap/models/j2_isotropic/model.py
:language: python
:pyobject: J2Isotropic3D.user_defined_return_mapping
:dedent: 4
```

| 項目 | 規定 |
|---|---|
| 戻り値 | `ReturnMappingResult`。全フィールドが必須である |
| `converged` | 必ず書く。手書きの反復は、許容値で終了したか、反復の上限で終了したかを、この値で報告する |
| 閉形式の解 | `n_iterations=1`、`residual_history=[f_trial, 0.0]`、`converged=True` |
| 結果の検査 | Δλ < 0 と非有限の値は、`converged=True` と報告されていても非収束として扱われる |
| 予測子を定義したモデル | `inc.trial` は、予測子が与えた全フィールドを持つ |

## 関連ページ

- モデルを書いて検査する手順: [モデルを書く](model-python.md)
- 手導出ブロックと近似の宣言: [検査の一覧](verification.md)
- Voigt 規約と `FlowVector`: [構成と規約](../concepts.md)
