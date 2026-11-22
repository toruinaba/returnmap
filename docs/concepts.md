# 構成と規約

returnmap の構成 (層、パッケージ、integrator の種類) と、全体に共通する規約 (Voigt 規約、流れ方向の表現、許容値) を示す。非収束の扱いは [非収束の扱い](guide/nonconvergence.md) にある。

## 1. 三層構造

```text
MaterialModel  ──定義──▶  yield_function / update_state / state_residual / flow
      │
      ▼
  integrator   ──1 ステップ──▶  stress_update(strain_inc, stress_n, state_n) → StressUpdateResult
      │
      ▼
    driver     ──履歴──▶  run(load) → DriverResult
                          iter_run(load) → DriverStep の列
```

| 層 | 役割 | 主な名前 |
|---|---|---|
| モデル | 構成則。降伏関数、状態の発展則、パラメータを持つ | `MaterialModel` |
| integrator | 1 ステップの応力更新 (弾性予測 → 判定 → リターンマッピング → 整合接線)。Abaqus による UMAT の呼び出し 1 回に相当する | `PythonIntegrator`、`FortranIntegrator` |
| driver | integrator を載荷履歴に沿って呼ぶ | `StrainDriver`、`StressDriver`、`MixedDriver` |

フィッティング (`fit_params`)、比較 (`compare_models`)、検査 (`returnmap.verification`)、サービス (`returnmap.service`) は、この 3 層の上に作られている。

演算は倍精度 (float64) で行う。自動微分は autograd による。

## 2. パッケージ

公開の import 元は、次の 7 つのパッケージである。各名前は、1 つのパッケージのみから import する。

| パッケージ | 内容 |
|---|---|
| `returnmap.core` | モデルの定義に用いるもの: `MaterialModel`、`ParamMeta`、状態変数の宣言 (`Implicit`、`Explicit`、`NTENS`、`SCALAR`)、次元 (`SOLID_3D`、`PLANE_STRAIN`、`PLANE_STRESS`、`UNIAXIAL_1D`)、`hand_block`、`hand_tangent`、`verified_against_fortran`、結果の型 |
| `returnmap.models` | 組み込みモデル (J2、J2 + 損傷、Armstrong-Frederick、Chaboche、Ohno-Wang、Yoshida-Uemori) |
| `returnmap.simulation` | integrator、driver、載荷履歴 (`FieldHistory`)、`compare_models`、`FortranModule` |
| `returnmap.fitting` | `fit_params`、`ParameterSpace`、`FitTarget`、`Observation` |
| `returnmap.verification` | 検査: `check_model`、`check_hand_blocks`、`check_solution`、`check_tangent`、`crosscheck`、`compare_jacobians`、`check_bindings`、UMAT の宣言 (`Umat`) |
| `returnmap.service` | Web API 向けの契約層: `ModelRegistry`、`run_simulation` |
| `returnmap.utils` | モデル用の補助関数 (滑らかな近似、主値分解) |

サブモジュール (`returnmap.core.state`、`returnmap.simulation.driver` など) は実装の置き場所であり、文書化されたインターフェースではない。

## 3. integrator の種類

| クラス | 解法 | 用途 |
|---|---|---|
| `PythonIntegrator(model)` | モデルが手書きの解法 (`user_defined_return_mapping`) を持てばそれ、無ければ自動微分 | 通常の計算 |
| `PythonNumericalIntegrator(model)` | 自動微分に固定 | 検査 (手書きの解法の基準) |
| `PythonAnalyticalIntegrator(model)` | 手書きの解法に固定 | 検査 (自動微分との比較) |
| `FortranIntegrator.from_model(fortran, subroutine, model)` | コンパイルした UMAT | Python のモデルとの比較 |
| `SubSteppingIntegrator(integrator)` | 上記のいずれかを包み、収束しなかった増分を分割する | 大きい増分 ([非収束の扱い](guide/nonconvergence.md)) |

| 設定 (Python の integrator) | 既定 | 内容 |
|---|---|---|
| `nr_jacobian` | `"exact"` | `"broyden"` は、ステップの初回のみ厳密なヤコビアンを計算し、以降を Broyden 更新で代用する。残差が減少しない場合は、厳密なヤコビアンを再計算する。整合接線は常に厳密であり、精度には影響しない。手書きの解法を持たないモデルの長い履歴とフィッティングに向く (実測で 2〜4 倍速い。反復は 1〜2 回増える) |
| `nr_restart` | `True` | 自動微分の解法が失敗したステップを、別の初期値から 1 回解き直す。`False` は解き直しを行わない ([非収束の扱い](guide/nonconvergence.md)) |
| `raise_on_nonconverged` | `True` | `False` は、非収束を例外ではなく値で返す ([非収束の扱い](guide/nonconvergence.md)) |
| `tol`、`max_iter` | `1e-10`、`50` | 応力更新の反復の許容値と上限 |

応力更新の未知数は、常に `[σ, Δλ, implicit の状態]` である。`Explicit` の状態は `update_state` で求められ、未知数に含まれない。

## 4. Voigt 規約

3D の応力とひずみは 6 成分であり、順序は `[11, 22, 33, 12, 13, 23]` (ABAQUS の UMAT の規約) である。せん断成分の扱いは、量の種類によって異なる。

| 種類 | 量 | せん断成分 |
|---|---|---|
| 応力系 | σ、偏差応力 s、背応力 α など | テンソル成分そのもの |
| ひずみ系 | ε、ε_p、Δε | 工学ひずみ (γ12 = 2ε12)。ABAQUS の `DSTRAN`、弾性剛性 `σ = C:ε`、driver の入力は、すべてこの扱いである |

ノルムや相当量の計算では、Mandel のスケーリング (せん断 ×√2) が内部で適用される。次元ごとの演算子 (`dev`、`vonmises`、`inner_product`、`strain_norm` など) は `StressDimension` が持ち、`MaterialModel` の同名のメソッドから用いる。

| 次元 | 成分数 | 成分 |
|---|---|---|
| `SOLID_3D` | 6 | 11、22、33、12、13、23 |
| `PLANE_STRAIN` | 4 | 11、22、33、12 |
| `PLANE_STRESS` | 3 | 11、22、12 |
| `UNIAXIAL_1D` | 1 | 11 |

## 5. 流れ方向の表現

塑性流れの方向 n は、同じ方向に対して、せん断成分が 2 倍異なる 2 つの Voigt 表現を持つ。両者の取り違えは、単軸と二軸の計算では結果に現れない。そのため、`model.flow(state)` は配列ではなく `FlowVector` を返し、用いる側が表現を指定する。

| 式 | 表現 | 用途 |
|---|---|---|
| `model.flow(state).strain_like` | 工学ひずみのせん断 | Δε_p = Δλ·n、`C @ n`、残差のひずみの次元の項 |
| `model.flow(state).stress_like` | テンソル成分 | 背応力など、応力系の内部変数の発展則 |
| `model.flow(state) * 2.0` など | — | `TypeError`。表現を指定しない演算は拒否される |

`flow` を定義する方法は [モデルのフック](guide/writing-models.md) の 5 節、型の詳細は API リファレンスの [`FlowVector`](reference/core.md) を参照。

## 6. 許容値

`tol` は、その対象の許容値を表す。

| 対象 | `tol` の意味 |
|---|---|
| integrator | 応力更新の反復の収束。判定は `tol × max(1, ‖σ_trial‖_∞)` に対して行い、単位 (MPa、Pa) に依存しない |
| driver (`StressDriver`、`MixedDriver`) | 外側の反復の収束。判定は `tol + rtol·‖σ_target‖_∞` に対して行う |
| 検査 (`check_hand_blocks`、`check_solution`、`check_tangent`、`compare_jacobians`、`check_bindings`) | 合否。誤差の定義は全検査で共通である ([検査の一覧](guide/verification.md)) |
| `crosscheck` | 合否。`stress_tol`、`state_tol`、`tangent_tol` で指定する。driver の外側の反復の設定は、`driver_factory` が作る driver が持つ |
| `fit_params` | 最適化の収束 |
