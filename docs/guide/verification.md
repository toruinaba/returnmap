# 検査の一覧

`returnmap.verification` の検査と、誤差の測り方を示す。モデルを段階ごとに検査する手順は [新しいモデルを足す](adding-a-model.md) にあり、このページは個々の検査の仕様を扱う。

検査はすべて関数であり、報告を返す。報告は、合否 `report.passed` と、内訳の表 `report.summary()` を持つ。テストは、どの検査でも `assert report.passed, report.summary()` の 1 行になる。

| 関数 | 比較するもの | 節 |
|---|---|---|
| `check_model(model)` | 以下の検査のうち、そのモデルに該当するものすべて | 1 |
| `check_hand_blocks(model, results)` | 手導出したヤコビアンのブロックと、自動微分 | 3 |
| `check_tangent(integrator, strain_inc, stress_n, state_n)` | 整合接線と、中心差分 | 4 |
| `check_solution(model, results)` | 手書きのリターンマッピングの解と、自動微分側の残差 | 5 |
| `crosscheck(integrator_a, integrator_b, load)` | 同じ載荷履歴を進めた 2 つの実装 | 6 |
| `compare_jacobians(model, result_a, result_b)` | 同じステップの 2 つの解における残差のヤコビアン | 6 |
| `check_bindings(model, fortran, cases)` | Fortran に対応づけたメソッドと、そのサブルーチン | 7 |

コードは [`examples/example_checks.py`](https://github.com/toruinaba/returnmap/blob/main/examples/example_checks.py) の各節と、組み込みモデルのソースである。

## 誤差の測り方

誤差の定義は、全検査で共通である。

| 項目 | 定義 |
|---|---|
| 誤差 | `max|候補 − 参照| / max|参照|`。配列 (ブロック、接線、応力ベクトル) 全体で測る |
| 合格 | `誤差 <= tol` |
| 成分ごとに割らない理由 | 厳密には 0 である成分に乗った丸め誤差を、100 % の誤差と数えないため |
| 分母に単位を持つ定数を置かない理由 | モデルを MPa で書いても Pa で書いても、同じ許容値が同じ意味を持つようにするため |
| NaN | 最悪の誤差として扱う |

比較する相手の配列が無い量は、次の尺度に対して測る。

| 量 | 尺度 |
|---|---|
| `check_solution` の残差の行 | その行を構成する項の大きさ (5 節) |
| 厳密に 0 のヤコビアンのブロック | その行と列から決まる大きさ (行を構成する項の大きさ ÷ 未知数の大きさ)。これより 12 桁以上小さいブロックを 0 と見なす |
| 0 の接線 (降伏棚、1 次元の完全塑性) | 弾性剛性。差分の接線が弾性剛性より 8 桁以上小さい場合に適用する |

## 1. まとめて検査する

`check_model` は、モデルが宣言しているものに応じて、該当する検査をすべて実行する。

```{literalinclude} ../../examples/example_checks.py
:language: python
:start-after: "# %% 1. まとめて検査する"
:end-before: "# %%"
```

```none
YUKinematic1D: every check that applies passes
  ok   converges                 71 steps over 2 histories, 57 plastic
  ok   tangent                   4.0e-09 (tol 1e-05) at 8 plastic steps
  —    hand blocks               no @hand_block on the model
  —    hand-written solution     no user_defined_return_mapping on the model
  —    hand-written vs autodiff  no user_defined_return_mapping on the model
  —    hand-written tangent      no user_defined_return_mapping on the model
  —    UMAT                      no UMAT declared for YUKinematic1D in returnmap.models.yu_kinematic.verification
  histories: uniaxial (24 plastic steps), reversed 2 % (33 plastic steps)
```

| 行 | 検査の内容 | 対象になる条件 |
|---|---|---|
| `converges` | 自動微分の経路が、すべての載荷履歴で最後まで収束する (Δλ ≥ 0、有限の値、塑性ステップがある) | 常に |
| `tangent` | 整合接線と中心差分の一致 (`check_tangent`、履歴ごとに数ステップ) | 常に |
| `hand blocks` | 手導出ブロックと自動微分の一致 (`check_hand_blocks`) | `@hand_block` がある |
| `hand-written solution` | 手書きの解が、自動微分側の残差を満たす (`check_solution`) | `user_defined_return_mapping` がある |
| `hand-written vs autodiff` | 手書きの経路と自動微分の経路の一致。応力・状態・接線を全ステップで比較する (`crosscheck`) | 同上 |
| `hand-written tangent` | 手書きの接線と中心差分の一致 | `user_defined_tangent` がある |
| `UMAT` | コンパイルした UMAT と Python のモデルの一致。同じ履歴を全ステップで比較する (`crosscheck`) | UMAT が宣言されており、ビルドしてある |

| 項目 | 内容 |
|---|---|
| `—` の行 | 宣言していない項目。不合格には数えない。自動微分で解くだけのモデルは、上の 2 行で検査され、合格する |
| 行の内訳 | `report["hand blocks"].reports` などから、元の検査の報告を取り出せる ([手書きの解法を照合する](model-hand-written.md) の 4 節) |
| UMAT の比較相手 | `PythonIntegrator(model)`。手書きのリターンマッピングがあればそれ、無ければ自動微分である |
| 未ビルドの UMAT | 行は `not built` になり、不合格には数えない ([UMAT を確かめる](model-umat.md)) |
| 検査できないこと | 定式化が意図どおりであるか。どの行も、同じモデルの 2 つの導出を比較する。全体が共有する仮定 (単位、符号の規約、定式化) は検出できない |

### 載荷履歴

履歴を渡さない場合、モデルの降伏ひずみに合わせた標準の履歴を用いる。

| 履歴 | 内容 |
|---|---|
| `uniaxial` | 単軸の引張 → 除荷 → 圧縮 |
| `multiaxial` | 直応力成分とせん断成分を同時に動かし、途中で方向を変える (多成分のモデルのみ)。転置した外積や、せん断の ×2 の誤りは、単軸の履歴では検出できない |

標準の履歴では機構が働かないモデルは、ファミリの `verification.py` に `histories(model)` を定義し、履歴を追加する。`check_model` は、モデルのモジュールの隣にあるこのモジュールを自動で読む (登録は不要である)。

```{literalinclude} ../../src/returnmap/models/yu_kinematic/verification.py
:language: python
:pyobject: histories
```

| 場合 | 方法 |
|---|---|
| 履歴の追加 | `histories(model)` が、名前 → 載荷履歴の dict を返す |
| 標準の履歴の置き換え | 標準と同じ名前 (`uniaxial`、`multiaxial`) で宣言する |
| 弾性域の無いモデル (下負荷面モデルなど) | 降伏ひずみが定まらないため、尺度のひずみを渡して標準の履歴を作る: `standard_histories(model, eps_y=...)` |
| 履歴を引数で渡す | `check_model(model, histories)`。渡した履歴のみを用いる |

`verification.py` の `MODELS` は、そのファミリのテスト (`tests/test_verification.py`) が `check_model` にかけるモデルの一覧である。

## 2. 残差のヤコビアン

`integrator.jacobian(result)` は、1 ステップの解における残差のヤコビアン (自動微分) を、名前付きのブロックで返す。手導出の参照として用いる。

```{literalinclude} ../../examples/example_checks.py
:language: python
:start-after: "# %% 2. 残差のヤコビアン"
:end-before: "# %%"
```

```none
['stress', 'dlambda', 'alpha', 'ep']
(6, 6)
[ 0.961 -0.48  -0.48   0.48   0.     0.   ]
```

行と列は、応力更新の未知数 `[σ, Δλ, implicit の状態]` である。`jac.part[行][列]` が 1 つのブロックを返す。

## 3. 手導出ブロックの宣言と照合

手導出したブロックは、モデル上で `@hand_block(行, 列)` を付けて宣言する。メソッド自体は変わらず、直接呼ぶことも、`verified_against_fortran` と重ねることもできる。次は `YUKinematicPS` の宣言である。

```{literalinclude} ../../src/returnmap/models/yu_kinematic/model.py
:language: python
:start-at: '@hand_block("theta", "stress")'
:end-at: "return - 2 / 3 * C_k * a * dlambda * np.eye(3)"
:dedent: 4
```

| 項目 | 内容 |
|---|---|
| 引数 | 名前で渡される: `state` (評価点の状態)、`dlambda`、`state_n`、`inc`。Fortran の引数リストに合わせた手書きのまま宣言できる |
| 戻り値 | そのブロックの配列 |

### 近似の宣言

意図した近似は、モデル上で宣言する。テストの許容値は変更しない (理由が残らないためである)。

| 宣言 | 意味 | 用途 |
|---|---|---|
| `approx=2e-2` | 全ステップで、自動微分との差をこの値まで許容する | 導出から省いた項が、どのステップにもある場合 |
| `approx_when="メソッド名"` | そのメソッドが真を返すステップのみを近似とする。他のステップは厳密な許容値で判定する | 切り替えを定数として微分している場合 (切り替わるステップのみ一致しない) |
| `approx_when` と `approx=` の併用 | 条件に該当するステップを、比較から除かず、その値で判定する | — |
| `note=` | 理由。どちらの宣言でも必須である | — |

`YUKinematicPS` は、停留曲面のゲートが切り替わるステップを近似として宣言している。条件のメソッドは、ブロックと同じく引数を名前で受け取り、真偽を返す。

```{literalinclude} ../../src/returnmap/models/yu_kinematic/model.py
:language: python
:start-at: "_GATE_NOTE = ("
:end-at: ")"
```

```{literalinclude} ../../src/returnmap/models/yu_kinematic/model.py
:language: python
:start-at: '@hand_block("theta", "beta",'
:end-at: "return 2 / 3 * C_k * a * dlambda * np.eye(3)"
:dedent: 4
```

```{literalinclude} ../../src/returnmap/models/yu_kinematic/model.py
:language: python
:pyobject: YUKinematic.stagnation_gate_in_transition
:dedent: 4
```

手書きの接線 (`user_defined_tangent`) の近似は、`@hand_tangent` で同様に宣言する。

```{literalinclude} ../../src/returnmap/models/yu_kinematic/model.py
:language: python
:start-at: '@hand_tangent(approx_when="stagnation_gate_in_transition", note=_GATE_NOTE)'
:end-at: "return self.calc_ddsdde("
:dedent: 4
```

接線の宣言を参照するのは `check_model` である。`check_tangent` と `crosscheck` は、integrator を受け取るためモデルの宣言を参照せず、渡された許容値で判定する。

### 照合

`check_hand_blocks` は、宣言済みのブロックを、履歴全体で自動微分と照合する。

```{literalinclude} ../../examples/example_checks.py
:language: python
:start-after: "# %% 3. 手導出ブロックの照合"
:end-before: "# %%"
```

出力の抜粋を示す (全 16 行のうち 3 行)。

```none
hand blocks: 16/16 derived, checked at 21 plastic steps — all within tolerance
  ok   [theta][beta] calc_ft_fb: 2.042e-17 (tol 1e-08, worst at point 1)  (approx at 3 of 21 points, 1.0e+00 there: the stagnation gate is differentiated as a constant; exact while it is 0 or 1, not on the step where it switches)
  ok   [dlambda][dlambda] calc_fy_fl: 0.000e+00 (tol 1e-08, worst at point 1)
```

| 項目 | 内容 |
|---|---|
| 集計 | 最大の誤差をブロックごとに集計する。1 つの誤りが他のブロックに隠れない |
| ブロックの誤差 | 履歴全体での最大差を、そのブロックの履歴全体での最大値で割った値。符号の誤りは 2、0 として省いたブロックは 1 になる |
| `approx_when` のステップ | 比較から除き、数を報告する (上の例では 21 ステップ中 3 ステップ) |
| 未導出のブロック | `not derived` として一覧に出る |
| 履歴 | 直ひずみとせん断を同時に動かすものを用いる。外積の転置や、せん断の ×2 の取り違えは、単軸や純せん断では検出できない |

`hand_jacobian(model, result)` (`returnmap.simulation`) は、宣言済みのブロックを手導出、残りを自動微分で埋めた行列を返す。

## 4. 接線の差分検査

`check_tangent` は、整合接線を、応力更新の中心差分と比較する。

```{literalinclude} ../../examples/example_checks.py
:language: python
:start-after: "# %% 4. 接線の差分検査"
:end-before: "# %%"
```

```none
tangent check: agrees with central differences
  ok   ddsdde: 8.525e-11 (tol 1e-05, worst at [0][1])
(6, 6)
```

| 項目 | 内容 |
|---|---|
| 引数 | `stress_update(strain_inc, stress_n, state_n)` と同じ順である |
| 対象 | `stress_update` を持つ integrator。手書きの接線 (`PythonAnalyticalIntegrator`) と UMAT (`FortranIntegrator`) にも用いられる |
| `report.rel_err_matrix` | 成分ごとの誤差 |

## 5. 手書きの解の検査

`check_solution` は、手書きのリターンマッピング (`user_defined_return_mapping`) が返した解を、自動微分側の残差に代入する。

```{literalinclude} ../../examples/example_checks.py
:language: python
:start-after: "# %% 5. 手書きの解の検査"
:end-before: "# %%"
```

```none
solution check: 21 plastic steps — every row satisfied
  ok   [stress]: 9.404e-17 (tol 1e-08, worst at point 22)
  ok   [dlambda]: 1.712e-16 (tol 1e-08, worst at point 20)
  ok   [theta]: 1.444e-13 (tol 1e-08, worst at point 16)
  ok   [beta]: 8.023e-17 (tol 1e-08, worst at point 16)
  ok   R (explicit state): 0.000e+00 (tol 1e-08, worst at point 1)
  ok   q (explicit state): 0.000e+00 (tol 1e-08, worst at point 1)
  ok   r (explicit state): 0.000e+00 (tol 1e-08, worst at point 1)
  ok   eps_eq (explicit state): 0.000e+00 (tol 1e-08, worst at point 1)
  ok   theta_max (explicit state): 0.000e+00 (tol 1e-08, worst at point 1)
```

| 項目 | 内容 |
|---|---|
| 判定の単位 | 残差の行ごと。手書きの残差の行が自動微分側と異なる場合、その行のみが 0 にならず、原因の行を特定できる |
| ステップの扱い | 各ステップを個別に判定する。誤差の累積の影響を受けない |
| 比較の対象 | 根のみ。行のスケールを変えた正しい書き換えは合格する。残差ベクトルを持たない閉形式の解にも用いられる |
| explicit の状態 | 報告された値を、`update_state` の結果と比較する |
| 動かなかった状態 | 履歴中に一度も変化しなかった状態は、`not exercised` として一覧に出る (その更新式は検査されていない) |
| 誤差の尺度 | 行を構成する項の大きさに対する相対値。ソルバが報告した残差ノルムの 10 倍までを許容する (閉形式の解には許容幅を与えない) |
| 非収束のステップ | 渡すと `ValueError` になる。検査の対象は解のみである |
| 分割されたステップ | `SubSteppingIntegrator` が分割したステップ (`n_substeps > 1`) は対象から除き、数を `left out` として報告する |

## 6. 2 つの実装の比較

`crosscheck` は、2 つの integrator を同じ載荷履歴で別々に進め、ステップごとに応力・状態・接線を比較する。用途は、自動微分と手書きの経路の比較、Python のモデルと UMAT の比較である。

```{literalinclude} ../../examples/example_checks.py
:language: python
:start-after: "# %% 6. 2 つの実装の比較"
```

```none
crosscheck: 35 steps — all agree
  ok   stress: 1.492e-16 (tol 1e-06)
  ok   state[ep]: 3.837e-16 (tol 1e-06)
  ok   tangent: 4.036e-16 (tol 1e-05)
```

```text
None
```

| 項目 | 内容 |
|---|---|
| driver | 載荷履歴から決まる。ひずみ履歴は `StrainDriver`、応力履歴は `StressDriver` である |
| `driver_factory` | 上記以外の driver を用いる場合に渡す。3D のモデルの単軸応力は `driver_factory=lambda i: MixedDriver(i, prescribed_strain_idx=[0])`。driver の設定 (外側の反復の許容値など) も、この関数が持つ |
| `report.first_failure` | 一致しなかった最初のステップ。2 つの結果 (`result_a`、`result_b`) を持つ。全ステップが一致した場合は `None` |
| `compare_jacobians(model, result_a, result_b)` | 2 つの解における残差のヤコビアンを、ブロックごとに比較する |

2 つの経路は独立に進むため、最初の不一致より後のステップは、すべて不一致になる。不一致の調べ方は [UMAT を確かめる](model-umat.md) の 6 節にある。手書きの解のどの行が原因かを特定する場合は、5 節の `check_solution` を用いる。

## 7. Fortran との対応の検査

`check_bindings` は、`verified_against_fortran` で Fortran のサブルーチンに対応づけたメソッドを、同じ引数で両側から呼び、戻り値を比較する。

| 引数 | 内容 |
|---|---|
| `model` | 対応づけを持つモデル |
| `fortran` | コンパイルしたモジュール (`FortranModule(モジュール名)`) |
| `cases` | メソッド名 → (Python 側の引数の組, Fortran 側の引数の組)。例: `{"elastic_stiffness": ((model.initial_state(),), (model.E, model.nu))}` |
| `tol` | 許容値。既定は `1e-10` |

対応づけたメソッドのうち、`cases` に含めなかったものは、`not checked` として一覧に出る。入力を手で用意する必要があるため、この検査は `check_model` に含まれない。

## 関連ページ

- 段階ごとの手順: [新しいモデルを足す](adding-a-model.md)
- フックの仕様: [モデルのフック](writing-models.md)
- 非収束の扱い: [概念](../concepts.md)
