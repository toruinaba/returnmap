# 検査の一覧

`returnmap.verification` の検査と、誤差の測り方を示す。モデルを段階ごとに検査する手順は [新しいモデルを足す](adding-a-model.md) にあり、このページは個々の検査の仕様を扱う。

検査はすべて関数であり、報告を返す。報告は、合否 `report.passed` と、内訳の表 `report.summary()` を持つ。テストは、どの検査でも `assert report.passed, report.summary()` の 1 行になる。

| 関数 | 比較するもの | 節 |
|---|---|---|
| `check_model(model)` | 以下の検査のうち、そのモデルに該当するものすべて | 1 |
| `check_range(model)` | 同じ検査を、パラメータの範囲の端と、乱数で作った載荷履歴で行う | 8 |
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
| 差分が分解できない大きさの接線 (ほぼ飽和した硬化) | 差分の丸め誤差 ÷ 許容値。差分の丸め誤差は、差分を取った応力の最大値の機械イプシロン 2 個分を、ひずみの刻みで割った値である。接線がこれより小さい場合、誤差は差分の丸め誤差の範囲内であれば合格になる |

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
  —    UMAT entry                no UMAT declared for YUKinematic1D
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
| `UMAT entry` | ABAQUS が呼ぶ入口 (`umat`) と、それが呼ぶ中核のサブルーチンの一致。同じ状態と増分を 1 ステップずつ与える | 同上。入口が宣言されている (既定は `umat`) |

| 項目 | 内容 |
|---|---|
| `—` の行 | 宣言していない項目。不合格には数えない。自動微分で解くだけのモデルは、上の 2 行で検査され、合格する |
| 行の内訳 | `report["hand blocks"].reports` などから、元の検査の報告を取り出せる ([手書きの解法を照合する](model-hand-written.md) の 4 節) |
| UMAT の比較相手 | `PythonIntegrator(model)`。手書きのリターンマッピングがあればそれ、無ければ自動微分である |
| UMAT の許容値 | 比較相手が手書きの解法の場合は `1e-10` (接線は `1e-9`)。同じ代数の移植であり、6 桁の一致では、単精度で書いた定数 (差は 1e-7 程度) を検出できない。比較相手が自動微分の場合は `stress_tol` (接線はその 10 倍) |
| 入口の比較相手 | 同じ UMAT の中核である。不合格の行から、誤りが解法にあるのか、`PROPS` と `STATEV` の展開にあるのかを区別できる |
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

## 8. パラメータの範囲と乱数の載荷

`check_range` は、`check_model` と同じ比較を、パラメータと載荷の 2 つの方向に広げて行う。`check_model` は 1 組のパラメータと、降伏ひずみに合わせた載荷履歴で検査するため、それ以外の条件でだけ現れる誤りを検出しない。

```{literalinclude} ../../examples/example_checks.py
:language: python
:start-after: "# %% 8. パラメータの範囲と乱数の載荷"
```

```none
J2Isotropic1D over its range: every case passes
  parameters — check_model at each end of every documented range
    ok   E=1000              5 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent)
    ok   E=1e+06             5 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent)
    ok   nu=0                5 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent)
    ok   nu=0.49             5 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent)
    ok   sigma_y0=1          5 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent)
    ok   sigma_y0=5000       5 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent)
    ok   H=0                 5 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent)
    ok   H=100000            5 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent)
  loading — random strain histories
    ok   random              36 of 36 increments solved, 16 plastic; solution 4.5e-16, hand-written vs autodiff 3.3e-16, tangent 4.2e-10
    ok   random large        24 of 24 increments solved, 23 plastic; solution 1.6e-15, hand-written vs autodiff 1.5e-15, tangent 2.3e-09
    ok   random single-step  6 of 6 increments solved, 6 plastic; solution 1.5e-14, hand-written vs autodiff 9.5e-15, tangent 2.6e-08
1.1e-12 (tol 1e-05) at 4 plastic steps
```

| 組 | 内容 | 既定 |
|---|---|---|
| `parameters` | パラメータを 1 つずつ、`param_meta` の範囲 (`bounds`) の下限と上限に置いたモデル。それぞれを `check_model` にかける | `parameter_variants(model)` |
| `loading` | 乱数で作ったひずみ履歴。全成分を動かし、6 つの点を直線で結ぶ | `random_histories(model)` |

乱数の履歴は、次の 3 つである。大きさの基準は、`check_model` がそのモデルに用いる履歴 (標準の履歴と、ファミリが宣言した履歴) の最大のひずみ成分である。

| 履歴 | 点の大きさ | 1 区間の増分の数 |
|---|---|---|
| `random` | 基準と同じ | 6 |
| `random large` | 基準の 4 倍 | 4 |
| `random single-step` | 基準の 20 倍 | 1 |

パラメータの組は、`check_model` と同じ基準で合否を決める。載荷の組は、誤った答えだけを不合格にする。

| 載荷の組の結果 | 扱い |
|---|---|
| 手書きの解が残差の根でない (`check_solution`) | 不合格 |
| 同じ状態から、手書きの経路と自動微分の経路がともに収束し、結果が異なる。応力、状態、接線を、ともに収束した全ステップで比較する | 不合格 |
| 接線が中心差分と一致しない (履歴ごとに数ステップ) | 不合格 |
| UMAT が、同じ状態と増分に対して、Python のモデルと異なる結果を返す。または、Python のモデルが解けない増分で収束する | 不合格 |
| ABAQUS の入口が、UMAT の中核と異なる結果を返す。または、解けなかった増分で、収束するか `STRESS`・`STATEV` を変更する | 不合格 |
| モデル自身の経路 (手書きの解法があればそれ、無ければ自動微分) が解けない増分 | 履歴から除き、直前の解けた状態から続ける。件数を報告する。解けない増分には、解法が例外を送出する増分 (特異なヤコビアン、算術の例外、モデル自身の反復の失敗を表す `RuntimeError`) を含む |
| 手書きの経路が解け、自動微分の経路が解けないステップ | 件数を報告する。手書きの解は残差に代入して確かめる |
| 差分を取るための応力更新が解けないステップ | その点の接線の比較を省き、件数を報告する |
| 手書きの接線を、宣言された近似のステップで比較する場合 | `check_model` と同じ扱いである (3 節)。差分を取る全ステップが近似に当たる場合は、不合格になる |
| 塑性の増分を 1 つも解けない履歴 | `—` (`not solved`)。その履歴は不合格に数えない |
| すべての履歴が `not solved` | 不合格。乱数の載荷では何も比較できていない |

| 項目 | 内容 |
|---|---|
| 報告 | `report.passed`、`report.summary()`。`report[名前]` が 1 つの組を返し、`.reports` に元の検査の報告を持つ |
| UMAT | 宣言されており、ビルドしてある場合に、すべての組で比較する。未ビルドの場合は、載荷の組の行に `UMAT not built` と表示し、Python の側だけで合否を決める |
| 実行時間 | パラメータの組 1 つにつき `check_model` 1 回分。完成したモデルに対して実行する検査であり、毎回のテストには含めない |
| 乱数の種 | `check_range(model, seed=...)`。同じモデルと種に対して、同じ履歴になる |
| パラメータの組の例外 | `check_model` が解法の例外で終わった組は、`could not be checked` として不合格になる |
| プログラムの誤り | `TypeError`、`AttributeError` などは、解けない増分として扱わず、そのまま送出する |
| 差分の参照の精度 | 差分を取る応力は、既定の許容値で解く。許容値は試行応力に比例するため、降伏ひずみの数十倍の増分では、正しい接線と差分の差が `1e-5` 程度になり、不合格になる場合がある (`YUKinematic1D` の `random large` の 1 点) |
| 接線の差分の刻み | `check_model` と `check_range` は、`1e-7` と、ステップの応力が弾性で対応するひずみ (`max|σ| / max|C|`) の 1 万分の 1 のうち、小さいほうを用いる。`max|C|` は弾性剛性の最大成分で、ヤング率の 1.3 倍程度である。降伏ひずみが 1e-3 程度のモデルでも、応力が小さいステップでは `1e-7` より小さくなる |

パラメータの範囲は、パラメータごとに独立に与えられている。1 つだけを端に置くと意味を失うモデル (境界面が降伏面の内側に入る、など) は、検査する組をファミリの `verification.py` に `variants(model)` として宣言する。`check_range` は、宣言があればそれを `parameter_variants(model)` の代わりに用いる。

| 場合 | 方法 |
|---|---|
| 組の宣言 | `variants(model)` が、名前 → モデルの dict を返す。`parameter_variants(model)` の結果から、意味を持つ組を選ぶ |
| 組を引数で渡す | `check_range(model, variants=...)`。空の dict を渡すと、パラメータの組を省く |
| 履歴を引数で渡す | `check_range(model, histories=...)`。空の dict を渡すと、載荷の組を省く |
| 範囲を持たないモデル | `param_meta` に `bounds` が無い場合、パラメータの組は無く、報告の `notes` にその旨が入る |

## 関連ページ

- 段階ごとの手順: [新しいモデルを足す](adding-a-model.md)
- フックの仕様: [モデルのフック](writing-models.md)
- 非収束の仕組み: [非収束の扱い](nonconvergence.md)
