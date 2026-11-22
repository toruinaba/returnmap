# フィッティングの設定

`fit_params` の探索空間、ターゲット、比較の方法、driver の指定を示す。試験の記録からパラメータを同定する手順は [試験データに合わせる](tutorial.md) にあり、このページは個々の設定を扱う。

コードは [`examples/example_fitting.py`](https://github.com/toruinaba/returnmap/blob/main/examples/example_fitting.py) の各節であり、上から順に実行できる (約 5 秒)。記録は、既知のパラメータの J2 モデル (`sigma_y0 = 250`、`H = 1000`) の応答にノイズを加えたものである。

## 1. 基本の呼び出し

呼び出しの形は 1 つであり、モデル、ターゲットの列、`ParameterSpace` を渡す。

```{literalinclude} ../../examples/example_fitting.py
:language: python
:start-after: "# %% 1. 基本の呼び出し"
:end-before: "# %%"
```

```none
{'E': 210000.0, 'nu': 0.3, 'sigma_y0': 250.1, 'H': 994.0}
残差 0.148、収束 True、評価 39 回、失敗 0 回
J2Isotropic1D
```

| 引数 | 内容 |
|---|---|
| 第 1 引数 (モデル) | テンプレート。クラス、次元、パラメータ以外の設定 (`J2IsotropicDamageUnilateral` の `c` など) を与える。各評価のモデルは `model.with_params(**params)` で作られ、テンプレートのパラメータの値は用いられない |
| 第 2 引数 (ターゲット) | `FitTarget` の列。`FitTarget.stress_strain(strain, stress)` は、ひずみ履歴と応力の記録からターゲットを 1 つ作る |
| 第 3 引数 (`ParameterSpace`) | モデルの全パラメータを、`free`、`fixed`、`derived` のいずれかで与える |
| `method` | `"L-BFGS-B"` (既定)、`"Nelder-Mead"`、`"differential_evolution"` |
| `max_iter`、`tol` | 最適化の設定。応力更新の反復の設定ではない |
| `nr_jacobian="broyden"` | 手書きの解法を持たないモデルの評価を速くする。既定は `"exact"` |
| `driver_factory` | driver の指定 (5 節) |
| `progress` | 評価ごとに呼ばれる関数。引数は (評価回数, その残差, 最良の残差) |

| 結果 | 内容 |
|---|---|
| `result.model` | 同定したパラメータを入れたモデル。クラス、次元、設定はテンプレートと同じである |
| `result.model_params` | モデルのパラメータ (`param_names` と同じキー) |
| `result.params` | 探索空間のすべての名前と値。`derived` の計算に用いた補助変数を含む (2 節) |
| `result.residual` | 最終の残差 |
| `result.success`、`result.n_evaluations` | 最適化の収束、評価回数 |
| `result.n_failed` | 非収束であった評価の数 (6 節) |

### 範囲

`free` の範囲 (上下限) は、探索の尺度を兼ねる。`"L-BFGS-B"` と `"Nelder-Mead"` は、各変数を範囲の幅で割った変数で探索する。硬化係数 (1e4 程度) と回復係数 (1e2 程度) のように、桁の異なるパラメータを同じ重みで動かすためである。

- 範囲には、ありうる値の幅を有限の上下限で与える。
- 上下限の無い変数は、初期値の大きさで割られる。
- 結果と履歴は、与えた単位のまま返る。

## 2. パラメータ間の制約

`derived` は、他の変数からパラメータを計算する関数を与える。`Y < B` のような順序の制約は、比を探索の変数とし、下限を 1.0 とすることで、不等式制約を用いずに表現できる。

```{literalinclude} ../../examples/example_fitting.py
:language: python
:start-after: "# %% 2. パラメータ間の制約"
:end-before: "# %%"
```

```text
['E', 'H', 'H_ratio', 'nu', 'sigma_y0']
['E', 'H', 'nu', 'sigma_y0']
```

```none
H_ratio 3.98、H 994.0
```

`result.params` は補助変数 (`H_ratio`) を含むため、`model.with_params(**result.params)` は `TypeError` になる。モデルの構築には `result.model` または `result.model_params` を用いる。

## 3. 複数の試験

ターゲットは 2 階層で指定する。階層を分ける理由は、計算コストの違いである。

| 階層 | 対応するもの | コスト |
|---|---|---|
| `FitTarget` | 1 回の解析 (載荷履歴 1 本と driver 1 つ) | 目的関数の評価 1 回につき、解析 1 回 |
| `Observation` | 1 本の記録 | 同じターゲットの複数の `Observation` は、解析の結果を共有する。追加のコストはほぼ無い |

```{literalinclude} ../../examples/example_fitting.py
:language: python
:start-after: "# %% 3. 複数の試験"
:end-before: "# %%"
```

```none
{'E': 210000.0, 'nu': 0.3, 'sigma_y0': 250.0, 'H': 1000.0}
```

| `Observation` の引数 | 内容 |
|---|---|
| `values` | 記録の値 |
| `field` | 比較する量。既定は応力。`"strain"`、または状態量の名前 (`"ep"` など) |
| `component` | テンソルの成分 |
| `weight` | 重み。量の桁の違い (応力とひずみなど) を調整する |
| `metric`、`x_values` | 比較の方法 (4 節) |

- 各 `Observation` の残差二乗和は、データ点数で割られる。点数の多い記録が支配することを防ぐためである。
- 既定の比較は、点対点の二乗平均である。記録と計算のサンプリングが一致している必要がある。
- ターゲットごとに driver を変える場合は、`FitTarget(..., driver_factory=...)` で指定する (5 節)。

## 4. サンプリングが一致しない記録

記録が細かく計測されており、計算を粗い増分で行う場合は、`metric="pcm"` を用いる。Partial Curve Mapping は、点対点ではなく、2 本の曲線の間の面積で比較する。点数と横軸の範囲は異なってよい。ヒステリシスループのような多価の曲線にも用いられる。

```{literalinclude} ../../examples/example_fitting.py
:language: python
:start-after: "# %% 4. サンプリングが一致しない記録"
:end-before: "# %%"
```

```none
H 1000
```

| 項目 | 内容 |
|---|---|
| 粗い増分による偏り | 粗い増分は、降伏の角を弦で切る。複数のパラメータの同時フィッティングでは、結果に偏りが生じる (実測: 30 増分で `H` に 12 % の誤差、400 増分で消える)。単一のパラメータの探索は偏らない |
| 記録のノイズ | PCM は、2 本の曲線を弧長で対応づける。ノイズのある点列を結んだ折れ線は、元の曲線より長い。そのため、記録の弧長は、点列から見積もったノイズの大きさの 10 倍の長さの尺度で測る。記録は間引かず、全点と比較する。フィルタや間引きを事前に行う必要は無い |
| ノイズによる限界 | 尺度は降伏の角や載荷反転の先端をまたぐため、弧長はわずかに短く評価される (応力範囲の 0.4 % のノイズで 1 % 前後)。真のパラメータでも PCM の値は 0 にならず、0.01〜0.05 程度になる。これより小さい差は判別できない |

詳細は `fitting/metrics.py` の docstring にある。

## 5. driver の指定

`driver_factory` は、integrator を受け取り、driver を返す関数である。

| 指定 | 動作 |
|---|---|
| 省略 | ターゲットの載荷履歴から決まる。ひずみ履歴は `StrainDriver`、応力履歴は `StressDriver` (`crosscheck`、`compare_models` と同じ規則) |
| `fit_params(..., driver_factory=...)` | 全ターゲットに適用される |
| `FitTarget(..., driver_factory=...)` | そのターゲットのみに適用される |

優先順位は、ターゲットの指定、`fit_params` の指定、載荷履歴から決まるもの、の順である。

次は、3D のモデルを単軸応力の記録に合わせる例である ([応答を計算する](simulation.md) の 5 節)。

```{literalinclude} ../../examples/example_fitting.py
:language: python
:start-after: "# %% 5. driver の指定"
```

```none
{'E': 210000.0, 'nu': 0.3, 'sigma_y0': 250.1, 'H': 994.0}
```

## 6. 非収束の評価

探索空間は、応力更新が解を持たないパラメータの値を含む。評価中の非収束は例外にならず、有限の大きな値 (1e30) として採点され、最適化は継続する。`inf` を用いないのは、`"L-BFGS-B"` の差分勾配を壊さないためである。

| 名前 | 内容 |
|---|---|
| `result.n_failed` | 非収束であった評価の数 |
| `result.failure_message`、`result.failure_params` | 最初の失敗の内容と、そのときのパラメータ |

全評価が失敗した場合に限り、例外になる。メッセージは次の 2 つを区別する。

| メッセージ | 状況 | 対処 |
|---|---|---|
| `fitting cannot start` | 初期値で計算できない | 初期値を見直す。記録を間引いて増分が大きくなった場合に生じやすい |
| `no point the optimiser tried was` | 初期値では計算できたが、探索した点がすべて失敗した | 範囲を狭める |

## 関連ページ

- 試験の記録からの同定の手順: [試験データに合わせる](tutorial.md)
- 同定したモデルの比較: [モデルを比べる](comparing.md)
- 非収束の仕組み: [非収束の扱い](nonconvergence.md)
