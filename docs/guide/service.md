# サービスに組み込む

Web API やジョブワーカーから returnmap を呼び出す手順を示す。`returnmap.service` は、Web のライブラリに依存しない契約層である。入出力は dict、配列、dataclass のみであり、HTTP やスキーマへの変換は組み込む側が行う。

コードは [`examples/example_service.py`](https://github.com/toruinaba/returnmap/blob/main/examples/example_service.py) の各節であり、上から順に実行できる。

## 1. レジストリ

`ModelRegistry` は、公開する名前とモデルの対応を持つ。全体で共有するインスタンスは無い。組み込む側が作り、公開するモデルだけを登録する。`default_registry()` は、組み込みの 1D モデルを登録したレジストリを返す。

```{literalinclude} ../../examples/example_service.py
:language: python
:start-after: "# %% 1. レジストリ"
:end-before: "# %%"
```

```text
['af_kinematic', 'j2_isotropic', 'ow_kinematic']
J2Isotropic1D 1 ('E', 'nu', 'sigma_y0', 'H')
E         MPa   210000.0  (1000.0, 1000000.0)  Young's modulus
nu                   0.3  (0.0, 0.49)  Poisson's ratio
sigma_y0  MPa      250.0  (1.0, 5000.0)  Initial yield stress
H         MPa     1000.0  (0.0, 100000.0)  Isotropic hardening modulus (linear)
```

| 名前 | 内容 |
|---|---|
| `registry.names()` | 登録されている名前 |
| `registry.get(名前)` | 登録されているモデル (テンプレート) |
| `registry.describe(名前)` | モデルの記述 (`ModelDescriptor`)。`class_name`、`param_names`、`param_meta`、`ntens`、`doc` を持つ。`GET /models` に相当する情報であり、`param_meta` (説明、単位、目安の値、範囲) から入力フォームを生成できる |

`param_meta` は、モデルのクラスに宣言する ([モデルのフック](writing-models.md) の 2 節)。リポジトリの間で情報が食い違わないようにするためである。

## 2. 計算の実行

`run_simulation` は、モデル、パラメータ、載荷履歴を受け取り、配列を返す。

```{literalinclude} ../../examples/example_service.py
:language: python
:start-after: "# %% 2. 計算の実行"
:end-before: "# %%"
```

```text
True 200 (200, 1) (200,)
最終の応力 253.8 MPa、進捗の通知 200 回 (最後は (200, 200))
```

| 引数 | 内容 |
|---|---|
| 第 1 引数 (モデル) | テンプレート (4 節) |
| `params` | パラメータ。`param_names` のすべてを、過不足なく与える |
| `load` | 載荷履歴の配列 (累積値) |
| `control` | `"strain"` (既定) または `"stress"` |
| `include_fields` | 返す状態量の名前。既定では応力とひずみのみを返す |
| `nr_jacobian` | 既定は `"exact"`。速度を優先する場合は `"broyden"` |
| `progress` | ステップごとに呼ばれる関数。引数は (完了したステップ数, 全ステップ数) |

| 結果 (`SimulationOutput`) | 内容 |
|---|---|
| `stress`、`strain` | 応力、ひずみ。形は `(ステップ数, 成分数)` |
| `fields[名前]` | `include_fields` で指定した状態量 |
| `converged`、`failed_step`、`n_steps` | 収束、収束しなかったステップの番号、完了したステップ数 |

## 3. 非収束と入力の誤り

パラメータを探索する利用者は、解の無い値も入力する。`run_simulation` は、途中で収束しなくなっても例外を送出せず、そこまでの履歴を返す。

```{literalinclude} ../../examples/example_service.py
:language: python
:start-after: "# %% 3. 非収束と入力の誤り"
:end-before: "# %%"
```

```text
False 1 1 (1, 1)
ValueError — J2Isotropic1D: parameter mismatch — missing ['H', 'nu', 'sigma_y0']
```

| 状況 | 動作 | API 層での扱いの例 |
|---|---|---|
| 計算が途中で収束しない | 値で返る。`converged=False`、`failed_step`、完了したステップまでの配列 | 200 で返し、途中までの曲線と非収束の位置を表示する |
| パラメータ名の過不足、未知の `include_fields` | `ValueError`。呼び出し側の誤りである | 400 で返す |

## 4. テンプレートの登録

レジストリに登録するものは、クラスではなく、モデルのインスタンス (テンプレート) である。`run_simulation` は、リクエストのパラメータを `model.with_params(**params)` でテンプレートに入れて計算する (`fit_params` と同じ方法)。テンプレート自体は変更されない。

クラスだけでは決まらない次の 2 つを、テンプレートが決める。

| 項目 | 内容 |
|---|---|
| 次元 | 次元を固定していないクラス (`J2Isotropic` など) を、どの次元で公開するか |
| パラメータ以外の設定 | `param_names` に無い設定 (`J2IsotropicDamageUnilateral` の `c` など) |

```{literalinclude} ../../examples/example_service.py
:language: python
:start-after: "# %% 4. テンプレートの登録"
```

```text
['j2_3d', 'j2_ps'] 6 3
```

- **テンプレートが持つパラメータの値は、計算に用いられない。** `run_simulation` の `params` は全パラメータを与える必要があり、テンプレートの値で補われることは無い。リクエストが指定していない値で計算が行われないようにするためである。テンプレートの値は、モデルを構築するための仮の値である。
- `default_registry()` のテンプレートは、各クラスの `param_meta` の目安の値で作られている。
- `ntens` は、`load` の 1 行の幅である。

### 非公開のモデルの登録

非公開のモデルのパッケージは、名前 → テンプレートの dict (`MODELS: dict[str, MaterialModel]`) を公開する。組み込む側は、それを import し、`registry.register(名前, モデル)` で自身のレジストリに登録する。

## 関連ページ

- 計算の基本: [応答を計算する](simulation.md)
- 非収束の仕組み: [非収束の扱い](nonconvergence.md)
