# UMAT を確かめる

手書きの解法を Fortran の UMAT に移植し、Python のモデルと全ステップで比較する手順を示す ([新しいモデルを足す](adding-a-model.md) の段階 C)。前提は [手書きの解法を照合する](model-hand-written.md) である。

コードは [`examples/model_family/`](https://github.com/toruinaba/returnmap/tree/main/examples/model_family) のものである。ビルドには gfortran が要る。gfortran の無いホストでは Docker を用いる (末尾の「Docker」)。

## 1. Fortran のソース

ソースはファミリの `fortran/` に置き、2 つのサブルーチンに分ける。

| サブルーチン | 役割 |
|---|---|
| 中核 (`voce_3d`) | 1 ステップの構成則の更新。f2py でビルドし、Python のモデルと比較する対象である |
| 入口 (`umat`) | ABAQUS が呼ぶサブルーチン。`PROPS` と `STATEV` を展開して中核を呼ぶ。中核が収束しなかった場合は、`PNEWDT` で増分の縮小を要求する |

```{literalinclude} ../../examples/model_family/src/mymodels/voce/fortran/voce_3d.f90
:language: fortran
:start-at: "subroutine voce_3d(E, nu, sigma_y0, Q, b, &"
:end-at: "intent(out) :: converged"
```

中核の引数は、次の順に並べる。この順であれば、呼び出し方を宣言に書く必要が無い。

| 区分 | 順序 |
|---|---|
| 入力 | パラメータ (`param_names` の順)、応力、状態変数 (宣言の順)、ひずみ増分 |
| 出力 | 応力、状態変数 (宣言の順)、接線。その後に、収束フラグなどの追加の出力 |

## 2. UMAT の宣言

UMAT を、ファミリの `verification.py` の `UMATS` に宣言する。ビルドと検査は、この宣言を参照する。

```{literalinclude} ../../examples/model_family/src/mymodels/voce/verification.py
:language: python
:start-at: "def _umat_failure(ret):"
```

| `Umat` の引数 | 内容 |
|---|---|
| 第 1 引数 | UMAT が移植したモデルのクラス。UMAT は 1 つの次元に対して書くため、次元を固定したクラスを指定する |
| `module` | ビルドする f2py モジュールの名前 |
| `sources` | Fortran のソース。モデルのモジュールがあるディレクトリからの相対パスで指定する |
| `libraries` | リンクするライブラリ (`["lapack", "blas"]` など) |
| `subroutine` | 1 ステップの更新を行うサブルーチンの名前。省略時は `module` と同じ名前 |
| `props` | UMAT が受け取るパラメータの名前と順序。`param_names` と順序が異なる場合、または一部のみを受け取る場合に指定する |
| `integrator_options` | `FortranIntegrator.from_model` に渡す追加の引数。主な用途は `parse_umat_failure` である |
| `entry` | ABAQUS が呼ぶ入口のサブルーチンの名前。既定は `umat`。入口を持たないモジュール、または検査しない場合は `None` |
| `statev` | 入口が `STATEV` に状態変数を並べる順序。テンソルは成分の数だけ場所を占める。省略時は宣言の順 (応力を除く) |

**UMAT が収束フラグを返す場合は、`parse_umat_failure` で読ませる。** 指定しない場合、UMAT が収束しなかったステップの値が、解として比較される。`parse_umat_failure` はサブルーチンの戻り値の組を受け取り、失敗の説明 (文字列) か、成功を表す `None` を返す。

## 3. ビルド

`returnmap build --family <パッケージ>` は、そのパッケージのモデルが宣言する UMAT をすべてビルドする。ソース、リンクするライブラリ、ABAQUS のユーティリティルーチンの代替 (スタブ) は、宣言とパッケージから取られる。

```bash
uv sync --extra dev --extra fortran           # meson、ninja を追加する
uv run returnmap build --family mymodels      # 宣言された UMAT をすべてビルドする
uv run returnmap list                         # ビルド済みのモジュールの一覧
uv run returnmap check --family mymodels --require-built   # 検査。未ビルドの UMAT を不合格に数える
```

UMAT を 1 つも宣言していないパッケージに対しては、`returnmap build --family` は何もビルドせずに終了する (終了コード 0)。UMAT を足す前と後で、同じコマンドの組を使える。

ソースを直接指定してビルドする場合は、次の指定を用いる。

```bash
uv run returnmap build --abaqus-stubs src/mymodels/voce/fortran/voce_3d.f90 --name voce_3d
```

| 指定 | 内容 |
|---|---|
| `--name` | モジュールの名前 (ソースを指定する場合は必須) |
| `--abaqus-stubs` | ABAQUS のユーティリティルーチン (`SINV`、`SPRINC`、`ROTSIG`) の代替を一緒にコンパイルする。UMAT のソースを ABAQUS の外でリンクするために要る |
| `--link lapack --link blas` | リンクするライブラリ。繰り返し指定できる |
| `--output-dir` | 出力先。既定はプロジェクト直下の `fortran/` |

スタブはパッケージに含まれており、`returnmap stubs` がそのパスを表示する (f2py やコンパイラを直接呼ぶ場合に用いる)。

## 4. 検査の報告

ビルドした環境では、`check_model` の報告の `UMAT` と `UMAT entry` の行が実行される。

```none
  ok   UMAT                      stress 2.7e-16, state 5.2e-16 (tol 1e-06), tangent 4.0e-16 (tol 1e-05) over 2 histories (voce_3d)
  ok   UMAT entry                0.0e+00 (tol 1e-12) against voce_3d over 2 histories (voce_3d.umat)
```

| 行 | 比較するもの | 方法 |
|---|---|---|
| `UMAT` | 中核のサブルーチンと、Python のモデル | 標準の載荷履歴を、それぞれ独立に進める。全ステップで応力・状態・接線を比較する |
| `UMAT entry` | 入口 (`umat`) と、中核のサブルーチン | 中核が進めた各ステップの開始時の応力・状態と、同じひずみ増分を、入口に ABAQUS の引数の形で与える。`PROPS` はパラメータ、`STATEV` は状態変数である。許容値は `1e-12` で、`check_model` の引数では変わらない (入口は中核の値を受け渡すのみである) |

入口の検査は、次の条件を前提とする。

| 項目 | 内容 |
|---|---|
| 引数の宣言 | 1 節の例と同じ `intent` で宣言する: `STRESS`・`STATEV`・`PNEWDT` は `intent(inout)`、`DDSDDE` とエネルギーの出力は `intent(out)`、ほかは `intent(in)` |
| 非収束の報告 | `PNEWDT` を 1 未満にする。このとき、`STRESS` と `STATEV` を変更しない |
| 入口に与えない量 | 全ひずみ (0)、時刻、温度、変形勾配 (単位行列)。これらを読む UMAT は、この検査の対象外である |
| `STATEV` の内容 | モデルの状態変数のみ。それ以外の場所を持つ UMAT は `entry=None` と宣言し、入口を個別に検査する |
| 入口を持たないモジュール | `entry=None` と宣言する。宣言しない場合、ビルドしたモジュールに `umat` が無いと、この行は不合格になる |

ビルドしていない環境では、2 つの行は `—` (`not built`) になり、不合格には数えない。gfortran の無い環境でも共通の検査が通るようにするためである。UMAT の検査を必ず実行する場所では、`fortran` マーカーのテストで行の状態を確認する。

```{literalinclude} ../../examples/model_family/src/mymodels/voce/tests/test_verification.py
:language: python
:start-at: "@pytest.mark.fortran"
```

| 項目 | 内容 |
|---|---|
| `fortran` マーカー | returnmap が配布する pytest プラグインが登録する。コンパイルしたモジュールが無い場合、テストは skip される |
| モジュールの import | 同じプラグインが、ビルド先 (`fortran/`) を `sys.path` に加える。`conftest.py` への記述は不要である |

## 5. 個別の比較

`check_model` が行う比較は、`crosscheck` で個別に実行できる。`Umat.integrator(model)` は、宣言から UMAT の integrator (`FortranIntegrator`) を作る。

```{literalinclude} ../../examples/model_family/umat_crosscheck.py
:language: python
:start-after: "# %% The UMAT as an integrator"
:end-before: "# %%"
```

```none
crosscheck: 24 steps — all agree
  ok   stress: 2.661e-16 (tol 1e-06)
  ok   state[ep]: 5.217e-16 (tol 1e-06)
  ok   tangent: 3.366e-16 (tol 1e-05)
```

## 6. 不一致の調査

不一致がある場合、`report.first_failure` が最初に一致しなかったステップを返す。このステップは、2 つの実装の結果 (`result_a`、`result_b`) をそのまま持つ。

次の例は、誤りのある UMAT の代わりとして、硬化の速さ `b` を変えた Python のモデルを UMAT と比較し、不一致を作っている。

```{literalinclude} ../../examples/model_family/umat_crosscheck.py
:language: python
:start-after: "# %% A disagreement, and where it comes from"
:end-before: "# %%"
```

```none
crosscheck: 24 of 24 steps disagree (first at step 0)
  FAIL stress: 4.646e-03 (tol 1e-06)
  FAIL state[ep]: 7.178e-03 (tol 1e-06)
  FAIL tangent: 8.560e-03 (tol 1e-05)
step 0, plastic: True
ep        2.1168e-04 (Python)   2.1232e-04 (UMAT)
stress 11    307.393 (Python)      307.314 (UMAT)
jacobian comparison: 4 blocks — 1 differ
  FAIL stress::stress: 1.029e-04 (tol 1e-08)
  ok   stress::dlambda: 2.375e-16 (tol 1e-08)
  ok   dlambda::stress: 1.463e-16 (tol 1e-08)
  ok   dlambda::dlambda: 0.000e+00 (tol 1e-08)
```

| 名前 | 内容 |
|---|---|
| `report.first_failure` | 最初に一致しなかったステップ。全ステップが一致した場合は `None` |
| `step.index` | ステップの番号 (0 始まり) |
| `step.result_a`、`step.result_b` | 2 つの実装の応力更新の結果。応力、状態、`dlambda` を直接比較できる |
| `compare_jacobians(model, result_a, result_b)` | 2 つの解における残差のヤコビアンを、ブロックごとに比較する |

既定では、2 つの実装はそれぞれ自身の結果から次のステップに進む。最初の不一致より後のステップは、すべて不一致になる。調査の対象は最初のステップである。

## 7. driver での実行

UMAT の integrator は、Python の integrator と同じく driver に渡せる。次の例は、3D の UMAT を混合制御で単軸応力にしている ([応答を計算する](simulation.md) の 5 節)。

```{literalinclude} ../../examples/model_family/umat_crosscheck.py
:language: python
:start-after: "# %% The UMAT under a driver"
```

```none
5.7e-14
```

## Docker

gfortran の無いホストでは、Docker のイメージの中でビルドとテストを行う。returnmap の [`Dockerfile`](https://github.com/toruinaba/returnmap/blob/main/Dockerfile) が雛形になる (gfortran と Python を含む)。returnmap 自身では、次のコマンドで実行する。

```bash
make docker-build      # イメージを作る
make docker-test       # イメージの中で UMAT をビルドし、fortran マーカーのテストを実行する
```

## 宣言を用いない場合

`Umat` の宣言を用いず、integrator を直接組み立てることもできる。

| 手順 | 呼び出し |
|---|---|
| モジュールの読み込み | `FortranModule(モジュール名)` |
| integrator の作成 | `FortranIntegrator.from_model(fortran, サブルーチン名, model)`。引数の規約は 1 節と同じである |
| `PROPS` の順序が異なる場合 | `param_fn=lambda: (model.E, model.nu, ...)` |
| 収束フラグの読み取り | `parse_umat_failure=...` (2 節) |

## 関連ページ

- 各検査の測定内容と許容値: [検査の一覧](verification.md)
- ビルド先のディレクトリの詳細: リポジトリの `fortran/README.md`
