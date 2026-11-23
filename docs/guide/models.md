# 組み込みモデルと次元選択

`returnmap.models` の組み込みモデルと、次元の選び方を示す。モデルを選んで計算する手順は [応答を計算する](simulation.md) の 2 節にある。

## 1. モデルの一覧

各ファミリは、次元を固定しない親クラス (物理) と、次元を固定したクラスを持つ。

| ファミリ | クラス (`*` は `3D` / `PS` / `1D`) | 内容 | 解法 | 定式化 |
|---|---|---|---|---|
| J2 等方硬化 | `J2Isotropic*` | 線形の等方硬化 | `3D` と `1D` は閉形式の解を持つ。`3D` は UMAT を持つ | [式](models/j2_isotropic.md) |
| J2 + Lemaitre 損傷 | `J2IsotropicDamage*` | 延性損傷 (弱連成) | 自動微分 | [式](models/j2_damage.md) |
| 同上 + crack closure | `J2IsotropicDamageCrackClosure*` | 損傷を進める力 `Y` を、主値で引張と圧縮に分ける | 自動微分 | 同上 |
| 同上 + 準片側 | `J2IsotropicDamageUnilateral*` | 応力も分ける。圧縮で剛性が回復する | 自動微分 | 同上 |
| Armstrong-Frederick | `AFKinematic*` | 移動硬化 | 自動微分 | [式](models/af_kinematic.md) |
| Chaboche (背応力 1 つ) | `ChabocheKinematic*` | 移動硬化と、飽和する等方硬化 | 自動微分 | [式](models/chaboche_single.md) |
| Chaboche (背応力 3 つ) | `ChabocheKinematicMulti*` | 背応力と等方硬化を、3 項ずつ重ねる | 自動微分 | [式](models/chaboche_multi.md) |
| Ohno-Wang 型 | `OWKinematic*` | 移動硬化 (implicit の発展則) | `3D` は手書きの解法 (スカラー 1 個の反復) と UMAT を持つ。ほかは自動微分 | [式](models/ow_kinematic.md) |
| Yoshida-Uemori | `YUKinematic*` | 二曲面と停留曲面 | `3D` と `PS` は手書きの解法、手導出のヤコビアン、UMAT を持つ | [式](models/yu_kinematic.md) |

定式化のページは、式、状態変数、パラメータ、飽和値などの性質、検証の範囲を、ファミリごとに示す。本文は、各ファミリのディレクトリにある `README.md` である。

```{toctree}
:hidden:

models/j2_isotropic
models/j2_damage
models/af_kinematic
models/chaboche_single
models/chaboche_multi
models/ow_kinematic
models/yu_kinematic
```

## 2. 次元

| クラス名の末尾 | 次元 | 成分数 | 想定する要素 |
|---|---|---|---|
| `3D` | `SOLID_3D` | 6 | C3D8 などの 3D ソリッド |
| (専用のクラスは無い) | `PLANE_STRAIN` | 4 | CPE4、CAX4 |
| `PS` | `PLANE_STRESS` | 3 | CPS4、シェル |
| `1D` | `UNIAXIAL_1D` | 1 | 1D のトラス、単軸試験のフィッティング |

| 項目 | 規定 |
|---|---|
| 次元を固定したクラス (`*3D`、`*PS`、`*1D`) | 別の `dimension=` を渡すと `TypeError` になる。閉形式の解のように、その次元でのみ成り立つ代数を持つためである |
| 親クラス (`J2Isotropic` など) | 構築時の `dimension=` で次元を選ぶ。省略時は 3D である |
| 平面ひずみ | 専用のクラスが無いため、親クラスに `dimension=PLANE_STRAIN` を渡す。例: `J2Isotropic(dimension=PLANE_STRAIN, E=..., nu=..., sigma_y0=..., H=...)` |
| 成分の順序とせん断の扱い | [構成と規約](../concepts.md) の 4 節 |

## 3. パラメータ

| 項目 | 規定 |
|---|---|
| 構築 | 全パラメータをキーワードで渡す。引数は `param_names` と同じである |
| 既定値 | 無い。UMAT の `PROPS` に既定値が無いことと同じ規則である。機能を無効にする値 (`p_D=0.0` など) も明示する |
| 名前の不足、誤り | その名前を示す `TypeError` になる |
| 意味、単位、目安の値、範囲 | 各クラスの `param_meta` と、定式化のページにある |
| 目安の値での構築 | `Model.at_defaults()` |

## 関連ページ

- モデルを選んで計算する手順: [応答を計算する](simulation.md)
- 自分のモデルを書く手順: [モデルを書く](model-python.md)
