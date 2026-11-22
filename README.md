# returnmap

<!-- docs:start -->
弾塑性の構成則を Python で記述して計算し、Abaqus の UMAT (Fortran) までを検証しながら進めるためのフレームワークである。降伏関数と硬化則を書けば、リターンマッピングと整合接線 (consistent tangent) は自動微分で解かれる。用途は 2 つある。

- **既存のモデルを使う** — 組み込みのモデルで応力-ひずみ応答を計算し、試験データに対してパラメータを同定し、モデルを比較する。
- **自分のモデルを UMAT にする** — Python のモデルを基準として、手で導いたリターンマッピングと Fortran の UMAT を、段階ごとに照合する。

主な機能:

- **整合接線の自動計算** — autograd による自動微分で求める。接線を手で導出する必要が無い。
- **汎用の応力更新** — 未知数は `[σ, Δλ, implicit の状態]` である。閉形式で更新できる状態は未知数に含まれない。Broyden 更新 (`nr_jacobian="broyden"`) で計算時間を短縮できる。
- **次元** — 3D、平面ひずみ、平面応力、1D を、共通の記述で扱う。
- **手導出の照合** — 手で導出したヤコビアンのブロックをモデル上で宣言し、自動微分とブロックごとに照合する。
- **組み込みモデル** — J2 等方硬化、J2 + Lemaitre 損傷、Armstrong-Frederick、Chaboche、Ohno-Wang、Yoshida-Uemori。
- **サービス向けの契約層** — Web API などへの組み込み用に、dict と配列で入出力する `ModelRegistry`、`run_simulation` を提供する。

## インストール

**使う** — PyPI には未公開なので、GitHub のタグから入れる (依存は numpy / scipy / autograd):

```bash
pip install "git+https://github.com/toruinaba/returnmap@v0.4.0"
# 図も描くなら matplotlib を足す:
pip install "returnmap[examples] @ git+https://github.com/toruinaba/returnmap@v0.4.0"
```

**開発する** — リポジトリを取ってきて、テストや Fortran のビルドまで含めた環境を作る:

```bash
git clone https://github.com/toruinaba/returnmap && cd returnmap
uv sync --extra dev          # テスト依存含む
uv sync --all-extras         # matplotlib / ドキュメント / Fortran ビルド依存含む

# Fortran モジュールのビルド (gfortran 必須、無ければ Docker: make docker-build)
uv run returnmap build --abaqus-stubs src/returnmap/models/j2_isotropic/fortran/j2_isotropic_3d.f90 --name j2_isotropic_3d
```

## クイックスタート

単軸の引張・圧縮の繰返しを、組み込みのモデルで計算する (`1D` のモデルは単軸応力である)。

```python
from returnmap.models import AFKinematic1D
from returnmap.simulation import FieldHistory, PythonIntegrator, StrainDriver

model = AFKinematic1D(E=206_000.0, nu=0.3, sigma_y0=250.0, C_k=20_000.0, gamma=100.0)
load = FieldHistory.cyclic_strain([0.015, -0.015, 0.015], n_per_segment=50)   # ±1.5 % の繰返し

result = StrainDriver(PythonIntegrator(model)).run(load)

print(f"ステップ数           : {len(result.stress)}")
print(f"最初の折返しでの応力 : {result.stress[49, 0]:.1f} MPa")
print(f"最後の累積塑性ひずみ : {result.history('ep')[-1]:.4f}")
```

- 各パラメータの説明、単位、目安の値は `AFKinematic1D.param_meta` にある。
- モデルの選択、載荷履歴、結果の参照は [応答を計算する](https://toruinaba.github.io/returnmap/guide/simulation.html)、試験データに対するパラメータの同定は [試験データに合わせる](https://toruinaba.github.io/returnmap/guide/tutorial.html) を参照。

文書の手順のページのコードは、`examples/` のスクリプトであり、そのまま実行できる。
<!-- docs:end -->

## ドキュメント

**<https://toruinaba.github.io/returnmap/>** で公開している。ガイド (日本語) と API リファレンス (英語、docstring から生成) がある。

**既存のモデルを使う**

| ページ | 内容 |
|---|---|
| [応答を計算する](https://toruinaba.github.io/returnmap/guide/simulation.html) | モデルの選択、載荷履歴の作成、結果の参照、応力制御と混合制御、非収束時の扱い |
| [試験データに合わせる](https://toruinaba.github.io/returnmap/guide/tutorial.html) | 試験の記録に対するパラメータの同定: 候補の比較、フィッティング、全点での仕上げ |
| [モデルを比べる](https://toruinaba.github.io/returnmap/guide/comparing.html) | `compare_models`: 応答、計算コスト、フィットの良さ、パラメータスタディ、解法の比較 |

**自分のモデルを UMAT にする**

| ページ | 内容 |
|---|---|
| [新しいモデルを足す](https://toruinaba.github.io/returnmap/guide/adding-a-model.html) | Python のモデルから検証済みの UMAT までの概要: 段階、リポジトリの構成、完了の条件 |
| [モデルを書く](https://toruinaba.github.io/returnmap/guide/model-python.html) | 段階 A: モデルの定義、検証の宣言、検査の報告 |
| [手書きの解法を照合する](https://toruinaba.github.io/returnmap/guide/model-hand-written.html) | 段階 B: 手書きのリターンマッピングと接線を、自動微分を基準に検査する |
| [UMAT を確かめる](https://toruinaba.github.io/returnmap/guide/model-umat.html) | 段階 C: UMAT の宣言、ビルド、Python のモデルとの比較、不一致の調査 |

**組み込む**

| ページ | 内容 |
|---|---|
| [サービスに組み込む](https://toruinaba.github.io/returnmap/guide/service.html) | `ModelRegistry`、`run_simulation`: Web API 向けの入出力 |

**解説**

| ページ | 内容 |
|---|---|
| [組み込みモデルと次元選択](https://toruinaba.github.io/returnmap/guide/models.html) | モデルの一覧と定式化 (式、パラメータ、性質)、次元の選び方 |
| [フィッティングの設定](https://toruinaba.github.io/returnmap/guide/fitting.html) | `fit_params` の探索空間、複数の試験、比較の方法 (PCM)、driver の指定 |
| [モデルのフック](https://toruinaba.github.io/returnmap/guide/writing-models.html) | フックの仕様: Explicit / Implicit、`flow`、滑らかな近似、弾塑性の判定、弾性予測子、手書きの解法 |
| [検査の一覧](https://toruinaba.github.io/returnmap/guide/verification.html) | 検査の仕様と誤差の測り方: `check_model`、手導出ブロックの照合、接線の差分検査、手書きの解の検査、2 つの実装の比較 |
| [構成と規約](https://toruinaba.github.io/returnmap/concepts.html) | 三層構造、7 つのパッケージ、integrator の種類、Voigt 規約、`FlowVector`、許容値 |
| [非収束の扱い](https://toruinaba.github.io/returnmap/guide/nonconvergence.html) | 収束の判定、解き直し、増分の分割、層ごとの受け取り方 |

**そのほか**

| ページ | 内容 |
|---|---|
| [開発](https://toruinaba.github.io/returnmap/development.html) | テスト、例と文書の作り方、リリース |
| [変更履歴](https://toruinaba.github.io/returnmap/changelog.html) | 版ごとの変更と、0.1.0 からの移行の対応表 |
| [設計メモ](https://toruinaba.github.io/returnmap/design/api-redesign.html) | API 再設計の方針と決定の記録 |

文書のソースは [`docs/`](docs/) と、各モデルの `README.md` (定式化) にある。`main` に入るたびに自動でビルドされ、公開される。手元でビルドするには:

```bash
uv sync --extra docs
make docs            # → docs/_build/html/index.html
```

## License

[MIT License](LICENSE)
