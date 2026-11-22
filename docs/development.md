# 開発

## テスト

```bash
make test          # 高速テスト (フレームワークの unit + integration と全モデル、slow/fortran 除く)
make test-models   # モデルのテストだけ (slow を含む)
make test-all      # 全スイート
uv run pytest src/returnmap/models/yu_kinematic -m "not fortran"    # 1 ファミリだけ
```

その他のターゲット (e2e / slow / docker) は `Makefile` を参照。テストは 2 箇所にある。

- **フレームワークのテスト** — `tests/{unit,integration}` が `src/returnmap` のモジュール構成に対応する (モジュールごとに unit か integration のどちらか)。複数のモデルにまたがるテストもここ (`tests/integration/models/`)。
- **モデルのテスト** — そのモデルのファミリの中、`src/returnmap/models/<family>/tests/` にある。pip で入れた先でも `pytest --pyargs returnmap.models.yu_kinematic` で実行できる。

### モデルのテスト

新しいファミリを足す手順と、1 つのディレクトリに入れる範囲 (1 つのモデルと、それを継承した変種) は [新しいモデルを足す](guide/adding-a-model.md) にある。ファミリのディレクトリは、そのモデルに属するものを全部持つ: モデル本体、検証の宣言 (`verification.py`)、UMAT のソース (`fortran/`)、テスト (`tests/`)。テストは公開の import パス (`returnmap.core` など) とそのファミリの中のものだけを import する — ディレクトリごと動かせるようにするためである。

| 名前 | 役割 |
|---|---|
| `README.md` | 定式化: 式、状態変数、パラメータ、性質、検証の範囲。[組み込みモデルと次元選択](guide/models.md) から読める |
| `verification.py` | 検証の宣言。`MODELS` (検査にかけるモデル)、必要なら `histories(model)` (標準の履歴に足す載荷履歴)、UMAT があれば `UMATS` (ソース、リンクするライブラリ) |
| `tests/test_verification.py` | 共通の検査 `check_model` を `MODELS` に回す。どのファミリも同じ数行で、**この名前は共通の検査専用** |
| `tests/test_model.py` | モデルの関数ごとの確認 (integrator を使わない) |
| `tests/test_<確かめること>.py` | そのモデルに固有の確認 |

「このモデルは検証されているか」は `test_verification.py` を見れば分かる。何が検査されるかは、モデルが宣言しているものから決まる ([検査の一覧](guide/verification.md) の 1 節)。共通の検査は 1 モデルあたり 0.1〜4 秒で、`make test` に入る (YU の 3D と平面応力は 5〜10 秒かかるので `slow`)。

2 つの独立した実装を比べるテスト (手書きの経路と自動微分、Python と UMAT) は、次の名前と基準値に揃える。外れるものは、そのモジュールの docstring に理由を書く。

| ファイル | 比べるもの |
|---|---|
| `test_analytical_vs_numerical.py` | 手書きのリターンマッピング vs フレームワークの NR |
| `test_numerical_vs_fortran.py` | Python vs コンパイルした UMAT (`@pytest.mark.fortran`。モジュールが無ければ skip) |

| 量 | 許容値 |
|---|---|
| 応力 | atol=1e-6 (絶対値、MPa) |
| Δλ、ep | atol=1e-10 (絶対値) |
| 接線 | 最大相対誤差 < 1e-5 |

## 例

`examples/` は、文書のページが取り込むコードのみを持つ。各スクリプトは、対応するページのコードそのものであり、上から順に実行できる。

| スクリプト | 対応するページ | 実行時間 |
|---|---|---|
| `examples/example_response.py` | [応答を計算する](guide/simulation.md) | 約 7 秒 |
| `examples/example_identification.py`、`examples/data/` | [試験データに合わせる](guide/tutorial.md) | 約 50 秒 |
| `examples/notebook_tuning.py` | 同上 (スライダーで調整する marimo のノートブック。`uv sync --extra notebook`) | — |
| `examples/example_compare_models.py` | [モデルを比べる](guide/comparing.md) | 約 15 秒 |
| `examples/model_family/` | [新しいモデルを足す](guide/adding-a-model.md) と、それに続く 3 ページ。モデルのリポジトリの雛形である。`report.py` は検査の報告を表示し、`umat_crosscheck.py` は UMAT との比較を個別に実行する (要 gfortran) | 数秒 |
| `examples/example_service.py` | [サービスに組み込む](guide/service.md) | 1 秒未満 |
| `examples/example_fitting.py` | [フィッティングの設定](guide/fitting.md) | 約 1 秒 |
| `examples/example_checks.py` | [検査の一覧](guide/verification.md) | 約 2 秒 |

## ドキュメント

ガイドは `docs/` の Markdown (日本語)、API リファレンスは docstring (英語) から Sphinx + MyST で生成する。

```bash
uv sync --extra docs
make docs            # docs/_build/html/index.html を生成
```

公開版は <https://toruinaba.github.io/returnmap/> にある。`main` に入るたびに、GitHub Actions (`.github/workflows/docs.yml`) が同じ `make docs` でビルドして公開する。プルリクエストではビルドだけを行い、公開はしない。

`make docs` は警告をエラーとして扱う。次のどれかがあるとビルドが失敗する:

- 解決できない相互参照 (改名・削除された名前へのリンク)
- `docs/` 内のページへのリンク切れ
- 公開パッケージの `__all__` に載っているのに、リファレンスに項目が無い名前 (docstring の無いクラスや関数、`docs/reference/` に説明の無い定数)

docstring の中で他の名前を参照するときは、公開パッケージのパスで書く (`` :class:`~returnmap.core.MaterialModel` ``)。

### ページの作り方

| 項目 | 規約 |
|---|---|
| 文体 | 技術文書の文体で書く。節の最初に範囲を 1〜2 文で述べ、コード、出力、仕様 (表) の順に並べる。見出しは名詞にする。設計の理由は、設計メモと docstring に書く |
| コード | ページに独自の Python のブロックを置かない。`literalinclude` で、実行されるファイルから取り込む。手順のページは `examples/` のスクリプトの節 (`# %%` の区切り)、仕様のページは組み込みモデルのソース (`src/returnmap/models/`) から取り込む |
| 出力 | スクリプトが表示するものを載せる。毎回同じ出力は `text` のブロック、実行ごとに変わる数値 (実行時間、フィッティングの結果、丸め誤差の大きさの値) を含む出力は `none` のブロックに書く |
| 図 | `docs/figures/` のスクリプトが `docs/guide/models/img/` に描く。`make docs-figures` で描き直す (約 1 分) |
| 例の追加 | `examples/` のスクリプトは、どれかのページから参照する。ページとスクリプトの対応は `tests/e2e/test_docs_examples.py` の `PAGES` に書く |

`tests/e2e/test_docs_examples.py` は、スクリプトを実行し、ページに載せた出力と比較する。

| ブロック | 比較の方法 |
|---|---|
| `text` | 各行が、スクリプトの出力の行と 1 文字ずつ一致する |
| `none` | 数値を伏せ、行、列、文言が一致する |

```bash
make test-e2e                 # ページのテストを含む (チュートリアルと図は slow、UMAT のページは fortran)
uv run pytest tests/e2e/test_docs_examples.py     # ページのテストのみ (slow を含む。約 2 分)
```
