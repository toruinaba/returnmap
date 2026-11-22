# returnmap

```{include} ../README.md
:start-after: <!-- docs:start -->
:end-before: <!-- docs:end -->
```

## どこから読むか

| 目的 | ページ |
|---|---|
| 既存のモデルで応答を計算する | [応答を計算する](guide/simulation.md) |
| 試験データに対してパラメータを同定する | [試験データに合わせる](guide/tutorial.md) |
| モデル、パラメータ、解法を比較する | [モデルを比べる](guide/comparing.md) |
| 自分のモデルを書き、UMAT まで進める | [新しいモデルを足す](guide/adding-a-model.md) |
| Web API などに組み込む | [サービスに組み込む](guide/service.md) |
| モデルの式、フックや検査の仕様、規約を調べる | 「解説」の各ページ |

手順のページのコードは、`examples/` のスクリプトであり、そのまま実行できる。

```{toctree}
:caption: 既存のモデルを使う
:maxdepth: 1

guide/simulation
guide/tutorial
guide/comparing
```

```{toctree}
:caption: 自分のモデルを UMAT にする
:maxdepth: 1

guide/adding-a-model
guide/model-python
guide/model-hand-written
guide/model-umat
```

```{toctree}
:caption: 組み込む
:maxdepth: 1

guide/service
```

```{toctree}
:caption: 解説
:maxdepth: 1

guide/models
guide/fitting
guide/writing-models
guide/verification
concepts
guide/nonconvergence
```

```{toctree}
:caption: リファレンス
:maxdepth: 1

reference/index
```

```{toctree}
:caption: 開発
:maxdepth: 1

development
changelog
design/api-redesign
```
