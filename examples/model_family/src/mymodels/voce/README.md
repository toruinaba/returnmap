# J2 + Voce 硬化

von Mises の降伏条件に、飽和する等方硬化 (Voce 則) を組み合わせたモデルである。

| クラス | 次元 | リターンマッピング | 接線 | UMAT |
|---|---|---|---|---|
| `J2Voce` | 構築時に `dimension=` で選ぶ | 自動微分 | 自動微分 | — |
| `J2Voce1D` | 単軸応力 | 自動微分 | 自動微分 | — |
| `J2Voce3D` | 3D ソリッド | 手書き | 手書き | `fortran/voce_3d.f90` |

`J2Voce1D` と `J2Voce3D` は、`J2Voce` の式をそのまま継承する。

## 式

### 出典

| 内容 | 出典 |
|---|---|
| 硬化則 | Voce, E. (1948). The relationship between stress and strain for homogeneous deformation. *Journal of the Institute of Metals*, 74, 537–562. |
| 後退 Euler による積分、半径方向のリターン、整合接線 | Simo, J. C., Hughes, T. J. R. (1998). *Computational Inelasticity*. Springer. 第 3 章 |

出典との差は無い。

### 記号と量の種別

| 記号 | 意味 | 種別 |
|---|---|---|
| $\boldsymbol\sigma$ | Cauchy 応力 | 応力系 |
| $\boldsymbol{s}$ | 偏差応力 | 応力系 |
| $\boldsymbol\varepsilon$、$\Delta\boldsymbol\varepsilon$ | 全ひずみ、その増分 | ひずみ系 |
| $\boldsymbol\varepsilon^p$ | 塑性ひずみ | ひずみ系 |
| $p$ | 相当塑性ひずみ | スカラー |
| $\mathbf{C}$ | 等方線形弾性の剛性 ($E$、$\nu$) | — |
| $\mu$ | せん断弾性係数 $E / (2(1+\nu))$ | — |

応力系の量のせん断成分はテンソル成分、ひずみ系の量のせん断成分は工学せん断ひずみ ($2\varepsilon_{12}$) である。成分の順序は $[11, 22, 33, 12, 13, 23]$ である。式中の $:$ はテンソルの二重縮約である。

### 構成式

微小ひずみ、速度非依存である。

$$
\boldsymbol\sigma = \mathbf{C} : (\boldsymbol\varepsilon - \boldsymbol\varepsilon^p)
$$

$$
f = \sigma_\mathrm{vm}(\boldsymbol\sigma) - \sigma_y(p),
\qquad
\sigma_\mathrm{vm} = \sqrt{\tfrac{3}{2}\,\boldsymbol{s}:\boldsymbol{s}},
\qquad
\sigma_y(p) = \sigma_{y0} + Q\left(1 - e^{-b p}\right)
$$

$$
\mathrm{d}\boldsymbol\varepsilon^p = \mathrm{d}\lambda\,\boldsymbol{n},
\qquad
\boldsymbol{n} = \frac{\partial f}{\partial\boldsymbol\sigma} = \frac{3}{2}\frac{\boldsymbol{s}}{\sigma_\mathrm{vm}},
\qquad
\mathrm{d}p = \mathrm{d}\lambda
$$

$$
f \le 0,\qquad \mathrm{d}\lambda \ge 0,\qquad \mathrm{d}\lambda\,f = 0
$$

降伏応力は $\sigma_{y0} + Q$ に飽和し、降伏直後の硬化係数は $Q b$ である。

### 時間積分

後退 Euler である。試行応力を $\boldsymbol\sigma^\mathrm{tr} = \boldsymbol\sigma_n + \mathbf{C}\,\Delta\boldsymbol\varepsilon$ とする。$f(\boldsymbol\sigma^\mathrm{tr}, p_n) \le 0$ のステップは弾性で、$\boldsymbol\sigma_{n+1} = \boldsymbol\sigma^\mathrm{tr}$、$p_{n+1} = p_n$ である。それ以外のステップは次を満たす。

$$
\boldsymbol\sigma_{n+1} = \boldsymbol\sigma^\mathrm{tr} - \Delta\lambda\,\mathbf{C}\,\boldsymbol{n}_{n+1},
\qquad
p_{n+1} = p_n + \Delta\lambda,
\qquad
f_{n+1} = 0
$$

## 状態変数

| 名前 | 種別 | 形 | 初期値 | 意味 |
|---|---|---|---|---|
| `stress` | — | テンソル (応力系) | 0 | Cauchy 応力 $\boldsymbol\sigma$ |
| `ep` | Explicit | スカラー | 0 | 相当塑性ひずみ $p$ |

自動微分の NR の未知数は $[\boldsymbol\sigma, \Delta\lambda]$ である。

## パラメータ

| 名前 | 意味 | 単位 | 目安の値 | 範囲 |
|---|---|---|---|---|
| `E` | ヤング率 | MPa | 210000 | 1000 〜 1000000 |
| `nu` | ポアソン比 | — | 0.3 | 0 〜 0.49 |
| `sigma_y0` | 初期降伏応力 | MPa | 250 | 1 〜 5000 |
| `Q` | 硬化の飽和値 | MPa | 150 | 0 〜 5000 |
| `b` | 飽和の速さ | — | 20 | 0 〜 1000 |

パラメータの宣言の順序は、この表の順序である。

## 性質

モデルが満たす性質である。「式」の節から決まり、実装によらない。

誤差は、性質が指すステップの値を並べた配列の全体で、`max|計算値 − 期待値| / max|期待値|` として測る。パラメータは目安の値である。「3D の単軸応力」は、`J2Voce3D` を 11 成分のひずみだけを与える混合制御 (ほかの成分の応力は 0) で載荷したものである。載荷の区間は等分割で、始点を含まない。

| 番号 | 対象 | 載荷 | 性質 | 許容値 |
|---|---|---|---|---|
| P1 | `J2Voce1D`、3D の単軸応力 | 軸ひずみ 0 → 0.01 → −0.01 → 0.01、各区間 20 ステップ | 塑性ステップの全部で $\lvert\sigma_{11}\rvert = \sigma_y(p)$ | 1e-8 |
| P2 | `J2Voce1D` | 軸ひずみ 0 → 1.0、200 ステップ | 最終ステップで $\sigma_{11} = \sigma_{y0} + Q$ | 1e-6 |
| P3 | `J2Voce1D`、3D の単軸応力 | 軸ひずみ 0 → 0.02、40 ステップ | 全ステップで $\varepsilon_{11} - \sigma_{11}/E = p$ | 1e-10 |
| P4 | `J2Voce3D` | 多軸のひずみ履歴 (`standard_histories` の `"multiaxial"`)。ひずみ制御 | 全ステップで $\tfrac{1}{3}\operatorname{tr}\boldsymbol\sigma = K \operatorname{tr}\boldsymbol\varepsilon$、$K = E / (3(1 - 2\nu))$ | 1e-10 |
| P5 | `J2Voce3D` | $2\varepsilon_{12}$ のみ 0 → 0.02、40 ステップ (ほかのひずみ成分は 0) | 塑性ステップの全部で $\lvert\sigma_{12}\rvert = \sigma_y(p)/\sqrt{3}$ | 1e-8 |
| P6 | 同上 | 同上 | 全ステップで $2\varepsilon_{12} - \sigma_{12}/\mu = \sqrt{3}\,p$ | 1e-10 |
| P7 | `J2Voce3D`、`J2Voce1D` | `standard_histories` の全部。ひずみ制御 | $Q = 0$ としたモデルが、`returnmap.models` の `J2Isotropic3D`、`J2Isotropic1D` ($H = 0$、ほかのパラメータは同じ、既定の経路) と、応力・`ep`・接線で一致する | 1e-8 |
| P8 | `J2Voce3D`、`J2Voce1D` | P1 の載荷 | 3D の単軸応力の $\sigma_{11}$ と $p$ が、`J2Voce1D` のものと一致する | 1e-8 |

`J2Voce3D` を対象とする性質は、手書きの経路と自動微分の経路の両方で成り立つ。

## 手書きの解法 (3D)

### 縮約

3D ソリッドでは $\mathbf{C}\,\boldsymbol{n} = 3\mu\,\boldsymbol{s}^\mathrm{tr}/q$ ($q = \sigma_\mathrm{vm}(\boldsymbol\sigma^\mathrm{tr})$) が成り立ち、リターンは試行偏差応力の向きに沿う。未知数は $\Delta\lambda$ だけになる。

$$
g(\Delta\lambda) = q - 3\mu\,\Delta\lambda - \sigma_y(p_n + \Delta\lambda) = 0
$$

$g$ は下に凸で単調に減り、塑性ステップでは $g(0) > 0$ である。$\Delta\lambda = 0$ からの Newton 法は、根を超えずに収束する。

### 反復

| 項目 | 内容 |
|---|---|
| 方法 | Newton 法 |
| 未知数 | $\Delta\lambda$ の 1 個 |
| 初期値 | $\Delta\lambda = 0$ |
| 更新 | $\Delta\lambda \leftarrow \Delta\lambda + g / (3\mu + H)$、$H = \mathrm{d}\sigma_y/\mathrm{d}p = Q\,b\,e^{-b(p_n + \Delta\lambda)}$ |
| 収束の判定 | $\lvert g\rvert < 10^{-12}\,q$。初期値でも判定する |
| 反復の上限 | 更新 25 回。その後の判定を満たさないステップは非収束である |
| 収束後の確認 | 無い。$g$ の性質から、根は $\Delta\lambda > 0$ である |

1 ステップの更新の回数の実測値は、標準の載荷履歴と「性質」の載荷で 2〜3 回である。

### 応力と接線

$$
\boldsymbol\sigma_{n+1} = \boldsymbol\sigma^\mathrm{tr} - \frac{3\mu\,\Delta\lambda}{q}\,\boldsymbol{s}^\mathrm{tr}
$$

$$
\mathbf{D} = \mathbf{I}_\mathrm{vol}\mathbf{C} + \theta\,\mathbf{I}_\mathrm{dev}\mathbf{C} - \beta\,\boldsymbol{s}^\mathrm{tr}\otimes\boldsymbol{s}^\mathrm{tr},
\qquad
\theta = 1 - \frac{3\mu\,\Delta\lambda}{q},
\qquad
\beta = \frac{9\mu^2\,(\sigma_y - H\,\Delta\lambda)}{(3\mu + H)\,q^3}
$$

$\sigma_y$ と $H$ は、ステップ終了時の $p_{n+1}$ での値である。接線は対称である。

## 近似

無い。手書きの解と接線は厳密である。検査の許容値と載荷履歴は、既定のものである。

## UMAT

`J2Voce3D` の手書きの解法の移植である。

| 項目 | 内容 |
|---|---|
| ソース | `fortran/voce_3d.f90` |
| f2py のモジュール名 | `voce_3d` |
| リンクするライブラリ | 無い |
| 中核のサブルーチン | `voce_3d`。1 ステップの更新を行う |
| 中核の入力 | パラメータ (表の順)、応力、`ep`、ひずみ増分 |
| 中核の出力 | 応力、`ep`、接線、収束フラグ (整数、1 = 収束、0 = 非収束) |
| ABAQUS の入口 | `umat`。`PROPS` と `STATEV` を展開して中核を呼ぶ |
| `PROPS` | 1 = `E`、2 = `nu`、3 = `sigma_y0`、4 = `Q`、5 = `b` |
| `STATEV` | 1 = `ep` |
| 反復と収束後の確認 | 手書きの解法と同じ |
| 非収束時の出力 | 中核と Python の手書きの経路は、最後の反復値を返す (収束フラグは 0) |
| 非収束時の入口 | `PNEWDT = 0.5` とし、`STRESS` と `STATEV` を変更しない |
| 収束時の入口 | `PNEWDT` を変更しない |
| 接線の対称性 | 対称 |
| 精度 | 倍精度。実数の定数は倍精度で書く |

## 検証

### 対象

| クラス | パラメータ |
|---|---|
| `J2Voce3D` | 目安の値 |
| `J2Voce1D` | 目安の値 |

### 報告

共通の検査の報告である。UMAT をビルドした環境 (コンテナ) で、`returnmap check --family mymodels --require-built --range` が表示したものである。

```none
J2Voce3D: every check that applies passes
  ok   converges                 59 steps over 2 histories, 47 plastic
  ok   tangent                   4.6e-10 (tol 1e-05) at 8 plastic steps
  —    hand blocks               no @hand_block on the model
  ok   hand-written solution     4.9e-13 (tol 1e-08) over 3 rows at 47 plastic steps
  ok   hand-written vs autodiff  stress 3.3e-12, state 4.2e-12 (tol 1e-06), tangent 4.8e-12 (tol 1e-05) over 2 histories
  ok   hand-written tangent      4.6e-10 (tol 1e-05) at 8 plastic steps
  ok   UMAT                      stress 2.7e-16, state 5.2e-16 (tol 1e-06), tangent 4.0e-16 (tol 1e-05) over 2 histories (voce_3d)
  ok   UMAT entry                0.0e+00 (tol 1e-06) against voce_3d over 2 histories (voce_3d.umat)
  histories: uniaxial (24 plastic steps), multiaxial (23 plastic steps)
J2Voce3D over its range: every case passes
  parameters — check_model at each end of every documented range
    ok   E=1000              7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   E=1e+06             7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   nu=0                7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   nu=0.49             7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   sigma_y0=1          7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   sigma_y0=5000       7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   Q=0                 7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   Q=5000              7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   b=0                 7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   b=1000              7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
  loading — random strain histories
    ok   random              36 of 36 increments solved, 35 plastic; solution 6.3e-13, hand-written vs autodiff 5.5e-12, tangent 2.4e-10, UMAT 4.2e-16, UMAT entry 0.0e+00
    ok   random large        24 of 24 increments solved, 24 plastic; solution 1.8e-13, hand-written vs autodiff 3.2e-11, tangent 1.1e-10, UMAT 8.1e-16, UMAT entry 0.0e+00
    ok   random single-step  6 of 6 increments solved, 6 plastic; solution 1.1e-14, hand-written vs autodiff 2.9e-12, tangent 6.0e-10, UMAT 2.4e-16, UMAT entry 0.0e+00

J2Voce1D: every check that applies passes
  ok   converges                 35 steps over 1 histories, 24 plastic
  ok   tangent                   7.2e-10 (tol 1e-05) at 4 plastic steps
  —    hand blocks               no @hand_block on the model
  —    hand-written solution     no user_defined_return_mapping on the model
  —    hand-written vs autodiff  no user_defined_return_mapping on the model
  —    hand-written tangent      no user_defined_return_mapping on the model
  —    UMAT                      no UMAT declared for J2Voce1D in mymodels.voce.verification
  —    UMAT entry                no UMAT declared for J2Voce1D
  histories: uniaxial (24 plastic steps)
J2Voce1D over its range: every case passes
  parameters — check_model at each end of every documented range
    ok   E=1000              2 checks pass (converges, tangent)
    ok   E=1e+06             2 checks pass (converges, tangent)
    ok   nu=0                2 checks pass (converges, tangent)
    ok   nu=0.49             2 checks pass (converges, tangent)
    ok   sigma_y0=1          2 checks pass (converges, tangent)
    ok   sigma_y0=5000       2 checks pass (converges, tangent)
    ok   Q=0                 2 checks pass (converges, tangent)
    ok   Q=5000              2 checks pass (converges, tangent)
    ok   b=0                 2 checks pass (converges, tangent)
    ok   b=1000              2 checks pass (converges, tangent)
  loading — random strain histories
    ok   random              36 of 36 increments solved, 16 plastic; tangent 3.1e-10
    ok   random large        24 of 24 increments solved, 23 plastic; tangent 2.6e-09
    ok   random single-step  6 of 6 increments solved, 6 plastic; tangent 3.4e-06

2 model(s) checked: every one passes
```

「性質」の P1〜P8 のテスト (`tests/test_properties.py`、18 件) は、全部合格である。

### 範囲

| 項目 | 状態 |
|---|---|
| 自動微分の経路の収束と接線、手書きの解と接線、UMAT の中核、ABAQUS の入口 | `tests/test_verification.py` (標準の載荷履歴、目安の値)。UMAT の 2 行は、ビルドした環境でだけ実行される |
| 性質 P1〜P8 | `tests/test_properties.py`。期待値は、テストの中で式から計算する |
| パラメータの範囲の端、乱数の載荷履歴 (降伏ひずみの 100 倍の増分まで) | テストには無い。上の報告のとおり、`returnmap check --range` で 1 回確かめた |
| 反復の上限に達するステップ | 確かめていない。標準の載荷履歴、性質の載荷、乱数の載荷履歴では生じない |
| ABAQUS での実行 | 確かめていない |
