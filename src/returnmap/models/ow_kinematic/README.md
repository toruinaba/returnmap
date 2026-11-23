# Ohno-Wang 型の移動硬化

Armstrong-Frederick の動的回復を、背応力の大きさで重み付けしたモデルである。離散化した発展則は、背応力について閉じた形で解けない。

| クラス | 次元 | リターンマッピング | 接線 | UMAT |
|---|---|---|---|---|
| `OWKinematic` | 構築時に `dimension=` で選ぶ | 自動微分 | 自動微分 | — |
| `OWKinematic1D` | 単軸応力 | 自動微分 | 自動微分 | — |
| `OWKinematic3D` | 3D ソリッド | 手書き | 手書き | `fortran/ow_kinematic_3d.f90` |
| `OWKinematicPS` | 平面応力 | 自動微分 | 自動微分 | — |

`OWKinematic1D`、`OWKinematic3D`、`OWKinematicPS` は、`OWKinematic` の式をそのまま継承する。returnmap では、implicit な状態変数を持つモデルの参照実装である。

## 式

### 出典

| 内容 | 出典 |
|---|---|
| 動的回復を背応力の大きさで制御する移動硬化則 | N. Ohno, J.-D. Wang (1993), *Kinematic hardening rules with critical state of dynamic recovery*, International Journal of Plasticity. |

出典との差は次のとおりである。このモデルは、出典の考え方を背応力 1 つ・べき指数 1 の形に簡略化したものである。複数の背応力の重ね合わせ、臨界状態 $\bar\alpha_i = r_i$ に対するべき乗、塑性ひずみ速度の背応力方向への射影は含まない。

### 記号と量の種別

| 記号 | 意味 | 種別 |
|---|---|---|
| $\boldsymbol\sigma$ | Cauchy 応力 | 応力系 |
| $\boldsymbol{s}$ | 偏差応力 | 応力系 |
| $\boldsymbol\alpha$ | 背応力 (偏差) | 応力系 |
| $\boldsymbol\xi$ | 有効応力 $\boldsymbol{s} - \boldsymbol\alpha$ | 応力系 |
| $\boldsymbol\varepsilon$、$\Delta\boldsymbol\varepsilon$ | 全ひずみ、その増分 | ひずみ系 |
| $\boldsymbol\varepsilon^p$ | 塑性ひずみ | ひずみ系 |
| $p$ | 相当塑性ひずみ | スカラー |
| $\mathbf{C}$ | 等方線形弾性の剛性 ($E$、$\nu$) | — |
| $\mu$ | せん断弾性係数 $E / (2(1+\nu))$ | — |

応力系の量のせん断成分はテンソル成分、ひずみ系の量のせん断成分は工学せん断ひずみ ($2\varepsilon_{12}$) である。成分の順序は、3D ソリッドでは $[11, 22, 33, 12, 13, 23]$、平面応力では $[11, 22, 12]$ である。式中の $:$ はテンソルの二重縮約である。

### 構成式

微小ひずみ、速度非依存である。

$$
\boldsymbol\sigma = \mathbf{C} : (\boldsymbol\varepsilon - \boldsymbol\varepsilon^p)
$$

$$
f = \bar\xi - \sigma_{y0},
\qquad
\boldsymbol\xi = \boldsymbol{s} - \boldsymbol\alpha,
\qquad
\bar\xi = \sqrt{\tfrac{3}{2}\,\boldsymbol\xi:\boldsymbol\xi}
$$

$$
\mathrm{d}\boldsymbol\varepsilon^p = \mathrm{d}\lambda\,\boldsymbol{n},
\qquad
\boldsymbol{n} = \frac{\partial f}{\partial\boldsymbol\sigma} = \frac{3}{2}\frac{\boldsymbol\xi}{\bar\xi},
\qquad
\mathrm{d}p = \mathrm{d}\lambda
$$

$$
\mathrm{d}\boldsymbol\alpha = \tfrac{2}{3}\,C_k\,\boldsymbol{n}\,\mathrm{d}p - \gamma\,\bar\alpha\,\boldsymbol\alpha\,\mathrm{d}p,
\qquad
\bar\alpha = \sqrt{\tfrac{3}{2}\,\boldsymbol\alpha:\boldsymbol\alpha}
$$

$$
f \le 0,\qquad \mathrm{d}\lambda \ge 0,\qquad \mathrm{d}\lambda\,f = 0
$$

### 時間積分

後退 Euler である。試行応力を $\boldsymbol\sigma^\mathrm{tr} = \boldsymbol\sigma_n + \mathbf{C}\,\Delta\boldsymbol\varepsilon$ とする。$f(\boldsymbol\sigma^\mathrm{tr}, \boldsymbol\alpha_n) \le 0$ のステップは弾性で、$\boldsymbol\sigma_{n+1} = \boldsymbol\sigma^\mathrm{tr}$、$\boldsymbol\alpha_{n+1} = \boldsymbol\alpha_n$、$p_{n+1} = p_n$ である。それ以外のステップは次を満たす。

$$
\boldsymbol\sigma_{n+1} = \boldsymbol\sigma^\mathrm{tr} - \Delta\lambda\,\mathbf{C}\,\boldsymbol{n}_{n+1}
$$

$$
\boldsymbol\alpha_{n+1} - \boldsymbol\alpha_n - \tfrac{2}{3}\,C_k\,\Delta\lambda\,\boldsymbol{n}_{n+1} + \gamma\,\Delta\lambda\,\bar\alpha_{n+1}\,\boldsymbol\alpha_{n+1} = \boldsymbol{0}
$$

$$
p_{n+1} = p_n + \Delta\lambda,
\qquad
f_{n+1} = 0
$$

## 状態変数

| 名前 | 種別 | 形 | 初期値 | 意味 |
|---|---|---|---|---|
| `stress` | Implicit | テンソル (応力系) | 0 | Cauchy 応力 $\boldsymbol\sigma$ |
| `alpha` | Implicit | テンソル (応力系) | 0 | 背応力 $\boldsymbol\alpha$ (偏差) |
| `ep` | Implicit | スカラー | 0 | 相当塑性ひずみ $p$ |

自動微分の NR の未知数は $[\boldsymbol\sigma, \Delta\lambda, \boldsymbol\alpha, p]$ である。

単軸応力では、`alpha` は偏差成分 $\alpha_{11}$ を持つ。単軸の背応力 (降伏曲面の中心の応力) は $a = \tfrac{3}{2}\alpha_{11}$ である。

平面応力 (returnmap の次元 `PLANE_STRESS`) では、`stress` は 11・22・12 成分を持ち、33・13・23 成分は 0 である。`alpha` は偏差テンソルの 11・22・12 成分を持ち、33 成分は $-(\alpha_{11} + \alpha_{22})$ である。

## パラメータ

| 名前 | 意味 | 単位 | 目安の値 | 範囲 |
|---|---|---|---|---|
| `E` | ヤング率 | MPa | 210000 | 1000 〜 1000000 |
| `nu` | ポアソン比 | — | 0.3 | 0 〜 0.49 |
| `sigma_y0` | 降伏応力 (一定) | MPa | 250 | 1 〜 5000 |
| `C_k` | 移動硬化係数 | MPa | 20000 | 0 〜 1000000 |
| `gamma` | 動的回復の係数 | 1/MPa | 100 | 0 〜 10000 |

パラメータの宣言の順序は、この表の順序である。$\gamma\,\bar\alpha$ が無次元になるので、$\gamma$ は応力の逆数の単位を持つ。

## 性質

モデルが満たす性質である。「式」の節から決まり、実装によらない。

誤差は、性質が指すステップの値を並べた配列の全体で、`max|計算値 − 期待値| / max|期待値|` として測る。パラメータは、目安の値と、目安の値の $\gamma$ を 1 に替えた値の 2 組である (P9 と P10 を除く)。「3D の単軸応力」は、`OWKinematic3D` を 11 成分のひずみだけを与える混合制御 (ほかの成分の応力は 0) で載荷したものである。単軸の背応力 $a$ は、1D では `alpha` の $\tfrac{3}{2}$ 倍、3D では `alpha` の 11 成分の $\tfrac{3}{2}$ 倍である。

| 番号 | 対象 | 載荷 | 性質 | 許容値 |
|---|---|---|---|---|
| P1 | `OWKinematic1D`、3D の単軸応力 | 軸ひずみ 0 → 0.01 → −0.01 → 0.01、各区間 20 ステップ | 塑性ステップの全部で $\lvert\sigma_{11} - a\rvert = \sigma_{y0}$ | 1e-8 |
| P2 | `OWKinematic1D` | 軸ひずみ 0 → 0.2、200 ステップ | 最終ステップで $\sigma_{11} = \sigma_{y0} + \sqrt{C_k/\gamma}$ | 1e-8 |
| P3 | `OWKinematic1D`、3D の単軸応力 | 軸ひずみ 0 → 0.02、40 ステップ | 全ステップで $\varepsilon_{11} - \sigma_{11}/E = p$ | 1e-10 |
| P4 | 同上 | 同上 | 塑性ステップの全部で $a_{n+1} - a_n = (p_{n+1} - p_n)\,(C_k - \gamma\,a_{n+1}^2)$ | 1e-6 |
| P5 | `OWKinematic3D` | 多軸のひずみ履歴 (`standard_histories` の `"multiaxial"`) | 全ステップで $\tfrac{1}{3}\operatorname{tr}\boldsymbol\sigma = K \operatorname{tr}\boldsymbol\varepsilon$、$K = E / (3(1 - 2\nu))$ | 1e-10 |
| P6 | 同上 | 同上 | 全ステップで $\alpha_{11} + \alpha_{22} + \alpha_{33} = 0$。誤差は `max|α11 + α22 + α33| / max|α の成分|` | 1e-10 |
| P7 | `OWKinematic3D` | $2\varepsilon_{12}$ のみ 0 → 0.02、40 ステップ (ほかのひずみ成分は 0) | 塑性ステップの全部で $\lvert\sigma_{12} - \alpha_{12}\rvert = \sigma_{y0}/\sqrt{3}$ | 1e-8 |
| P8 | 同上 | 同上 | 全ステップで $2\varepsilon_{12} - \sigma_{12}/\mu = \sqrt{3}\,p$ | 1e-10 |
| P9 | `OWKinematic3D`、`OWKinematic1D` | `standard_histories` の全部 | 目安の値の $\gamma$ を 0 に替えたモデルが、`returnmap.models` の `AFKinematic3D`、`AFKinematic1D` (`gamma = 0`、ほかのパラメータは同じ、既定の経路) と、応力・`alpha`・`ep`・接線で一致する | 1e-8 |
| P10 | `OWKinematic1D` | P3 の載荷 | 目安の値の $\gamma$ を 0 に替えたモデルで、塑性ステップの全部の接線が $E\,C_k / (E + C_k)$ | 1e-8 |
| P11 | `OWKinematic3D`、`OWKinematic1D` | P1 の載荷 | 3D の単軸応力の $\sigma_{11}$、$p$、$a$ が、`OWKinematic1D` のものと一致する | 1e-6 |

`OWKinematic3D` を対象とする性質は、手書きの経路と自動微分の経路の両方で成り立つ。

平面応力のクラスの性質は、次のとおりである。パラメータは、上と同じ 2 組である。

| 番号 | 対象 | 載荷 | 性質 | 許容値 |
|---|---|---|---|---|
| P12 | `OWKinematicPS` | `OWKinematicPS` から作った `standard_histories` の全部。ひずみ制御 | 全ステップで、応力と `alpha` (11・22・12 成分)、`ep` が、`OWKinematic3D` (既定の経路) に 11・22・12 成分の同じひずみを与え、33・13・23 成分の応力を 0 とする混合制御で載荷した結果と一致する | 1e-6 |

## 手書きの解法 (3D)

### 縮約

3D ソリッドでは、偏差の方向 $\boldsymbol{n}$ に対して $\mathbf{C}\,\boldsymbol{n} = 2\mu\,\boldsymbol{n}$ が成り立つ。塑性ステップの未知数は、次のスカラー $w$ の 1 個になる。

$$
w = \frac{1}{1 + \gamma\,\Delta\lambda\,\bar\alpha_{n+1}},
\qquad
0 < w \le 1
$$

背応力の式は $\boldsymbol\alpha_{n+1} = w\,(\boldsymbol\alpha_n + \tfrac{2}{3} C_k\,\Delta\lambda\,\boldsymbol{n})$ と書ける。これと $\boldsymbol{s}_{n+1} = \boldsymbol{s}^\mathrm{tr} - 2\mu\,\Delta\lambda\,\boldsymbol{n}$、$\boldsymbol\xi_{n+1} = \tfrac{2}{3}\sigma_{y0}\,\boldsymbol{n}$ ($f_{n+1} = 0$) から $\tfrac{2}{3}\,(\sigma_{y0} + (3\mu + C_k w)\,\Delta\lambda)\,\boldsymbol{n} = \boldsymbol{s}^\mathrm{tr} - w\,\boldsymbol\alpha_n$ が得られ、$\boldsymbol{n}$ と $\Delta\lambda$ は $w$ の閉じた式になる。$\boldsymbol{s}^\mathrm{tr}$ は試行応力の偏差である。

$$
\boldsymbol\eta = \boldsymbol{s}^\mathrm{tr} - w\,\boldsymbol\alpha_n,
\qquad
q = \sqrt{\tfrac{3}{2}\,\boldsymbol\eta:\boldsymbol\eta},
\qquad
\boldsymbol{n} = \frac{3}{2}\frac{\boldsymbol\eta}{q},
\qquad
c = \boldsymbol{n}:\boldsymbol\alpha_n
$$

$$
M = 3\mu + C_k\,w,
\qquad
\Delta\lambda = \frac{q - \sigma_{y0}}{M},
\qquad
\boldsymbol{b} = \boldsymbol\alpha_n + \tfrac{2}{3}\,C_k\,\Delta\lambda\,\boldsymbol{n},
\qquad
A = \sqrt{\tfrac{3}{2}\,\boldsymbol{b}:\boldsymbol{b}}
$$

$\bar\alpha_{n+1} = w\,A$ なので、$w$ の定義が残りの 1 つの式である。

$$
g(w) = w - 1 + \gamma\,\Delta\lambda(w)\,w^2 A(w) = 0
$$

$g$ は無次元である。$w$ についての微分は次のとおりである。

$$
\frac{\mathrm{d}g}{\mathrm{d}w} = 1 + 2\gamma\,w\,\Delta\lambda\,A
- \gamma\,w^2 \left[ \frac{(A + k_1)\,(c + C_k\,\Delta\lambda)}{M} + \frac{k_2\,(\bar\alpha_n^2 - c^2)}{q} \right]
$$

$$
k_1 = \frac{C_k\,\Delta\lambda\,(c + C_k\,\Delta\lambda)}{A},
\qquad
k_2 = \frac{C_k\,\Delta\lambda^2}{A},
\qquad
\bar\alpha_n^2 = \tfrac{3}{2}\,\boldsymbol\alpha_n:\boldsymbol\alpha_n
$$

$A = 0$ のときは $k_1 = k_2 = 0$ とする。この式は、次の微分を $g$ に代入したものである。

$$
\frac{\mathrm{d}q}{\mathrm{d}w} = -c,
\qquad
\frac{\mathrm{d}c}{\mathrm{d}w} = \frac{c^2 - \bar\alpha_n^2}{q},
\qquad
\frac{\mathrm{d}\Delta\lambda}{\mathrm{d}w} = -\frac{c + C_k\,\Delta\lambda}{M},
\qquad
\Delta\lambda\,\frac{\mathrm{d}A}{\mathrm{d}w} = k_1\,\frac{\mathrm{d}\Delta\lambda}{\mathrm{d}w} + k_2\,\frac{\mathrm{d}c}{\mathrm{d}w}
$$

根の位置について、次が成り立つ。

| 条件 | 性質 |
|---|---|
| 塑性ステップ ($f(\boldsymbol\sigma^\mathrm{tr}, \boldsymbol\alpha_n) > 0$) | $q(1) > \sigma_{y0}$ なので $\Delta\lambda(1) > 0$、$g(1) = \gamma\,\Delta\lambda(1)\,A(1) \ge 0$ |
| $q(w) \le \sigma_{y0}$ ($\Delta\lambda(w) \le 0$) となる $w < 1$ | $g(w) \le w - 1 < 0$ |
| $w \to 0$ | $g \to -1$ |
| $\gamma = 0$ | $g(1) = 0$。根は $w = 1$ で、$\boldsymbol\alpha_{n+1} = \boldsymbol{b}$ (線形の移動硬化) |

したがって根は区間 $(0, 1]$ にあり、根では $\Delta\lambda > 0$ である。

### 反復

| 項目 | 内容 |
|---|---|
| 方法 | 区間で保護した Newton 法。Newton 法の更新が区間の外に出るときは、区間の中点を使う |
| 未知数 | $w$ の 1 個 |
| 初期値 | $w = 1$、区間 $(w_\mathrm{lo}, w_\mathrm{hi}) = (0, 1)$ |
| 更新 (1) 区間 | $g > 0$ なら $w_\mathrm{hi} \leftarrow w$、それ以外は $w_\mathrm{lo} \leftarrow w$ |
| 更新 (2) 候補 | $w_\mathrm{N} = w - g / (\mathrm{d}g/\mathrm{d}w)$ |
| 更新 (3) 採否 | $w_\mathrm{lo} < w_\mathrm{N} < w_\mathrm{hi}$ なら $w \leftarrow w_\mathrm{N}$、それ以外は $w \leftarrow (w_\mathrm{lo} + w_\mathrm{hi})/2$ |
| $q \le \sigma_{y0}$ の点 | $g = -1$ として扱う ($g < 0$ の側)。Newton 法の候補は作らず、次の点は区間の中点である。$q$ が NaN の点も同じ扱いである |
| $\mathrm{d}g/\mathrm{d}w = 0$ の点 | Newton 法の候補は作らず、次の点は区間の中点である |
| 収束の判定 | $\lvert g\rvert < 10^{-12}$。初期値でも判定する |
| 反復の上限 | 更新 100 回。超えたステップは非収束である |
| 収束後の確認 | $\Delta\lambda < 0$、または $\Delta\lambda$・$\boldsymbol\sigma_{n+1}$・$\boldsymbol\alpha_{n+1}$ のいずれかが有限でないステップは非収束である |
| 非収束時の出力 | 応力・`alpha`・`ep` はステップ開始時の値、$\Delta\lambda = 0$、接線は弾性剛性 |
| 報告する残差 | 初期値と各更新の後の $\lvert g\rvert$ ($q \le \sigma_{y0}$ の点は 1) |

1 ステップの更新の回数の実測値は次のとおりである (塑性ステップ。標準の載荷履歴と、「性質」の P1、P3、P7 の載荷)。

| パラメータ | 更新の回数 |
|---|---|
| 目安の値 ($\gamma = 100$) | 3〜6 回 |
| $\gamma = 1$ | 2〜3 回 |
| $\gamma = 0$ | 0 回 |

これらの載荷では、中点への切り替えと $q \le \sigma_{y0}$ の点は生じない。

### 応力と接線

収束した $w$ での $\boldsymbol{n}$、$\Delta\lambda$、$\boldsymbol{b}$ から、次のとおり更新する。

$$
\boldsymbol\sigma_{n+1} = \boldsymbol\sigma^\mathrm{tr} - 2\mu\,\Delta\lambda\,\boldsymbol{n},
\qquad
\boldsymbol\alpha_{n+1} = w\,\boldsymbol{b},
\qquad
p_{n+1} = p_n + \Delta\lambda
$$

接線は次のとおりである。

$$
\mathbf{D} = \mathbf{C}
- \frac{3\mu\,\Delta\lambda}{q}\,\mathbf{I}_\mathrm{dev}\mathbf{C}
- 4\mu^2 \left( \frac{1}{M} - \frac{\Delta\lambda}{q} \right) \boldsymbol{n}\otimes\boldsymbol{n}
- \frac{2\mu}{\mathrm{d}g/\mathrm{d}w}\,\boldsymbol{h}\otimes\boldsymbol{G}
$$

$$
\boldsymbol{h} = \left[ \frac{2\mu\,(c + C_k\,\Delta\lambda)}{M} - \frac{2\mu\,\Delta\lambda\,c}{q} \right] \boldsymbol{n} + \frac{3\mu\,\Delta\lambda}{q}\,\boldsymbol\alpha_n
$$

$$
\boldsymbol{G} = \gamma\,w^2 \left[ \frac{A + k_1}{M} - \frac{k_2\,c}{q} \right] \boldsymbol{n} + \frac{3}{2}\,\gamma\,w^2\,\frac{k_2}{q}\,\boldsymbol\alpha_n
$$

| 項 | 内容 |
|---|---|
| 第 1〜3 項 | $w$ を固定したときの接線 |
| 第 4 項 | $w$ の変化の寄与。$\boldsymbol{G} = \partial g/\partial\boldsymbol{s}^\mathrm{tr}$ ($w$ を固定)、$\mathrm{d}w = -(\boldsymbol{G}:\mathrm{d}\boldsymbol{s}^\mathrm{tr})/(\mathrm{d}g/\mathrm{d}w)$、$\boldsymbol{h} = \partial\boldsymbol\sigma_{n+1}/\partial w$ ($\boldsymbol{s}^\mathrm{tr}$ を固定) |
| 成分 | $\boldsymbol{n}$、$\boldsymbol\alpha_n$、$\boldsymbol{h}$、$\boldsymbol{G}$ は応力系の成分 (せん断はテンソル成分) で、$\otimes$ は成分どうしの積 $a_i b_j$ である。行は応力の成分、列はひずみの成分 (工学せん断ひずみ) に対応する。偏差テンソル $\boldsymbol{t}$ に対して $\boldsymbol{t}:\mathrm{d}\boldsymbol{s}^\mathrm{tr} = 2\mu \sum_i t_i\,\mathrm{d}\varepsilon_i$ が成り立つので、せん断成分に係数は付かない |
| 評価点 | すべての量は収束した $w$ での値である。UMAT は反復の $w$ を使う。Python の接線は、$w = 1/(1 + \gamma\,\Delta\lambda\,\bar\alpha_{n+1})$ として収束した状態から求める |
| 対称性 | $\gamma > 0$ では、$\boldsymbol{h}\otimes\boldsymbol{G}$ のため一般に非対称である |
| $\gamma = 0$ | $\boldsymbol{G} = \boldsymbol{0}$ で、第 4 項は消える |
| 非収束のステップ | 弾性剛性 $\mathbf{C}$ |

## 近似

無い。手書きの解と接線は厳密である。検査の許容値と載荷履歴は、既定のものである。

## UMAT

`OWKinematic3D` の手書きの解法の移植である。

| 項目 | 内容 |
|---|---|
| ソース | `fortran/ow_kinematic_3d.f90` |
| f2py のモジュール名 | `ow_kinematic_3d` |
| 中核のサブルーチン | `ow_kinematic_3d`。1 ステップの更新を行う |
| 中核の入力 | パラメータ (表の順)、応力、`alpha`、`ep`、ひずみ増分 |
| 中核の出力 | 応力、`alpha`、`ep`、接線、収束フラグ (整数、1 = 収束、0 = 非収束) |
| ABAQUS の入口 | `umat`。`PROPS` と `STATEV` を展開して中核を呼ぶ |
| `PROPS` | 1 = `E`、2 = `nu`、3 = `sigma_y0`、4 = `C_k`、5 = `gamma` |
| `STATEV` | 1〜6 = `alpha` (成分の順序は上記)、7 = `ep` |
| 反復の上限と収束の判定 | 手書きの解法と同じ |
| 非収束の条件 | 反復の上限を超えた場合、$\Delta\lambda < 0$ の場合、解が有限でない場合 |
| 非収束時の中核の出力 | 応力・`alpha`・`ep` は入力の値、接線は弾性剛性、収束フラグは 0 |
| 非収束時の入口 | `PNEWDT = 0.5` とし、`STRESS` と `STATEV` を変更しない |
| 収束時の入口 | `PNEWDT` を変更しない |
| 精度 | 倍精度。実数の定数は倍精度で書く |

## 検証

### 対象

| クラス | パラメータ |
|---|---|
| `OWKinematic3D` | 目安の値。目安の値の `gamma` を 1 に替えた値 |
| `OWKinematicPS` | 同上 |
| `OWKinematic1D` | 同上 |

目安の値 (`gamma = 100`) では、背応力は 2〜3 ステップで飽和する。`gamma = 1` では、発展則の過渡が載荷履歴の全体に広がる。

### 報告

共通の検査の報告である。UMAT をビルドした環境 (コンテナ) で、`returnmap check --family returnmap.models.ow_kinematic --require-built --range` が表示したものである。

```none
OWKinematic3D: every check that applies passes
  ok   converges                 59 steps over 2 histories, 47 plastic
  ok   tangent                   9.6e-10 (tol 1e-05) at 8 plastic steps
  —    hand blocks               no @hand_block on the model
  ok   hand-written solution     7.4e-14 (tol 1e-08) over 4 rows at 47 plastic steps
  ok   hand-written vs autodiff  stress 1.7e-11, state 2.9e-09 (tol 1e-06), tangent 2.5e-11 (tol 1e-05) over 2 histories
  ok   hand-written tangent      4.4e-10 (tol 1e-05) at 8 plastic steps
  ok   UMAT                      stress 3.0e-16, state 1.2e-15 (tol 1e-10), tangent 2.4e-15 (tol 1e-09) over 2 histories (ow_kinematic_3d)
  ok   UMAT entry                0.0e+00 (tol 1e-12) against ow_kinematic_3d over 2 histories (ow_kinematic_3d.umat)
  histories: uniaxial (24 plastic steps), multiaxial (23 plastic steps)
OWKinematic3D over its range: 1 case(s) failed — E=1000
  parameters — check_model at each end of every documented range
    FAIL E=1000              converges: uniaxial: stopped at step 5; multiaxial: stopped at step 0; tangent: no plastic step to check a tangent at; hand-written vs autodiff: stress 9.4e-02, state 7.3e+02 (tol 1e-06), tangent 3.8e-01 (tol 1e-05); not converged at uniaxial step 5, multiaxial step 0 over 2 histories
    ok   E=1e+06             7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   nu=0                7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   nu=0.49             7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   sigma_y0=1          7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   sigma_y0=5000       7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   C_k=0               7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   C_k=1e+06           7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   gamma=0             7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   gamma=10000         7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
  loading — random strain histories
    ok   random              36 of 36 increments solved, 35 plastic; solution 2.3e-14, hand-written vs autodiff 2.8e-09, tangent 2.4e-10, UMAT 2.6e-15, UMAT entry 0.0e+00
    ok   random large        24 of 24 increments solved, 24 plastic; solution 1.1e-14, hand-written vs autodiff 8.6e-10, tangent 1.3e-10, UMAT 2.8e-15, UMAT entry 0.0e+00
    ok   random single-step  6 of 6 increments solved, 6 plastic; solution 8.4e-15, hand-written vs autodiff 2.3e-11, tangent 7.5e-10, UMAT 4.2e-16, UMAT entry 0.0e+00

OWKinematic3D (gamma=1):
OWKinematic3D: every check that applies passes
  ok   converges                 59 steps over 2 histories, 47 plastic
  ok   tangent                   4.3e-10 (tol 1e-05) at 8 plastic steps
  —    hand blocks               no @hand_block on the model
  ok   hand-written solution     6.6e-13 (tol 1e-08) over 4 rows at 47 plastic steps
  ok   hand-written vs autodiff  stress 3.9e-12, state 3.2e-10 (tol 1e-06), tangent 9.4e-12 (tol 1e-05) over 2 histories
  ok   hand-written tangent      4.3e-10 (tol 1e-05) at 8 plastic steps
  ok   UMAT                      stress 2.7e-16, state 7.8e-16 (tol 1e-10), tangent 2.4e-14 (tol 1e-09) over 2 histories (ow_kinematic_3d)
  ok   UMAT entry                0.0e+00 (tol 1e-12) against ow_kinematic_3d over 2 histories (ow_kinematic_3d.umat)
  histories: uniaxial (24 plastic steps), multiaxial (23 plastic steps)
OWKinematic3D over its range: every case passes
  parameters — check_model at each end of every documented range
    ok   E=1000              7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   E=1e+06             7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   nu=0                7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   nu=0.49             7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   sigma_y0=1          7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   sigma_y0=5000       7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   C_k=0               7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   C_k=1e+06           7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   gamma=0             7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
    ok   gamma=10000         7 checks pass (converges, tangent, hand-written solution, hand-written vs autodiff, hand-written tangent, UMAT, UMAT entry)
  loading — random strain histories
    ok   random              36 of 36 increments solved, 35 plastic; solution 2.6e-13, hand-written vs autodiff 3.7e-11, tangent 2.3e-10, UMAT 2.5e-14, UMAT entry 0.0e+00
    ok   random large        24 of 24 increments solved, 24 plastic; solution 7.1e-15, hand-written vs autodiff 2.5e-09, tangent 8.6e-11, UMAT 5.7e-15, UMAT entry 0.0e+00
    ok   random single-step  6 of 6 increments solved, 6 plastic; solution 1.3e-14, hand-written vs autodiff 8.5e-13, tangent 5.9e-10, UMAT 1.3e-15, UMAT entry 0.0e+00

OWKinematicPS: every check that applies passes
  ok   converges                 59 steps over 2 histories, 46 plastic
  ok   tangent                   2.1e-09 (tol 1e-05) at 8 plastic steps
  —    hand blocks               no @hand_block on the model
  —    hand-written solution     no user_defined_return_mapping on the model
  —    hand-written vs autodiff  no user_defined_return_mapping on the model
  —    hand-written tangent      no user_defined_return_mapping on the model
  —    UMAT                      no UMAT declared for OWKinematicPS in returnmap.models.ow_kinematic.verification
  —    UMAT entry                no UMAT declared for OWKinematicPS
  histories: uniaxial (24 plastic steps), multiaxial (22 plastic steps)
OWKinematicPS over its range: 2 case(s) failed — E=1000, C_k=1e+06
  parameters — check_model at each end of every documented range
    FAIL E=1000              converges: uniaxial: stopped at step 5; multiaxial: stopped at step 1; tangent: no plastic step to check a tangent at
    ok   E=1e+06             2 checks pass (converges, tangent)
    ok   nu=0                2 checks pass (converges, tangent)
    ok   nu=0.49             2 checks pass (converges, tangent)
    ok   sigma_y0=1          2 checks pass (converges, tangent)
    ok   sigma_y0=5000       2 checks pass (converges, tangent)
    ok   C_k=0               2 checks pass (converges, tangent)
    FAIL C_k=1e+06           converges: multiaxial: stopped at step 1
    ok   gamma=0             2 checks pass (converges, tangent)
    ok   gamma=10000         2 checks pass (converges, tangent)
  loading — random strain histories
    ok   random              36 of 36 increments solved, 33 plastic; tangent 9.1e-10
    ok   random large        24 of 24 increments solved, 24 plastic; tangent 4.0e-10
    ok   random single-step  6 of 6 increments solved, 6 plastic; tangent 1.3e-09

OWKinematicPS (gamma=1):
OWKinematicPS: every check that applies passes
  ok   converges                 59 steps over 2 histories, 46 plastic
  ok   tangent                   1.5e-09 (tol 1e-05) at 8 plastic steps
  —    hand blocks               no @hand_block on the model
  —    hand-written solution     no user_defined_return_mapping on the model
  —    hand-written vs autodiff  no user_defined_return_mapping on the model
  —    hand-written tangent      no user_defined_return_mapping on the model
  —    UMAT                      no UMAT declared for OWKinematicPS in returnmap.models.ow_kinematic.verification
  —    UMAT entry                no UMAT declared for OWKinematicPS
  histories: uniaxial (24 plastic steps), multiaxial (22 plastic steps)
OWKinematicPS over its range: 1 case(s) failed — E=1000
  parameters — check_model at each end of every documented range
    FAIL E=1000              converges: multiaxial: stopped at step 1
    ok   E=1e+06             2 checks pass (converges, tangent)
    ok   nu=0                2 checks pass (converges, tangent)
    ok   nu=0.49             2 checks pass (converges, tangent)
    ok   sigma_y0=1          2 checks pass (converges, tangent)
    ok   sigma_y0=5000       2 checks pass (converges, tangent)
    ok   C_k=0               2 checks pass (converges, tangent)
    ok   C_k=1e+06           2 checks pass (converges, tangent)
    ok   gamma=0             2 checks pass (converges, tangent)
    ok   gamma=10000         2 checks pass (converges, tangent)
  loading — random strain histories
    ok   random              36 of 36 increments solved, 33 plastic; tangent 5.2e-10
    ok   random large        24 of 24 increments solved, 24 plastic; tangent 9.2e-11
    ok   random single-step  6 of 6 increments solved, 6 plastic; tangent 3.1e-08

OWKinematic1D: every check that applies passes
  ok   converges                 35 steps over 1 histories, 24 plastic
  ok   tangent                   1.1e-06 (tol 1e-05) at 4 plastic steps
  —    hand blocks               no @hand_block on the model
  —    hand-written solution     no user_defined_return_mapping on the model
  —    hand-written vs autodiff  no user_defined_return_mapping on the model
  —    hand-written tangent      no user_defined_return_mapping on the model
  —    UMAT                      no UMAT declared for OWKinematic1D in returnmap.models.ow_kinematic.verification
  —    UMAT entry                no UMAT declared for OWKinematic1D
  histories: uniaxial (24 plastic steps)
OWKinematic1D over its range: 2 case(s) failed — E=1000, gamma=10000
  parameters — check_model at each end of every documented range
    FAIL E=1000              converges: uniaxial: stopped at step 5; tangent: no plastic step to check a tangent at
    ok   E=1e+06             2 checks pass (converges, tangent)
    ok   nu=0                2 checks pass (converges, tangent)
    ok   nu=0.49             2 checks pass (converges, tangent)
    ok   sigma_y0=1          2 checks pass (converges, tangent)
    ok   sigma_y0=5000       2 checks pass (converges, tangent)
    ok   C_k=0               2 checks pass (converges, tangent)
    ok   C_k=1e+06           2 checks pass (converges, tangent)
    ok   gamma=0             2 checks pass (converges, tangent)
    FAIL gamma=10000         tangent: 1.9e-05 (tol 1e-05) at 4 plastic steps
  loading — random strain histories
    ok   random              36 of 36 increments solved, 16 plastic; tangent 9.6e-09
    ok   random large        24 of 24 increments solved, 23 plastic; tangent 7.7e-07
    ok   random single-step  6 of 6 increments solved, 6 plastic; tangent 1.7e-06

OWKinematic1D (gamma=1):
OWKinematic1D: every check that applies passes
  ok   converges                 35 steps over 1 histories, 24 plastic
  ok   tangent                   9.0e-10 (tol 1e-05) at 4 plastic steps
  —    hand blocks               no @hand_block on the model
  —    hand-written solution     no user_defined_return_mapping on the model
  —    hand-written vs autodiff  no user_defined_return_mapping on the model
  —    hand-written tangent      no user_defined_return_mapping on the model
  —    UMAT                      no UMAT declared for OWKinematic1D in returnmap.models.ow_kinematic.verification
  —    UMAT entry                no UMAT declared for OWKinematic1D
  histories: uniaxial (24 plastic steps)
OWKinematic1D over its range: 1 case(s) failed — gamma=10000
  parameters — check_model at each end of every documented range
    ok   E=1000              2 checks pass (converges, tangent)
    ok   E=1e+06             2 checks pass (converges, tangent)
    ok   nu=0                2 checks pass (converges, tangent)
    ok   nu=0.49             2 checks pass (converges, tangent)
    ok   sigma_y0=1          2 checks pass (converges, tangent)
    ok   sigma_y0=5000       2 checks pass (converges, tangent)
    ok   C_k=0               2 checks pass (converges, tangent)
    ok   C_k=1e+06           2 checks pass (converges, tangent)
    ok   gamma=0             2 checks pass (converges, tangent)
    FAIL gamma=10000         tangent: 1.9e-05 (tol 1e-05) at 4 plastic steps
  loading — random strain histories
    ok   random              36 of 36 increments solved, 16 plastic; tangent 2.2e-08
    ok   random large        24 of 24 increments solved, 23 plastic; tangent 2.9e-10
    ok   random single-step  6 of 6 increments solved, 6 plastic; tangent 2.0e-07

6 model(s) checked: 5 failed — OWKinematic3D, OWKinematicPS, OWKinematicPS (gamma=1), OWKinematic1D, OWKinematic1D (gamma=1)
```

「性質」の P1〜P12 のテスト (`tests/test_properties.py`) は、全部合格である。

### 範囲

| 項目 | 状態 |
|---|---|
| 自動微分の経路の収束と接線、手書きの解と接線、UMAT の中核、ABAQUS の入口 | `tests/test_verification.py` (標準の載荷履歴。「対象」の 6 モデル)。UMAT の 2 行は、ビルドした環境でだけ実行される |
| 性質 P1〜P12 | `tests/test_properties.py`。期待値は、テストの中で式から計算する |
| パラメータの範囲の端、乱数の載荷履歴 | テストには無い。上の報告のとおり、`returnmap check --range` で確かめた |
| 移動硬化係数がヤング率の数倍を超える組 | 自動微分の経路が収束せず、報告の `E=1000` と、平面応力の `C_k=1e+06` が不合格になる。`C_k` が `E` の 2 倍までは収束し、4 倍では収束しない。パラメータの範囲は 1 つずつ与えてあり、この組み合わせを除いていない。3D の手書きの経路は、この組でも解ける (不合格になるのは、自動微分の経路の収束の行と、それと比べる行である)。[#160](https://github.com/toruinaba/returnmap/issues/160) |
| 単軸の `gamma=10000` | 接線が中心差分と 1.9e-5 の差になり、不合格になる (許容値は 1e-5)。接線は正しく (0.095 MPa)、差分を取る応力が NR の許容値までしか解かれていないためである (許容値を 1e-13 にすると、差は 9.4e-7 になる)。[#157](https://github.com/toruinaba/returnmap/issues/157) |
| 手書きの解法の、実行されない分岐 ($q \le \sigma_{y0}$ の点、区間の中点への切り替え) | 標準の載荷履歴と「性質」の載荷では生じない。乱数の載荷履歴では、中点への切り替えが生じる |
| 反復の上限に達するステップ | 確かめていない |
| ABAQUS での実行 | 確かめていない |
