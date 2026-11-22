"""サービスに組み込む — レジストリにモデルを登録し、plain な入出力で計算する。

文書の「サービスに組み込む」(docs/guide/service.md) のコード。文書のページは、
このファイルの各節 (``# %%`` の区切り) をそのまま載せている。

    uv run python examples/example_service.py
"""

# %% 1. レジストリ
from returnmap.service import default_registry

registry = default_registry()                       # 組み込みの 1D モデルを登録したレジストリ
print(registry.names())

description = registry.describe("j2_isotropic")     # GET /models に相当する情報
print(description.class_name, description.ntens, description.param_names)
for name in description.param_names:
    meta = description.param_meta[name]
    print(f"{name:9s} {meta.unit:4s} {meta.default:>9}  {meta.bounds}  {meta.doc}")

# %% 2. 計算の実行
import numpy as np

from returnmap.service import run_simulation

params = {name: meta.default for name, meta in description.param_meta.items()}   # 入力フォームの初期値

progress = []
out = run_simulation(
    registry.get("j2_isotropic"),
    params,
    load=np.linspace(0.0, 5e-3, 200),               # ひずみ履歴 (control="stress" で応力制御)
    include_fields=["ep"],                          # 状態量は指定したものだけを返す
    progress=lambda done, total: progress.append((done, total)),
)
print(out.converged, out.n_steps, out.stress.shape, out.fields["ep"].shape)
print(f"最終の応力 {out.stress[-1, 0]:.1f} MPa、進捗の通知 {len(progress)} 回 (最後は {progress[-1]})")

# %% 3. 非収束と入力の誤り
unsolvable = dict(E=206e3, nu=0.3, sigma_y0=250.0, C_k=1e9, gamma=1e12)
partial = run_simulation(registry.get("ow_kinematic"), unsolvable, np.linspace(0.0, 0.05, 10))
print(partial.converged, partial.failed_step, partial.n_steps, partial.stress.shape)

try:
    run_simulation(registry.get("j2_isotropic"), {"E": 206e3}, np.linspace(0.0, 1e-3, 5))
except ValueError as error:
    print(type(error).__name__, "—", error)

# %% 4. テンプレートの登録
from returnmap.core import PLANE_STRESS
from returnmap.models import J2Isotropic
from returnmap.service import ModelRegistry

values = {"E": 206e3, "nu": 0.3, "sigma_y0": 250.0, "H": 1000.0}
own = ModelRegistry()
own.register("j2_3d", J2Isotropic(**values))                             # 次元の既定は 3D
own.register("j2_ps", J2Isotropic(dimension=PLANE_STRESS, **values))     # 同じクラスを別の次元で
print(own.names(), own.describe("j2_3d").ntens, own.describe("j2_ps").ntens)
