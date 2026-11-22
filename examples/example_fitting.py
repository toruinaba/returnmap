"""フィッティングの設定 — fit_params の探索空間、ターゲット、比較の方法、driver を指定する。

文書の「フィッティングの設定」(docs/guide/fitting.md) のコード。文書のページは、
このファイルの各節 (``# %%`` の区切り) をそのまま載せている。

    uv run python examples/example_fitting.py

実行時間は約 5 秒。記録は、既知のパラメータの J2 モデル (sigma_y0 = 250、H = 1000) の応答に
ノイズを加えて作る。フィッティングの結果は、この値に戻る。
"""

# %% 1. 基本の呼び出し
import numpy as np

from returnmap.fitting import FitTarget, ParameterSpace, fit_params
from returnmap.models import J2Isotropic1D
from returnmap.simulation import FieldHistory, PythonIntegrator, StrainDriver

true_model = J2Isotropic1D(E=210_000.0, nu=0.3, sigma_y0=250.0, H=1_000.0)
rng = np.random.default_rng(42)


def record(load):
    """true_model の応答に、標準偏差 0.5 MPa のノイズを加えた応力の記録。"""
    stress = StrainDriver(PythonIntegrator(true_model)).run(load).stress[:, 0]
    return stress + rng.normal(0.0, 0.5, len(stress))


load = FieldHistory.cyclic_strain([5e-3, -5e-3, 5e-3, -5e-3], n_per_segment=12)
strain_exp, stress_exp = load.data[:, 0], record(load)

space = ParameterSpace(
    free={                                          # 名前: (初期値, (下限, 上限))
        "sigma_y0": (180.0, (50.0, 600.0)),
        "H": (500.0, (0.0, 10_000.0)),
    },
    fixed={"E": 210_000.0, "nu": 0.3},
)
result = fit_params(
    true_model,                                     # テンプレート (パラメータの値は用いない)
    [FitTarget.stress_strain(strain_exp, stress_exp)],
    space,
    method="L-BFGS-B",
    max_iter=200,
)
print({name: round(value, 1) for name, value in result.model_params.items()})
print(f"残差 {result.residual:.3f}、収束 {result.success}、評価 {result.n_evaluations} 回、失敗 {result.n_failed} 回")
print(type(result.model).__name__)                  # 同定したパラメータを入れたモデル

# %% 2. パラメータ間の制約
ratio_space = ParameterSpace(
    free={
        "sigma_y0": (180.0, (50.0, 600.0)),
        "H_ratio": (2.0, (1.0, 20.0)),              # H / sigma_y0 >= 1
    },
    fixed={"E": 210_000.0, "nu": 0.3},
    derived={"H": lambda p: p["sigma_y0"] * p["H_ratio"]},
)
by_ratio = fit_params(true_model, [FitTarget.stress_strain(strain_exp, stress_exp)], ratio_space)
print(sorted(by_ratio.params))                      # 探索空間のすべての名前
print(sorted(by_ratio.model_params))                # モデルのパラメータ
print(f"H_ratio {by_ratio.params['H_ratio']:.2f}、H {by_ratio.model_params['H']:.1f}")

# %% 3. 複数の試験
from returnmap.fitting import Observation
from returnmap.simulation import FieldType, StressDriver

# 試験 1: ひずみ制御。同じ解析から、応力と累積塑性ひずみの 2 つの記録を比較する
load_1 = FieldHistory.cyclic_strain([3e-3, -3e-3], n_per_segment=10)
run_1 = StrainDriver(PythonIntegrator(true_model)).run(load_1)

# 試験 2: 応力制御。ひずみの記録を比較する
load_2 = FieldHistory.cyclic_stress([260.0, -260.0], n_per_segment=10)
run_2 = StressDriver(PythonIntegrator(true_model)).run(load_2)

targets = [
    FitTarget(
        load=load_1,
        observations=[
            Observation(values=run_1.stress[:, 0], component=0),
            Observation(values=run_1.history("ep"), field="ep", weight=1e8),
        ],
        name="strain-controlled",
    ),
    FitTarget(
        load=load_2,
        observations=[Observation(values=run_2.strain[:, 0], field="strain", component=0, weight=1e8)],
        name="stress-controlled",
    ),
]
multi = fit_params(true_model, targets, space)
print({name: round(value, 1) for name, value in multi.model_params.items()})

# %% 4. サンプリングが一致しない記録
dense = FieldHistory.cyclic_strain([5e-3], n_per_segment=200)       # 記録: 200 点
coarse = FieldHistory.cyclic_strain([5e-3], n_per_segment=25)       # 計算: 25 ステップ
stress_dense = StrainDriver(PythonIntegrator(true_model)).run(dense).stress[:, 0]

pcm_target = FitTarget(
    load=coarse,
    observations=[Observation(
        values=stress_dense,
        x_values=dense.data[:, 0],                  # 記録の横軸 (pcm では必須)
        metric="pcm",
        component=0,
    )],
)
by_pcm = fit_params(
    true_model, [pcm_target],
    ParameterSpace(free={"H": (500.0, (0.0, 10_000.0))},
                   fixed={"E": 210_000.0, "nu": 0.3, "sigma_y0": 250.0}),
)
print(f"H {by_pcm.model_params['H']:.0f}")

# %% 5. driver の指定
from returnmap.models import J2Isotropic3D
from returnmap.simulation import MixedDriver

# 3D のモデルを、単軸応力の記録に合わせる
by_3d = fit_params(
    J2Isotropic3D(**true_model.params),
    [FitTarget.stress_strain(strain_exp, stress_exp)],
    space,
    driver_factory=lambda integrator: MixedDriver(integrator, prescribed_strain_idx=[0]),
)
print({name: round(value, 1) for name, value in by_3d.model_params.items()})
