"""検査を個別に実行する — check_model が行う検査を、1 つずつ呼び出す。

文書の「検査の一覧」(docs/guide/verification.md) のコード。文書のページは、
このファイルの各節 (``# %%`` の区切り) をそのまま載せている。

    uv run python examples/example_checks.py

実行時間は約 6 秒。対象は組み込みのモデルである: 手導出ブロックと手書きの解法を持つ
``YUKinematicPS``、閉形式の解を持つ ``J2Isotropic3D`` と ``J2Isotropic1D``、implicit の
状態を持つ ``OWKinematic3D``。
"""

# %% 1. まとめて検査する
from returnmap.models import YUKinematic1D
from returnmap.verification import check_model

report = check_model(YUKinematic1D.at_defaults())
assert report.passed, report.summary()
print(report.summary())

# %% 2. 残差のヤコビアン
import numpy as np

from returnmap.models import OWKinematic3D
from returnmap.simulation import PythonIntegrator

ow = OWKinematic3D.at_defaults()
integrator = PythonIntegrator(ow)
strain_inc = np.array([3e-3, 0.0, 0.0, 1e-3, 0.0, 0.0])          # 直ひずみとせん断
result = integrator.stress_update(strain_inc, np.zeros(6), ow.initial_state())

jac = integrator.jacobian(result)                               # 自動微分の残差ヤコビアン
print(list(jac.part))                                           # 行 (= 列) の名前
print(jac.part["alpha"]["stress"].shape)                        # ∂R_alpha/∂σ
print(np.round(jac.part["dlambda"]["stress"], 3))               # ∂f/∂σ

# %% 3. 手導出ブロックの照合
from returnmap.models import YUKinematicPS
from returnmap.simulation import StrainDriver
from returnmap.verification import check_hand_blocks, standard_histories

yu = YUKinematicPS.at_defaults()
load = standard_histories(yu)["multiaxial"]                     # 直ひずみとせん断を同時に動かす履歴

results = StrainDriver(PythonIntegrator(yu)).run(load).step_results
report = check_hand_blocks(yu, results)
assert report.passed, report.summary()
print(report.summary())

# %% 4. 接線の差分検査
from returnmap.models import J2Isotropic3D
from returnmap.verification import check_tangent

j2 = J2Isotropic3D.at_defaults()
report = check_tangent(PythonIntegrator(j2), strain_inc, np.zeros(6), j2.initial_state())
assert report.passed, report.summary()
print(report.summary())
print(report.rel_err_matrix.shape)                              # 成分ごとの誤差

# %% 5. 手書きの解の検査
from returnmap.simulation import PythonAnalyticalIntegrator
from returnmap.verification import check_solution

results = StrainDriver(PythonAnalyticalIntegrator(yu)).run(load).step_results   # 手書きの経路
report = check_solution(yu, results)
assert report.passed, report.summary()
print(report.summary())

# %% 6. 2 つの実装の比較
from returnmap.simulation import FieldHistory, FieldType, PythonNumericalIntegrator
from returnmap.verification import crosscheck, generate_strain_history

load = FieldHistory(FieldType.STRAIN, generate_strain_history(j2))
report = crosscheck(PythonNumericalIntegrator(j2), PythonAnalyticalIntegrator(j2), load)
assert report.passed, report.summary()
print(report.summary())
print(report.first_failure)                                     # 一致しなかった最初のステップ

# %% 8. パラメータの範囲と乱数の載荷
from returnmap.models import J2Isotropic1D
from returnmap.verification import check_range

report = check_range(J2Isotropic1D.at_defaults())
assert report.passed, report.summary()
print(report.summary())
print(report["H=0"].reports["model"]["hand-written tangent"].detail)   # 1 つの組の、1 つの行
