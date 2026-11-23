"""The properties of the formulation — ``README.md``, section "性質", P1 to P12.

One test per property.  Each follows from the equations and holds whatever
the implementation, so the expected values are written here from the
equations, not taken from the model's own methods.  The error is
``max|computed − expected| / max|expected|`` over the whole array of the steps
the property names, and the tolerances are the README's.  The parameters are
the documented typical ones and the typical ones with ``gamma = 1`` (P9 and
P10 have their own: ``gamma = 0``).
"""

import numpy as np
import pytest

from returnmap.models import (
    AFKinematic1D, AFKinematic3D, OWKinematic1D, OWKinematic3D, OWKinematicPS,
)
from returnmap.simulation import (
    FieldHistory,
    MixedDriver,
    PythonAnalyticalIntegrator,
    PythonIntegrator,
    PythonNumericalIntegrator,
    StrainDriver,
)
from returnmap.verification import standard_histories


def rel_err(computed, expected) -> float:
    """``max|computed − expected| / max|expected|`` over the whole array."""
    computed = np.asarray(computed, dtype=float)
    expected = np.asarray(expected, dtype=float)
    return float(np.max(np.abs(computed - expected)) / np.max(np.abs(expected)))


def plastic_steps(result) -> np.ndarray:
    """Mask of the plastic steps of a run."""
    return np.array([r.is_plastic is True for r in result.step_results])


def tangents(result) -> np.ndarray:
    """Consistent tangent of every step, shape ``(N, ntens, ntens)``."""
    return np.array([np.asarray(r.ddsdde, dtype=float) for r in result.step_results])


def uniaxial_backstress(result) -> np.ndarray:
    """a = (3/2) α11: the stress at the centre of the yield surface, per step."""
    return 1.5 * result.history("alpha")[:, 0]


# The two parameter sets of the README: the typical values, and the typical
# values with gamma = 1.
PARAMS = {"typical": {}, "gamma=1": {"gamma": 1.0}}

# A property of OWKinematic3D holds on both of its paths: the hand-written
# return mapping and tangent, and the autodiff ones.
PATHS_3D = {
    "hand-written": PythonAnalyticalIntegrator,
    "autodiff": PythonNumericalIntegrator,
}

SOLID_CASES = [
    pytest.param(integrator_cls, params, id=f"OWKinematic3D {path}, {name}")
    for path, integrator_cls in PATHS_3D.items()
    for name, params in PARAMS.items()
]


# ---------------------------------------------------------------------------
# Uniaxial stress: OWKinematic1D, and OWKinematic3D with only ε11 prescribed
# ---------------------------------------------------------------------------

def uniaxial_1d(params):
    """``OWKinematic1D`` driven by its one strain component."""
    def run(axial: FieldHistory):
        model = OWKinematic1D.at_defaults(**params)
        return model, StrainDriver(PythonNumericalIntegrator(model)).run(axial)
    return run


def uniaxial_3d(integrator_cls, params):
    """``OWKinematic3D`` under uniaxial stress: ε11 prescribed, the other stresses zero."""
    def run(axial: FieldHistory):
        model = OWKinematic3D.at_defaults(**params)
        driver = MixedDriver(integrator_cls(model), prescribed_strain_idx=[0])
        return model, driver.run(axial)
    return run


BAR_CASES = [
    pytest.param(uniaxial_1d(params), id=f"OWKinematic1D, {name}")
    for name, params in PARAMS.items()
]
# The autodiff path under mixed control takes a few seconds a case (an outer
# iteration on the lateral strains around every stress update): opt-in.
UNIAXIAL_STRESS_3D = [
    pytest.param(
        uniaxial_3d(integrator_cls, params), id=f"OWKinematic3D {path}, {name}",
        marks=[pytest.mark.slow] if path == "autodiff" else [],
    )
    for path, integrator_cls in PATHS_3D.items()
    for name, params in PARAMS.items()
]
UNIAXIAL_STRESS = BAR_CASES + UNIAXIAL_STRESS_3D

# Axial strain 0 → 0.01 → −0.01 → 0.01, 20 steps per leg (P1, P11).
REVERSED = FieldHistory.cyclic_strain([0.01, -0.01, 0.01], n_per_segment=20)
# Axial strain 0 → 0.02 in 40 steps (P3, P4, P10).
MONOTONIC = FieldHistory.cyclic_strain([0.02], n_per_segment=40)


# ---------------------------------------------------------------------------
# P1 — on the yield surface at every plastic step
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("run", UNIAXIAL_STRESS)
def test_p1_axial_stress_is_one_yield_stress_from_the_backstress(run):
    model, result = run(REVERSED)
    plastic = plastic_steps(result)
    assert plastic.any()
    stress_11 = result.stress[:, 0]
    a = uniaxial_backstress(result)
    assert rel_err(
        np.abs(stress_11 - a)[plastic], np.full(plastic.sum(), model.sigma_y0)
    ) <= 1e-8


# ---------------------------------------------------------------------------
# P2 — saturation of the backstress
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("run", BAR_CASES)
def test_p2_axial_stress_saturates_at_sigma_y0_plus_sqrt_ck_over_gamma(run):
    model, result = run(FieldHistory.cyclic_strain([0.2], n_per_segment=200))
    saturated = model.sigma_y0 + np.sqrt(model.C_k / model.gamma)
    assert rel_err(result.stress[-1, 0], saturated) <= 1e-8


# ---------------------------------------------------------------------------
# P3 — the plastic part of the axial strain is p
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("run", UNIAXIAL_STRESS)
def test_p3_plastic_axial_strain_is_p(run):
    model, result = run(MONOTONIC)
    plastic_strain_11 = result.strain[:, 0] - result.stress[:, 0] / model.E
    assert rel_err(plastic_strain_11, result.history("ep")) <= 1e-10


# ---------------------------------------------------------------------------
# P4 — the discrete evolution law of the uniaxial backstress
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("run", UNIAXIAL_STRESS)
def test_p4_backstress_increment_follows_the_backward_euler_law(run):
    model, result = run(MONOTONIC)
    plastic = plastic_steps(result)
    assert plastic.any()
    a = uniaxial_backstress(result)
    p = result.history("ep")
    # Increments over each step; the run starts from a = 0, p = 0.
    da = np.diff(a, prepend=0.0)
    dp = np.diff(p, prepend=0.0)
    expected = dp * (model.C_k - model.gamma * a ** 2)
    assert rel_err(da[plastic], expected[plastic]) <= 1e-6


# ---------------------------------------------------------------------------
# P5, P6 — the multiaxial strain history
# ---------------------------------------------------------------------------

def multiaxial(integrator_cls, params):
    model = OWKinematic3D.at_defaults(**params)
    load = standard_histories(model)["multiaxial"]
    return model, StrainDriver(integrator_cls(model)).run(load)


@pytest.mark.parametrize("integrator_cls, params", SOLID_CASES)
def test_p5_mean_stress_is_bulk_modulus_times_volumetric_strain(integrator_cls, params):
    model, result = multiaxial(integrator_cls, params)
    assert plastic_steps(result).any()
    K = model.E / (3.0 * (1.0 - 2.0 * model.nu))
    mean_stress = result.stress[:, :3].sum(axis=1) / 3.0
    volumetric_strain = result.strain[:, :3].sum(axis=1)
    assert rel_err(mean_stress, K * volumetric_strain) <= 1e-10


@pytest.mark.parametrize("integrator_cls, params", SOLID_CASES)
def test_p6_backstress_is_deviatoric(integrator_cls, params):
    _, result = multiaxial(integrator_cls, params)
    assert plastic_steps(result).any()
    alpha = result.history("alpha")
    trace = alpha[:, :3].sum(axis=1)
    assert np.max(np.abs(trace)) / np.max(np.abs(alpha)) <= 1e-10


# ---------------------------------------------------------------------------
# P7, P8 — simple shear: 2ε12 from 0 to 0.02 in 40 steps, nothing else
# ---------------------------------------------------------------------------

def simple_shear(integrator_cls, params):
    model = OWKinematic3D.at_defaults(**params)
    load = FieldHistory.cyclic_strain([0.02], n_per_segment=40, component=3, ntens=6)
    return model, StrainDriver(integrator_cls(model)).run(load)


@pytest.mark.parametrize("integrator_cls, params", SOLID_CASES)
def test_p7_shear_stress_is_yield_stress_over_sqrt3_from_the_backstress(integrator_cls, params):
    model, result = simple_shear(integrator_cls, params)
    plastic = plastic_steps(result)
    assert plastic.any()
    stress_12 = result.stress[:, 3]
    alpha_12 = result.history("alpha")[:, 3]
    assert rel_err(
        np.abs(stress_12 - alpha_12)[plastic],
        np.full(plastic.sum(), model.sigma_y0 / np.sqrt(3.0)),
    ) <= 1e-8


@pytest.mark.parametrize("integrator_cls, params", SOLID_CASES)
def test_p8_plastic_shear_strain_is_sqrt3_p(integrator_cls, params):
    model, result = simple_shear(integrator_cls, params)
    mu = model.E / (2.0 * (1.0 + model.nu))
    plastic_gamma_12 = result.strain[:, 3] - result.stress[:, 3] / mu
    assert rel_err(plastic_gamma_12, np.sqrt(3.0) * result.history("ep")) <= 1e-10


# ---------------------------------------------------------------------------
# P9 — without dynamic recovery, the model is Armstrong-Frederick with gamma = 0
# ---------------------------------------------------------------------------

P9_CASES = [
    pytest.param(OWKinematic3D, AFKinematic3D, PythonAnalyticalIntegrator,
                 id="OWKinematic3D hand-written"),
    pytest.param(OWKinematic3D, AFKinematic3D, PythonNumericalIntegrator,
                 id="OWKinematic3D autodiff"),
    pytest.param(OWKinematic1D, AFKinematic1D, PythonNumericalIntegrator,
                 id="OWKinematic1D"),
]


@pytest.mark.parametrize("ow_cls, reference_cls, integrator_cls", P9_CASES)
def test_p9_gamma_zero_is_af_kinematic_with_gamma_zero(ow_cls, reference_cls, integrator_cls):
    ow = ow_cls.at_defaults(gamma=0.0)
    reference = reference_cls(
        E=ow.E, nu=ow.nu, sigma_y0=ow.sigma_y0, C_k=ow.C_k, gamma=0.0)
    histories = standard_histories(ow)
    assert histories
    for name, load in histories.items():
        got = StrainDriver(integrator_cls(ow)).run(load)
        expected = StrainDriver(PythonIntegrator(reference)).run(load)
        assert plastic_steps(expected).any(), name
        assert rel_err(got.stress, expected.stress) <= 1e-8, name
        assert rel_err(got.history("alpha"), expected.history("alpha")) <= 1e-8, name
        assert rel_err(got.history("ep"), expected.history("ep")) <= 1e-8, name
        assert rel_err(tangents(got), tangents(expected)) <= 1e-8, name


# ---------------------------------------------------------------------------
# P10 — without dynamic recovery, the uniaxial tangent is that of linear hardening
# ---------------------------------------------------------------------------

def test_p10_gamma_zero_tangent_is_the_linear_hardening_modulus():
    model, result = uniaxial_1d({"gamma": 0.0})(MONOTONIC)
    plastic = plastic_steps(result)
    assert plastic.any()
    expected = model.E * model.C_k / (model.E + model.C_k)
    tangent_11 = tangents(result)[:, 0, 0]
    assert rel_err(tangent_11[plastic], np.full(plastic.sum(), expected)) <= 1e-8


# ---------------------------------------------------------------------------
# P11 — the 3D model under uniaxial stress is the 1D model
# ---------------------------------------------------------------------------

P11_CASES = [
    pytest.param(uniaxial_3d(integrator_cls, params), uniaxial_1d(params),
                 id=f"OWKinematic3D {path}, {name}")
    for path, integrator_cls in PATHS_3D.items()
    for name, params in PARAMS.items()
]


@pytest.mark.parametrize("run_solid, run_bar", P11_CASES)
def test_p11_uniaxial_stress_in_3d_is_the_1d_model(run_solid, run_bar):
    _, solid = run_solid(REVERSED)
    _, bar = run_bar(REVERSED)
    assert len(solid.step_results) == len(bar.step_results) == len(REVERSED.data)
    assert rel_err(solid.stress[:, 0], bar.stress[:, 0]) <= 1e-6
    assert rel_err(solid.history("ep"), bar.history("ep")) <= 1e-6
    assert rel_err(uniaxial_backstress(solid), uniaxial_backstress(bar)) <= 1e-6


# ---------------------------------------------------------------------------
# P12 — plane stress is the 3D model with the out-of-plane stresses held at zero
# ---------------------------------------------------------------------------

IN_PLANE = [0, 1, 3]            # 11, 22, 12 in the 3D Voigt order


@pytest.mark.parametrize("params", PARAMS.values(), ids=PARAMS.keys())
def test_p12_plane_stress_is_the_3d_model_with_zero_out_of_plane_stress(params):
    plane = OWKinematicPS.at_defaults(**params)
    solid = OWKinematic3D.at_defaults(**params)
    histories = standard_histories(plane)
    assert set(histories) == {"uniaxial", "multiaxial"}
    for name, load in histories.items():
        got = StrainDriver(PythonNumericalIntegrator(plane)).run(load)
        expected = MixedDriver(PythonIntegrator(solid), prescribed_strain_idx=IN_PLANE).run(load)
        assert len(got.step_results) == len(expected.step_results) == len(load.data), name
        assert plastic_steps(got).any(), name
        assert rel_err(got.stress, expected.stress[:, IN_PLANE]) <= 1e-6, name
        assert rel_err(got.history("alpha"), expected.history("alpha")[:, IN_PLANE]) <= 1e-6, name
        assert rel_err(got.history("ep"), expected.history("ep")) <= 1e-6, name
