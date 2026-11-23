"""The properties of the formulation — ``README.md``, section "性質", P1 to P8.

One test per property.  Each follows from the equations and holds whatever
the implementation, so the expected values are written here from the
equations, not taken from the model's own methods.  The error is
``max|computed − expected| / max|expected|`` over the steps the property
names, at the documented typical parameters, and the tolerances are the
README's.
"""

import numpy as np
import pytest

from returnmap.models import J2Isotropic1D, J2Isotropic3D
from returnmap.simulation import (
    FieldHistory,
    MixedDriver,
    PythonAnalyticalIntegrator,
    PythonIntegrator,
    PythonNumericalIntegrator,
    StrainDriver,
)
from returnmap.verification import standard_histories
from mymodels.voce import J2Voce1D, J2Voce3D


def rel_err(computed, expected) -> float:
    """``max|computed − expected| / max|expected|`` over the whole array."""
    computed = np.asarray(computed, dtype=float)
    expected = np.asarray(expected, dtype=float)
    return float(np.max(np.abs(computed - expected)) / np.max(np.abs(expected)))


def sigma_y(model, p):
    """Voce yield stress σ_y(p) = σ_y0 + Q (1 − exp(−b p))."""
    return model.sigma_y0 + model.Q * (1.0 - np.exp(-model.b * np.asarray(p, dtype=float)))


def plastic_steps(result) -> np.ndarray:
    """Mask of the plastic steps of a run."""
    return np.array([r.is_plastic is True for r in result.step_results])


def tangents(result) -> np.ndarray:
    """Consistent tangent of every step, shape ``(N, ntens, ntens)``."""
    return np.array([np.asarray(r.ddsdde, dtype=float) for r in result.step_results])


# A property of J2Voce3D holds on both of its paths: the hand-written return
# mapping and tangent, and the autodiff ones.
PATHS_3D = [
    pytest.param(PythonAnalyticalIntegrator, id="J2Voce3D hand-written"),
    pytest.param(PythonNumericalIntegrator, id="J2Voce3D autodiff"),
]


# ---------------------------------------------------------------------------
# Uniaxial stress: J2Voce1D, and J2Voce3D with only the 11 strain prescribed
# ---------------------------------------------------------------------------

def uniaxial_1d(axial: FieldHistory):
    """``J2Voce1D`` driven by its one strain component."""
    model = J2Voce1D.at_defaults()
    return model, StrainDriver(PythonNumericalIntegrator(model)).run(axial)


def uniaxial_3d(integrator_cls):
    """``J2Voce3D`` under uniaxial stress: ε11 prescribed, the other stresses zero."""
    def run(axial: FieldHistory):
        model = J2Voce3D.at_defaults()
        driver = MixedDriver(integrator_cls(model), prescribed_strain_idx=[0])
        return model, driver.run(axial)
    return run


UNIAXIAL_STRESS_3D = [
    pytest.param(uniaxial_3d(PythonAnalyticalIntegrator), id="J2Voce3D hand-written"),
    pytest.param(uniaxial_3d(PythonNumericalIntegrator), id="J2Voce3D autodiff"),
]
UNIAXIAL_STRESS = [pytest.param(uniaxial_1d, id="J2Voce1D")] + UNIAXIAL_STRESS_3D

# Axial strain 0 → 0.01 → −0.01 → 0.01, 20 steps per leg (P1, P8).
REVERSED = FieldHistory.cyclic_strain([0.01, -0.01, 0.01], n_per_segment=20)


# ---------------------------------------------------------------------------
# P1 — on the yield surface at every plastic step
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("run", UNIAXIAL_STRESS)
def test_p1_axial_stress_is_the_yield_stress_on_plastic_steps(run):
    model, result = run(REVERSED)
    plastic = plastic_steps(result)
    assert plastic.any()
    stress_11 = result.stress[:, 0]
    p = result.history("ep")
    assert rel_err(np.abs(stress_11[plastic]), sigma_y(model, p[plastic])) <= 1e-8


# ---------------------------------------------------------------------------
# P2 — saturation
# ---------------------------------------------------------------------------

def test_p2_axial_stress_saturates_at_sigma_y0_plus_q():
    model, result = uniaxial_1d(FieldHistory.cyclic_strain([1.0], n_per_segment=200))
    assert rel_err(result.stress[-1, 0], model.sigma_y0 + model.Q) <= 1e-6


# ---------------------------------------------------------------------------
# P3 — the plastic part of the axial strain is p
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("run", UNIAXIAL_STRESS)
def test_p3_plastic_axial_strain_is_p(run):
    model, result = run(FieldHistory.cyclic_strain([0.02], n_per_segment=40))
    plastic_strain_11 = result.strain[:, 0] - result.stress[:, 0] / model.E
    assert rel_err(plastic_strain_11, result.history("ep")) <= 1e-10


# ---------------------------------------------------------------------------
# P4 — the volumetric response stays elastic
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("integrator_cls", PATHS_3D)
def test_p4_mean_stress_is_bulk_modulus_times_volumetric_strain(integrator_cls):
    model = J2Voce3D.at_defaults()
    result = StrainDriver(integrator_cls(model)).run(standard_histories(model)["multiaxial"])
    assert plastic_steps(result).any()
    K = model.E / (3.0 * (1.0 - 2.0 * model.nu))
    mean_stress = result.stress[:, :3].sum(axis=1) / 3.0
    volumetric_strain = result.strain[:, :3].sum(axis=1)
    assert rel_err(mean_stress, K * volumetric_strain) <= 1e-10


# ---------------------------------------------------------------------------
# P5, P6 — simple shear: γ12 from 0 to 0.02 in 40 steps, nothing else
# ---------------------------------------------------------------------------

def simple_shear(integrator_cls):
    model = J2Voce3D.at_defaults()
    load = FieldHistory.cyclic_strain([0.02], n_per_segment=40, component=3, ntens=6)
    return model, StrainDriver(integrator_cls(model)).run(load)


@pytest.mark.parametrize("integrator_cls", PATHS_3D)
def test_p5_shear_stress_is_the_yield_stress_over_sqrt3_on_plastic_steps(integrator_cls):
    model, result = simple_shear(integrator_cls)
    plastic = plastic_steps(result)
    assert plastic.any()
    stress_12 = result.stress[:, 3]
    p = result.history("ep")
    assert rel_err(
        np.abs(stress_12[plastic]), sigma_y(model, p[plastic]) / np.sqrt(3.0)
    ) <= 1e-8


@pytest.mark.parametrize("integrator_cls", PATHS_3D)
def test_p6_plastic_shear_strain_is_sqrt3_p(integrator_cls):
    model, result = simple_shear(integrator_cls)
    mu = model.E / (2.0 * (1.0 + model.nu))
    plastic_gamma_12 = result.strain[:, 3] - result.stress[:, 3] / mu
    assert rel_err(plastic_gamma_12, np.sqrt(3.0) * result.history("ep")) <= 1e-10


# ---------------------------------------------------------------------------
# P7 — without hardening, the model is J2 perfect plasticity
# ---------------------------------------------------------------------------

P7_CASES = [
    pytest.param(J2Voce3D, J2Isotropic3D, PythonAnalyticalIntegrator,
                 id="J2Voce3D hand-written"),
    pytest.param(J2Voce3D, J2Isotropic3D, PythonNumericalIntegrator,
                 id="J2Voce3D autodiff"),
    pytest.param(J2Voce1D, J2Isotropic1D, PythonNumericalIntegrator, id="J2Voce1D"),
]


@pytest.mark.parametrize("voce_cls, reference_cls, integrator_cls", P7_CASES)
def test_p7_q_zero_is_j2_isotropic_with_h_zero(voce_cls, reference_cls, integrator_cls):
    voce = voce_cls.at_defaults(Q=0.0)
    reference = reference_cls(E=voce.E, nu=voce.nu, sigma_y0=voce.sigma_y0, H=0.0)
    histories = standard_histories(voce)
    assert histories
    for name, load in histories.items():
        got = StrainDriver(integrator_cls(voce)).run(load)
        expected = StrainDriver(PythonIntegrator(reference)).run(load)
        assert plastic_steps(expected).any(), name
        assert rel_err(got.stress, expected.stress) <= 1e-8, name
        assert rel_err(got.history("ep"), expected.history("ep")) <= 1e-8, name
        assert rel_err(tangents(got), tangents(expected)) <= 1e-8, name


# ---------------------------------------------------------------------------
# P8 — the 3D model under uniaxial stress is the 1D model
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("run", UNIAXIAL_STRESS_3D)
def test_p8_uniaxial_stress_in_3d_is_the_1d_model(run):
    _, solid = run(REVERSED)
    _, bar = uniaxial_1d(REVERSED)
    assert len(solid.step_results) == len(bar.step_results) == len(REVERSED.data)
    assert rel_err(solid.stress[:, 0], bar.stress[:, 0]) <= 1e-8
    assert rel_err(solid.history("ep"), bar.history("ep")) <= 1e-8
