"""Umat — the declaration of a model's Fortran port, and what reads it.

No compiler is needed: the "compiled module" is a Python stand-in that
answers the way an f2py module would, in a throw-away package laid out like a
model family.
"""

import sys
import textwrap

import pytest

from returnmap.verification import Umat, check_model
from returnmap.verification.umat import declared_umats, umat_of

FAMILY = '''
from returnmap.models import J2Isotropic3D


class PortedJ2(J2Isotropic3D):
    """A model with a UMAT."""


class UnportedJ2(J2Isotropic3D):
    """A model of the same family without one."""
'''

VERIFICATION = '''
from returnmap.verification import Umat

from {package}.model import PortedJ2

UMATS = [
    Umat(PortedJ2, module="{module}", sources=["fortran/ported_j2.f90"],
         libraries=["lapack"], subroutine="j2"),
]
'''

# What an f2py module looks like from Python: a callable per subroutine.
# The J2 step is computed in Python; ``SCALE`` lets a test make it wrong.
STAND_IN = '''
import numpy as np
from returnmap.models import J2Isotropic3D
from returnmap.simulation import PythonIntegrator

SCALE = {scale}


def j2(E, nu, sigma_y0, H, stress_n, ep, dstran):
    model = J2Isotropic3D(E=E, nu=nu, sigma_y0=sigma_y0, H=H)
    state_n = {{"stress": np.asarray(stress_n), "ep": np.array(ep)}}
    r = PythonIntegrator(model).stress_update(np.asarray(dstran), np.asarray(stress_n), state_n)
    return SCALE * np.asarray(r.stress), float(r.state["ep"]), np.asarray(r.ddsdde)
'''

# The ABAQUS entry, as f2py exposes one declared like the shipped UMATs:
# STRESS, STATEV and PNEWDT are updated in place, DDSDDE and the energies are
# returned.  ``PROPS_ORDER`` is where it reads E, nu, sigma_y0, H from.
ENTRY = '''

PROPS_ORDER = {props_order}


def umat(stress, statev, stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname,
         ndi, nshr, props, coords, drot, pnewdt, celent, dfgrd0, dfgrd1,
         noel, npt, layer, kspt, kstep, kinc):
    E, nu, sigma_y0, H = (props[i] for i in PROPS_ORDER)
    stress_out, ep_out, ddsdde = j2(E, nu, sigma_y0, H, stress.copy(), statev[0], dstran)
    stress[:] = stress_out
    statev[0] = ep_out
    zeros = np.zeros(len(stress))
    return ddsdde, 0.0, 0.0, 0.0, 0.0, zeros, zeros, 0.0
'''


@pytest.fixture
def family(tmp_path, monkeypatch, request):
    """A package ``<name>`` with model.py and verification.py, importable."""
    package = f"fam_{request.node.name.split('[')[0][-24:].strip('_')}"
    module = f"{package}_umat"
    root = tmp_path / package
    (root / "sub").mkdir(parents=True)
    (root / "__init__.py").write_text("")
    (root / "model.py").write_text(textwrap.dedent(FAMILY))
    (root / "verification.py").write_text(
        textwrap.dedent(VERIFICATION).format(package=package, module=module))
    monkeypatch.syspath_prepend(str(tmp_path))
    yield package, module, tmp_path
    for name in [n for n in sys.modules if n == package or n.startswith(f"{package}.")
                 or n == module]:
        del sys.modules[name]


def build(tmp_path, module, scale=1.0, *, entry=True, props_order=(0, 1, 2, 3)):
    text = textwrap.dedent(STAND_IN).format(scale=scale)
    if entry:
        text += textwrap.dedent(ENTRY).format(props_order=props_order)
    (tmp_path / f"{module}.py").write_text(text)


def models_of(package):
    import importlib
    return importlib.import_module(f"{package}.model")


# ---------------------------------------------------------------------------
# The declaration
# ---------------------------------------------------------------------------

class TestDeclaration:
    def test_sources_are_relative_to_the_model_module(self, family):
        package, module, tmp_path = family
        umat = umat_of(models_of(package).PortedJ2.at_defaults())
        assert isinstance(umat, Umat)
        assert umat.module == module and list(umat.libraries) == ["lapack"]
        assert umat.source_paths() == [tmp_path / package / "fortran" / "ported_j2.f90"]

    def test_only_the_class_it_names_has_a_umat(self, family):
        package, _, _ = family
        assert umat_of(models_of(package).UnportedJ2.at_defaults()) is None

    def test_a_model_outside_a_package_has_none(self):
        from returnmap.models import J2Isotropic3D

        class Local(J2Isotropic3D):
            pass

        assert umat_of(Local.at_defaults()) is None

    def test_declared_umats_of_a_family_and_of_a_package_of_families(self, family):
        package, module, _ = family
        assert [u.module for u in declared_umats(package)] == [module]
        modules = [u.module for u in declared_umats("returnmap.models")]
        assert modules == [
            "j2_isotropic_3d", "ow_kinematic_3d", "yu_kinematic_3d", "yu_kinematic_ps"]
        assert declared_umats("returnmap.models.af_kinematic") == []

    def test_the_reference_umats_point_at_their_sources(self):
        for umat in declared_umats("returnmap.models"):
            assert all(path.is_file() for path in umat.source_paths()), umat.module

    def test_built_means_importable(self, family):
        package, module, tmp_path = family
        umat = umat_of(models_of(package).PortedJ2.at_defaults())
        assert not umat.is_built()
        build(tmp_path, module)
        assert umat.is_built()


# ---------------------------------------------------------------------------
# props: the order the UMAT takes the parameters in
# ---------------------------------------------------------------------------

REORDERED = '''
import numpy as np
from returnmap.models import J2Isotropic3D
from returnmap.simulation import PythonIntegrator


def j2(sigma_y0, H, E, stress_n, ep, dstran):          # no nu: fixed inside
    model = J2Isotropic3D(E=E, nu=0.3, sigma_y0=sigma_y0, H=H)
    state_n = {"stress": np.asarray(stress_n), "ep": np.array(ep)}
    r = PythonIntegrator(model).stress_update(np.asarray(dstran), np.asarray(stress_n), state_n)
    return np.asarray(r.stress), float(r.state["ep"]), np.asarray(r.ddsdde)
'''


class TestProps:
    def test_parameters_are_passed_in_the_declared_order(self, family):
        """A UMAT that takes (sigma_y0, H, E) and no nu."""
        package, module, tmp_path = family
        (tmp_path / f"{module}.py").write_text(textwrap.dedent(REORDERED))
        PortedJ2 = models_of(package).PortedJ2
        model = PortedJ2.at_defaults()
        umat = Umat(PortedJ2, module=module, sources=[], subroutine="j2",
                    props=["sigma_y0", "H", "E"])

        from returnmap.simulation import PythonIntegrator
        from returnmap.verification import crosscheck, standard_histories

        load = standard_histories(model)["multiaxial"]
        report = crosscheck(PythonIntegrator(model), umat.integrator(model), load)
        assert report.passed, report.summary()
        # In param_names order the same module is called with the wrong values.
        with pytest.raises(Exception):
            in_default_order = Umat(PortedJ2, module=module, sources=[], subroutine="j2")
            crosscheck(PythonIntegrator(model), in_default_order.integrator(model), load)

    def test_a_name_that_is_not_a_parameter_is_refused_at_declaration(self):
        from returnmap.models import J2Isotropic3D

        with pytest.raises(ValueError, match=r"props names \['yield_stress'\]"):
            Umat(J2Isotropic3D, module="m", sources=[], props=["E", "yield_stress"])


# ---------------------------------------------------------------------------
# check_model
# ---------------------------------------------------------------------------

class TestCheckModel:
    def test_not_declared(self, family):
        package, _, _ = family
        report = check_model(models_of(package).UnportedJ2.at_defaults())
        assert report["UMAT"].status == "not declared"
        assert report.passed

    def test_declared_but_not_built_does_not_count_against_the_model(self, family):
        package, module, _ = family
        report = check_model(models_of(package).PortedJ2.at_defaults())
        check = report["UMAT"]
        assert check.status == "not built" and check.passed and report.passed
        assert check.detail == (
            f"module {module} is not built — returnmap build --family {package}")
        assert "—    UMAT" in report.summary()

    def test_a_built_umat_is_compared_over_every_history(self, family):
        package, module, tmp_path = family
        build(tmp_path, module)
        report = check_model(models_of(package).PortedJ2.at_defaults())
        check = report["UMAT"]
        assert check.status == "passed", report.summary()
        assert set(check.reports) == {"uniaxial", "multiaxial"}
        assert check.detail.endswith(f"over 2 histories ({module})")

    def test_a_port_is_held_to_round_off_not_to_six_digits(self, family):
        """A constant typed in single precision moves the stress by 1e-7."""
        package, module, tmp_path = family
        build(tmp_path, module, scale=1.0 + 1e-7)
        report = check_model(models_of(package).PortedJ2.at_defaults())
        check = report["UMAT"]
        assert check.status == "failed", report.summary()
        assert "(tol 1e-10)" in check.detail and "(tol 1e-09)" in check.detail

    def test_a_umat_of_a_model_solved_by_autodiff_alone_is_held_to_stress_tol(
            self, af_module, monkeypatch):
        """There is no hand-written algebra for it to be a port of: its own
        iteration and the autodiff one agree as far as each was carried."""
        import returnmap.verification.model as check_module
        from returnmap.models import AFKinematic3D

        module = af_module()
        with open(__import__("importlib").util.find_spec(module).origin, "a") as stand_in:
            stand_in.write(
                "\n\n_exact = core\n\n"
                "def core(*args):\n"
                "    stress, alpha, ep, ddsdde = _exact(*args)\n"
                "    return stress * (1.0 + 1e-8), alpha, ep, ddsdde\n")
        declared = Umat(AFKinematic3D, module=module, sources=[], subroutine="core",
                        statev=["ep", "alpha"], entry=None)
        monkeypatch.setattr(check_module, "umat_of", lambda model: declared)
        report = check_model(AFKinematic3D.at_defaults())
        check = report["UMAT"]
        assert check.status == "passed", report.summary()
        assert "(tol 1e-06)" in check.detail
        assert check_model(AFKinematic3D.at_defaults(), stress_tol=1e-9)["UMAT"].status == "failed"

    def test_a_umat_that_disagrees_fails_only_its_own_row(self, family):
        package, module, tmp_path = family
        build(tmp_path, module, scale=1.001)          # stress 0.1 % too large
        report = check_model(models_of(package).PortedJ2.at_defaults())
        assert [c.name for c in report.checks if c.status == "failed"] == ["UMAT"]
        assert not report.passed
        assert "stress" in report["UMAT"].detail


# ---------------------------------------------------------------------------
# The ABAQUS entry
# ---------------------------------------------------------------------------

class TestEntryRow:
    def test_not_declared_without_a_umat(self, family):
        package, _, _ = family
        report = check_model(models_of(package).UnportedJ2.at_defaults())
        assert report["UMAT entry"].status == "not declared"

    def test_not_built_does_not_count_against_the_model(self, family):
        package, _, _ = family
        report = check_model(models_of(package).PortedJ2.at_defaults())
        assert report["UMAT entry"].status == "not built" and report.passed

    def test_a_faithful_entry_agrees_with_the_core(self, family):
        package, module, tmp_path = family
        build(tmp_path, module)
        report = check_model(models_of(package).PortedJ2.at_defaults())
        check = report["UMAT entry"]
        assert check.status == "passed", report.summary()
        assert set(check.reports) == {"uniaxial", "multiaxial"}
        assert check.detail.endswith(f"against j2 over 2 histories ({module}.umat)")
        assert report.passed

    def test_an_entry_that_reads_props_in_another_order_fails_only_its_own_row(self, family):
        package, module, tmp_path = family
        build(tmp_path, module, props_order=(0, 1, 3, 2))      # sigma_y0 and H swapped
        report = check_model(models_of(package).PortedJ2.at_defaults())
        assert [c.name for c in report.checks if c.status == "failed"] == ["UMAT entry"]

    def test_the_entry_is_held_to_round_off_of_the_core(self, family):
        """An entry that returns the core's stress to seven digits only."""
        package, module, tmp_path = family
        build(tmp_path, module)
        with open(tmp_path / f"{module}.py", "a") as stand_in:
            stand_in.write(
                "\n\n_exact = umat\n\n"
                "def umat(stress, *args):\n"
                "    out = _exact(stress, *args)\n"
                "    stress *= 1.0 + 1e-7\n"
                "    return out\n")
        report = check_model(models_of(package).PortedJ2.at_defaults())
        check = report["UMAT entry"]
        assert check.status == "failed"
        assert "(tol 1e-12)" in check.detail

    def test_the_entry_is_compared_with_the_core_not_with_the_model(self, family):
        """A core that is wrong fails its own row; the entry that wraps it faithfully passes."""
        package, module, tmp_path = family
        build(tmp_path, module, scale=1.001)
        report = check_model(models_of(package).PortedJ2.at_defaults())
        assert report["UMAT"].status == "failed"
        assert report["UMAT entry"].status == "passed"

    def test_a_built_module_without_the_entry_fails(self, family):
        package, module, tmp_path = family
        build(tmp_path, module, entry=False)
        report = check_model(models_of(package).PortedJ2.at_defaults())
        check = report["UMAT entry"]
        assert check.status == "failed"
        assert "has no subroutine 'umat'" in check.detail and "entry=None" in check.detail

    def test_no_entry_declared_is_not_checked(self, family, monkeypatch):
        package, module, tmp_path = family
        build(tmp_path, module, entry=False)
        PortedJ2 = models_of(package).PortedJ2
        declared = Umat(PortedJ2, module=module, sources=[], subroutine="j2", entry=None)
        import returnmap.verification.model as check_module

        monkeypatch.setattr(check_module, "umat_of", lambda model: declared)
        report = check_model(PortedJ2.at_defaults())
        assert report["UMAT entry"].status == "not declared"
        assert report["UMAT"].status == "passed" and report.passed

    def test_an_entry_with_another_argument_list_fails_with_the_reason(self, family):
        package, module, tmp_path = family
        build(tmp_path, module, entry=False)
        with open(tmp_path / f"{module}.py", "a") as stand_in:
            stand_in.write("\n\ndef umat(stress, statev, props):\n    return None\n")
        report = check_model(models_of(package).PortedJ2.at_defaults())
        check = report["UMAT entry"]
        assert check.status == "failed"
        assert "could not be called with the ABAQUS argument list" in check.detail


AF_STAND_IN = '''
import numpy as np
from returnmap.models import AFKinematic3D
from returnmap.simulation import PythonIntegrator

GIVE_UP_ABOVE = {give_up_above}


def core(E, nu, sigma_y0, C_k, gamma, stress_n, alpha, ep, dstran):
    model = AFKinematic3D(E=E, nu=nu, sigma_y0=sigma_y0, C_k=C_k, gamma=gamma)
    state_n = {{"stress": np.asarray(stress_n), "alpha": np.asarray(alpha), "ep": np.array(ep)}}
    r = PythonIntegrator(model).stress_update(np.asarray(dstran), np.asarray(stress_n), state_n)
    return np.asarray(r.stress), np.asarray(r.state["alpha"]), float(r.state["ep"]), np.asarray(r.ddsdde)


def umat(stress, statev, stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname,
         ndi, nshr, props, coords, drot, pnewdt, celent, dfgrd0, dfgrd1,
         noel, npt, layer, kspt, kstep, kinc):
    """STATEV: 1 = ep, 2-7 = alpha."""
    zeros = np.zeros(len(stress))
    stress_out, alpha, ep, ddsdde = core(*props, stress.copy(), statev[1:7].copy(), statev[0], dstran)
    if np.max(np.abs(dstran)) > GIVE_UP_ABOVE:
        pnewdt[...] = 0.5                       # STRESS and STATEV left as they were
        return ddsdde, 0.0, 0.0, 0.0, 0.0, zeros, zeros, 0.0
    stress[:] = stress_out
    statev[0] = ep
    statev[1:7] = alpha
    return ddsdde, 0.0, 0.0, 0.0, 0.0, zeros, zeros, 0.0
'''


@pytest.fixture
def af_module(tmp_path, monkeypatch, request):
    """An importable stand-in for a compiled AF UMAT; ``build(give_up_above)`` writes it."""
    module = f"af_{request.node.name.split('[')[0][-24:].strip('_')}"
    monkeypatch.syspath_prepend(str(tmp_path))

    def build(give_up_above=1.0):
        (tmp_path / f"{module}.py").write_text(
            textwrap.dedent(AF_STAND_IN).format(give_up_above=give_up_above))
        return module

    yield build
    sys.modules.pop(module, None)


class TestStatev:
    """Where the entry keeps the state variables, when not in declaration order."""

    @staticmethod
    def run(umat, model):
        from returnmap.verification import crosscheck, standard_histories

        load = standard_histories(model)["multiaxial"]
        core = umat.integrator(model)
        return crosscheck(core, umat.entry_integrator(model), load, lockstep=True)

    def test_the_declared_layout_is_used_for_packing_and_unpacking(self, af_module):
        from returnmap.models import AFKinematic3D

        model = AFKinematic3D.at_defaults()
        umat = Umat(AFKinematic3D, module=af_module(), sources=[], subroutine="core",
                    statev=["ep", "alpha"])
        report = self.run(umat, model)
        assert report.passed, report.summary()
        assert report.max_stress_rel_err == 0.0

    def test_the_declaration_order_is_the_default_and_is_wrong_for_this_module(self, af_module):
        from returnmap.models import AFKinematic3D

        model = AFKinematic3D.at_defaults()
        umat = Umat(AFKinematic3D, module=af_module(), sources=[], subroutine="core")
        assert not self.run(umat, model).passed

    def test_statev_has_to_hold_each_state_variable_once(self):
        from returnmap.models import AFKinematic3D

        with pytest.raises(ValueError, match=r"statev names \['ep'\]"):
            Umat(AFKinematic3D, module="m", sources=[], statev=["ep"])
        with pytest.raises(ValueError, match="each state variable of the model once"):
            Umat(AFKinematic3D, module="m", sources=[], statev=["ep", "alpha", "stress"])

    def test_a_lowered_pnewdt_is_a_step_that_did_not_converge(self, af_module):
        import numpy as np
        from returnmap.models import AFKinematic3D

        model = AFKinematic3D.at_defaults()
        umat = Umat(AFKinematic3D, module=af_module(give_up_above=1e-3), sources=[],
                    subroutine="core", statev=["ep", "alpha"])
        entry = umat.entry_integrator(model, raise_on_nonconverged=False)
        stress_n, state_n = np.zeros(6), model.initial_state()
        small = entry.stress_update(np.array([5e-4, 0, 0, 0, 0, 0.0]), stress_n, state_n)
        large = entry.stress_update(np.array([5e-3, 0, 0, 0, 0, 0.0]), stress_n, state_n)
        assert small.converged and not large.converged
        np.testing.assert_array_equal(large.stress, stress_n)        # handed back untouched
        with pytest.raises(RuntimeError, match="PNEWDT = 0.5"):
            umat.entry_integrator(model).stress_update(
                np.array([5e-3, 0, 0, 0, 0, 0.0]), stress_n, state_n)

    def test_an_entry_that_is_not_declared_cannot_be_asked_for(self, af_module):
        from returnmap.models import AFKinematic3D

        umat = Umat(AFKinematic3D, module=af_module(), sources=[], subroutine="core", entry=None)
        assert not umat.has_entry()
        with pytest.raises(AttributeError, match="declares no ABAQUS entry"):
            umat.entry_integrator(AFKinematic3D.at_defaults())
