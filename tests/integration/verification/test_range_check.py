"""check_range — a model at the ends of its parameter ranges and under random loading."""

import numpy as np
import pytest

from returnmap.core import ParamMeta, ReturnMappingResult
from returnmap.models import AFKinematic3D, J2Isotropic1D, J2Isotropic3D
from returnmap.simulation import FieldType
from returnmap.verification import (
    RangeCase, RangeReport, check_model, check_range, parameter_variants,
    random_histories, standard_histories,
)


class _WrongOnLargeSteps(J2Isotropic3D):
    """Closed form with the plastic correction 0.1 % too large — but only on a
    step whose trial stress is more than ten times the yield stress, which no
    standard history takes."""

    def user_defined_return_mapping(self, inc):
        rm = super().user_defined_return_mapping(inc)
        if self.vonmises(inc.trial["stress"]) < 10.0 * self.sigma_y0:
            return rm
        stress = inc.trial["stress"] + 1.001 * (rm.stress - inc.trial["stress"])
        return ReturnMappingResult(
            stress=stress, state={**rm.state, "stress": stress},
            dlambda=rm.dlambda, n_iterations=1, residual_history=[1.0, 0.0],
            converged=True,
        )


class _WrongWithoutHardening(J2Isotropic3D):
    """Hand-written tangent 1 % too stiff when ``H`` is exactly zero — the lower
    end of the documented range, not the typical value."""

    def user_defined_tangent(self, new, dlambda, inc):
        scale = 1.01 if self.H == 0.0 else 1.0
        return scale * super().user_defined_tangent(new, dlambda, inc)


class _GivesUpOnLargeSteps(J2Isotropic3D):
    """Right wherever it answers; reports "not converged" on a step whose trial
    stress is more than ten times the yield stress."""

    def user_defined_return_mapping(self, inc):
        rm = super().user_defined_return_mapping(inc)
        if self.vonmises(inc.trial["stress"]) < 10.0 * self.sigma_y0:
            return rm
        return ReturnMappingResult(
            stress=rm.stress, state=rm.state, dlambda=rm.dlambda, n_iterations=1,
            residual_history=[1.0, 1.0], converged=False,
        )


class _NoBounds(J2Isotropic3D):
    param_meta = {
        name: ParamMeta(doc=meta.doc, unit=meta.unit, default=meta.default)
        for name, meta in J2Isotropic3D.param_meta.items()
    }


# ---------------------------------------------------------------------------
# The histories and the variants
# ---------------------------------------------------------------------------

class TestRandomHistories:
    def test_three_step_sizes_in_every_component(self):
        model = J2Isotropic3D.at_defaults()
        histories = random_histories(model)
        assert list(histories) == ["random", "random large", "random single-step"]
        assert [len(h.data) for h in histories.values()] == [36, 24, 6]
        for load in histories.values():
            assert load.type is FieldType.STRAIN
            assert load.data.shape[1] == model.ntens
            assert np.all(np.abs(load.data).max(axis=0) > 0.0)        # no idle component

    def test_default_amplitude_is_that_of_the_histories_check_model_uses(self):
        model = J2Isotropic3D.at_defaults()
        amplitude = max(np.abs(h.data).max() for h in standard_histories(model).values())
        explicit = random_histories(model, amplitude)
        for name, load in random_histories(model).items():
            np.testing.assert_array_equal(load.data, explicit[name].data)

    def test_the_same_seed_gives_the_same_histories(self):
        model = AFKinematic3D.at_defaults()
        first, again, other = (random_histories(model, seed=s) for s in (3, 3, 4))
        for name in first:
            np.testing.assert_array_equal(first[name].data, again[name].data)
            assert not np.array_equal(first[name].data, other[name].data)

    def test_step_sizes_grow_from_one_history_to_the_next(self):
        steps = [
            np.abs(np.diff(load.data, axis=0, prepend=0.0)).max()
            for load in random_histories(J2Isotropic3D.at_defaults()).values()
        ]
        assert steps[0] < steps[1] < steps[2]

    def test_one_component_model(self):
        histories = random_histories(J2Isotropic1D.at_defaults())
        assert all(load.data.shape[1] == 1 for load in histories.values())

    def test_amplitude_must_be_positive(self):
        with pytest.raises(ValueError, match="amplitude must be positive"):
            random_histories(J2Isotropic3D.at_defaults(), 0.0)


class TestParameterVariants:
    def test_each_end_of_each_documented_range(self):
        model = J2Isotropic3D.at_defaults()
        variants = parameter_variants(model)
        assert list(variants) == [
            "E=1000", "E=1e+06", "nu=0", "nu=0.49", "sigma_y0=1", "sigma_y0=5000",
            "H=0", "H=100000",
        ]
        assert variants["H=0"].H == 0.0
        assert variants["H=0"].E == model.E                 # one parameter at a time
        assert type(variants["H=0"]) is J2Isotropic3D

    def test_a_bound_the_model_already_sits_at_gives_no_variant(self):
        variants = parameter_variants(J2Isotropic3D.at_defaults(H=0.0))
        assert "H=0" not in variants and "H=100000" in variants

    def test_no_documented_range_gives_no_variant(self):
        assert parameter_variants(_NoBounds.at_defaults()) == {}


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------

class TestCheckRange:
    def test_a_sound_model_passes_every_case(self):
        report = check_range(J2Isotropic3D.at_defaults())
        assert isinstance(report, RangeReport)
        assert report.passed, report.summary()
        assert [c.group for c in report.cases] == ["parameters"] * 8 + ["loading"] * 3
        assert all(isinstance(case, RangeCase) for case in report.cases)
        assert str(report) == report.summary()
        assert report.summary().startswith("J2Isotropic3D over its range: every case passes")

    def test_a_parameter_case_is_check_model_at_the_variant(self):
        model = J2Isotropic3D.at_defaults()
        report = check_range(model, histories={})
        case = report["H=0"]
        assert case.group == "parameters"
        own = check_model(model.with_params(H=0.0))
        assert [c.status for c in case.reports["model"].checks] == [c.status for c in own.checks]

    def test_a_defect_at_the_end_of_a_range_fails_that_variant_only(self):
        model = _WrongWithoutHardening.at_defaults()
        assert check_model(model).passed                    # invisible at the typical values
        report = check_range(model, histories={})
        assert not report.passed
        assert [c.name for c in report.cases if c.status == "failed"] == ["H=0"]
        assert "hand-written tangent" in report["H=0"].detail
        assert "1 case(s) failed — H=0" in report.summary()

    def test_a_wrong_answer_on_a_large_step_fails_the_loading_case(self):
        model = _WrongOnLargeSteps.at_defaults()
        assert check_model(model).passed                    # the standard histories do not reach it
        report = check_range(model, variants={})
        assert report["random"].status == "passed"
        for name in ("random large", "random single-step"):
            case = report[name]
            assert case.status == "failed", report.summary()
            assert "solution" in case.detail
            assert not case.reports["solution"].passed
            assert not all(
                s.passed for s in case.reports["hand-written vs autodiff"].steps
            )

    def test_an_increment_that_is_not_solved_is_counted_not_failed(self):
        report = check_range(_GivesUpOnLargeSteps.at_defaults(), variants={})
        assert report.passed, report.summary()
        large = report["random large"]
        assert large.status == "passed"
        solved, total = (int(n) for n in large.detail.split(" increments")[0].split(" of "))
        assert 0 < solved < total
        single = report["random single-step"]
        assert single.status == "not solved"
        assert "nothing to compare" in single.detail

    def test_a_model_solved_by_autodiff_alone(self):
        report = check_range(AFKinematic3D.at_defaults(), variants={})
        assert report.passed, report.summary()
        assert set(report["random"].reports) == {"tangent"}

    def test_no_documented_range_is_said_and_does_not_fail(self):
        report = check_range(_NoBounds.at_defaults())
        assert report.passed
        assert [c.group for c in report.cases] == ["loading"] * 3
        assert any("documents no parameter range" in note for note in report.notes)
        assert "note:" in report.summary()

    def test_nothing_to_run_is_not_a_pass(self):
        report = check_range(J2Isotropic3D.at_defaults(), variants={}, histories={})
        assert not report.passed
        assert "no case was run" in report.summary()

    def test_given_variants_and_histories_are_the_only_ones_used(self):
        model = J2Isotropic3D.at_defaults()
        report = check_range(
            model,
            variants={"soft": model.with_params(E=50_000.0)},
            histories={"one": random_histories(model, seed=7)["random"]},
        )
        assert [c.name for c in report.cases] == ["soft", "one"]
        assert report.passed, report.summary()

    def test_variants_a_family_declares_replace_the_default_ones(self, monkeypatch):
        import returnmap.verification.range_check as module

        model = J2Isotropic3D.at_defaults()
        declared = {"stiff": model.with_params(H=50_000.0)}
        monkeypatch.setattr(module, "_declared_variants", lambda model: declared)
        report = check_range(model, histories={})
        assert [c.name for c in report.cases] == ["stiff"]
        assert report.notes == []

    def test_no_declaration_means_the_documented_ranges(self):
        import returnmap.verification.range_check as module

        assert module._declared_variants(J2Isotropic3D.at_defaults()) is None

    def test_an_unknown_case_name_lists_the_cases(self):
        report = check_range(J2Isotropic3D.at_defaults(), variants={})
        with pytest.raises(KeyError, match="random large"):
            report["H=0"]

    def test_a_class_is_not_a_model(self):
        with pytest.raises(TypeError, match="at_defaults"):
            check_range(J2Isotropic3D)

    def test_a_stress_history_is_refused(self):
        from returnmap.simulation import FieldHistory

        model = J2Isotropic3D.at_defaults()
        load = FieldHistory(FieldType.STRESS, np.ones((3, 6)))
        with pytest.raises(ValueError, match="not a strain history"):
            check_range(model, variants={}, histories={"stress": load})


# ---------------------------------------------------------------------------
# The UMAT, under random loading
# ---------------------------------------------------------------------------

class _StandIn:
    """What ``umat_of`` returns, with a Python integrator in place of the
    compiled module: the comparison is the same, and needs no compiler."""

    def __init__(self, port, *, reads_its_flag=True, entry=None):
        self._port = port
        self._entry = entry
        self.integrator_options = {"parse_umat_failure": None} if reads_its_flag else {}

    def is_built(self):
        return True

    def has_entry(self):
        return self._entry is not None

    def integrator(self, model, **options):
        from returnmap.simulation import PythonIntegrator

        return PythonIntegrator(self._port(model), **options)

    def entry_integrator(self, model, **options):
        from returnmap.simulation import PythonIntegrator

        return PythonIntegrator(self._entry(model), **options)


class TestUmatUnderRandomLoading:
    @staticmethod
    def declare(monkeypatch, stand_in):
        import returnmap.verification.range_check as module

        monkeypatch.setattr(module, "umat_of", lambda model: stand_in)

    def test_a_faithful_port_agrees_on_every_history(self, monkeypatch):
        self.declare(monkeypatch, _StandIn(lambda model: model))
        report = check_range(J2Isotropic3D.at_defaults(), variants={})
        assert report.passed, report.summary()
        for case in report.cases:
            assert case.reports["UMAT"].passed
            assert "UMAT 0.0e+00" in case.detail

    def test_a_port_with_another_constant_fails(self, monkeypatch):
        self.declare(monkeypatch, _StandIn(lambda model: model.with_params(H=1.01 * model.H)))
        report = check_range(J2Isotropic3D.at_defaults(), variants={})
        assert [c.status for c in report.cases] == ["failed"] * 3
        assert all(case.detail.startswith("UMAT failed") for case in report.cases)

    def test_a_port_that_converges_where_the_model_gives_up_fails(self, monkeypatch):
        sound = J2Isotropic3D.at_defaults()
        self.declare(monkeypatch, _StandIn(lambda model: sound))
        report = check_range(_GivesUpOnLargeSteps.at_defaults(), variants={})
        large = report["random large"]
        assert large.status == "failed"
        assert "the Python model did not solve" in large.detail
        assert report["random"].status == "passed"          # nothing was given up on there

    def test_without_a_flag_to_read_the_given_up_increments_are_not_compared(self, monkeypatch):
        sound = J2Isotropic3D.at_defaults()
        self.declare(monkeypatch, _StandIn(lambda model: sound, reads_its_flag=False))
        report = check_range(_GivesUpOnLargeSteps.at_defaults(), variants={})
        assert report["random large"].status == "passed"

    def test_a_umat_that_is_not_built_is_said_and_does_not_fail(self, monkeypatch):
        stand_in = _StandIn(lambda model: model)
        stand_in.is_built = lambda: False
        self.declare(monkeypatch, stand_in)
        report = check_range(J2Isotropic3D.at_defaults(), variants={})
        assert report.passed
        assert all("UMAT not built" in case.detail for case in report.cases)
        assert all("UMAT" not in case.reports for case in report.cases)


class TestEntryUnderRandomLoading:
    declare = staticmethod(TestUmatUnderRandomLoading.declare)

    def test_a_faithful_entry_agrees_with_the_core(self, monkeypatch):
        self.declare(monkeypatch, _StandIn(lambda model: model, entry=lambda model: model))
        report = check_range(J2Isotropic3D.at_defaults(), variants={})
        assert report.passed, report.summary()
        assert all(case.reports["UMAT entry"].passed for case in report.cases)

    def test_no_entry_means_no_comparison(self, monkeypatch):
        self.declare(monkeypatch, _StandIn(lambda model: model))
        report = check_range(J2Isotropic3D.at_defaults(), variants={})
        assert all("UMAT entry" not in case.reports for case in report.cases)

    def test_an_entry_that_unpacks_another_constant_fails_its_own_part(self, monkeypatch):
        self.declare(monkeypatch, _StandIn(
            lambda model: model, entry=lambda model: model.with_params(H=1.01 * model.H)))
        report = check_range(J2Isotropic3D.at_defaults(), variants={})
        for case in report.cases:
            assert case.status == "failed"
            assert case.detail.startswith("UMAT entry failed")
            assert case.reports["UMAT"].passed

    def test_an_entry_that_moves_on_where_the_model_gives_up_fails(self, monkeypatch):
        sound = J2Isotropic3D.at_defaults()
        self.declare(monkeypatch, _StandIn(
            lambda model: model, reads_its_flag=False, entry=lambda model: sound))
        report = check_range(_GivesUpOnLargeSteps.at_defaults(), variants={})
        large = report["random large"]
        assert large.status == "failed"
        assert "converged or changed STRESS / STATEV" in large.detail


@pytest.mark.fortran
def test_the_compiled_umat_is_compared_in_every_case():
    """With ``j2_isotropic_3d`` built: every variant and every random history."""
    pytest.importorskip("j2_isotropic_3d", reason="j2_isotropic_3d is not built")
    report = check_range(J2Isotropic3D.at_defaults())
    assert report.passed, report.summary()
    for case in report.cases:
        if case.group == "parameters":
            assert case.reports["model"]["UMAT"].status == "passed", case.name
            assert case.reports["model"]["UMAT entry"].status == "passed", case.name
        else:
            assert case.reports["UMAT"].passed, case.name
            assert case.reports["UMAT entry"].passed, case.name


def test_the_tangent_is_differenced_with_a_smaller_step_where_strains_are_small():
    """A yield strain two hundred times smaller is checked as closely: the
    difference step shrinks with it instead of staying at 1e-7."""
    typical = check_model(J2Isotropic3D.at_defaults())
    small = check_model(J2Isotropic3D.at_defaults(sigma_y0=1.0))
    assert small.passed, small.summary()
    errors = [
        max(r.max_rel_err for r in report["tangent"].reports.values())
        for report in (typical, small)
    ]
    assert errors[1] < 100.0 * errors[0]
