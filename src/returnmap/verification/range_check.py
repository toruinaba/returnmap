"""A model over the range it documents, not only at one point of it.

:func:`~returnmap.verification.check_model` checks a model at one set of
parameters, on loading histories sized to its yield strain.  A defect that
only shows elsewhere — a branch of a hand-written solver the standard
histories never reach, two ``PROPS`` swapped in a UMAT whose typical values
happen to be equal, a mechanism that is idle at the typical parameters —
passes it.  :func:`check_range` runs the same comparisons in two more
directions:

==============  ===========================================================
parameters      the model at each end of the range every parameter
                documents (``param_meta`` bounds), one parameter at a time,
                each through :func:`~returnmap.verification.check_model`
loading         random strain histories in every component, at three step
                sizes up to a hundred times the usual one
==============  ===========================================================

::

    from returnmap.verification import check_range

    report = check_range(model)
    assert report.passed, report.summary()

It costs one ``check_model`` per parameter variant, so it is a check to run
when a model is finished rather than in every test run.

Parameters that are not independent
-----------------------------------
The bounds of ``param_meta`` are given parameter by parameter, and moving one
of them to its bound while the others keep their typical values can leave the
model without meaning — a bounding surface inside the yield surface, say.
A family whose parameters are tied this way declares the variants worth
checking, as ``variants(model)`` in the module named ``verification`` next to
the model (where ``histories(model)`` is declared too)::

    def variants(model):
        return {
            name: variant for name, variant in parameter_variants(model).items()
            if variant.B > variant.Y
        }

:func:`check_range` uses the declared variants in place of
:func:`parameter_variants` when there are any.

What counts as a failure under random loading
----------------------------------------------
Large random increments are harder than anything a model is normally asked
to solve, and a solver that gives up on one is not thereby wrong.  A loading
case therefore fails only on a *wrong answer*:

* a hand-written solution that is not a root of the framework residual
  (:func:`~returnmap.verification.check_solution`);
* a step on which the hand-written and the autodiff return mapping both
  converge, from the same state, to different results;
* a tangent that disagrees with central differences;
* a compiled UMAT that, given the same state and increment, returns something
  else than the Python model — or converges where the Python model does not.

Increments the model's own path (the hand-written return mapping when there
is one, the autodiff one otherwise) does not solve are left out, the history
continues from the last solved state, and the report says how many there
were.  Steps the autodiff path does not solve while the hand-written one
does are counted the same way: the hand-written solution is still checked
against the residual, which needs no Newton iteration to converge.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np

from returnmap.core.material import MaterialModel
from returnmap.core.material.hand_block import hand_tangent_record
from returnmap.core.result import ReturnMappingResult, StressUpdateResult
from returnmap.simulation.integrator import (
    PythonAnalyticalIntegrator, PythonIntegrator, PythonNumericalIntegrator,
)
from returnmap.simulation.types import FieldHistory, FieldType
from returnmap.verification.crosscheck import crosscheck
from returnmap.verification.jacobian import check_hand_blocks
from returnmap.verification.model import (
    ModelReport, _check_tangents, _family_declaration, _model_histories, _worst,
    check_model,
)
from returnmap.verification.solution import check_solution
from returnmap.verification.umat import umat_of

if TYPE_CHECKING:
    from returnmap.simulation.integrator import StressIntegrator

RangeStatus = Literal["passed", "failed", "not solved"]

# (name, size relative to the amplitude, steps per leg): the first is of the
# size of the standard histories, the last takes each leg in one increment.
_RANDOM_SHAPES = (("random", 1.0, 6), ("random large", 4.0, 4), ("random single-step", 20.0, 1))
_N_LEGS = 6
_SOLVER_ERRORS = (np.linalg.LinAlgError, FloatingPointError, ZeroDivisionError)


@dataclass
class RangeCase:
    """One line of a :class:`RangeReport`.

    Attributes
    ----------
    name : str
        The parameter variant (``"gamma=0"``) or the loading history
        (``"random large"``).
    group : {"parameters", "loading"}
    status : {"passed", "failed", "not solved"}
        ``"not solved"`` is a loading history on which the model's own path
        solved no plastic increment, so that nothing was compared; it does
        not count against the model.
    detail : str
        What was compared and how far apart the two sides were, or what
        failed.
    reports : dict
        The reports behind the line: ``{"model": ModelReport}`` for a
        parameter variant; for a loading history the underlying checks by
        name (``"solution"``, ``"hand-written vs autodiff"``, ``"tangent"``,
        ``"hand blocks"``, ``"UMAT"``).
    """

    name: str
    group: Literal["parameters", "loading"]
    status: RangeStatus
    detail: str
    reports: dict = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return self.status != "failed"


@dataclass
class RangeReport:
    """Result of :func:`check_range`.

    Attributes
    ----------
    model : str
        Class name of the model that was checked.
    cases : list[RangeCase]
        One entry per parameter variant, then one per loading history.
    notes : list[str]
        What was left out and why — a model that documents no parameter
        range, for instance.
    """

    model: str
    cases: list[RangeCase]
    notes: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        """At least one case was run and none failed."""
        return bool(self.cases) and all(case.passed for case in self.cases)

    def __getitem__(self, name: str) -> RangeCase:
        for case in self.cases:
            if case.name == name:
                return case
        raise KeyError(f"no case named {name!r}; have {[c.name for c in self.cases]}")

    def summary(self) -> str:
        """A table: one row per case, parameter variants first."""
        failed = [c.name for c in self.cases if c.status == "failed"]
        if not self.cases:
            head = f"{self.model} over its range: no case was run"
        elif failed:
            head = (f"{self.model} over its range: {len(failed)} case(s) failed — "
                    f"{', '.join(failed)}")
        else:
            head = f"{self.model} over its range: every case passes"
        marks = {"passed": "ok  ", "failed": "FAIL", "not solved": "—   "}
        titles = {
            "parameters": "parameters — check_model at each end of every documented range",
            "loading": "loading — random strain histories",
        }
        width = max((len(c.name) for c in self.cases), default=0)
        lines = [head]
        for group, title in titles.items():
            cases = [c for c in self.cases if c.group == group]
            if not cases:
                continue
            lines.append(f"  {title}")
            lines += [
                f"    {marks[c.status]} {c.name.ljust(width)}  {c.detail}" for c in cases
            ]
        lines += [f"  note: {note}" for note in self.notes]
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()


def random_histories(
    model: MaterialModel, amplitude: float | None = None, *, seed: int = 0,
) -> dict[str, FieldHistory]:
    """Random strain histories, at three step sizes.

    Each history is a walk through six random points of strain space — every
    component is driven, in no fixed proportion — taken in straight legs:

    ``"random"``
        Points of the size of ``amplitude``, six increments a leg: the size
        of step the standard histories take, in directions they do not.
    ``"random large"``
        Four times as far, four increments a leg.
    ``"random single-step"``
        Twenty times as far, each leg in one increment.

    Parameters
    ----------
    model : MaterialModel
    amplitude : float or None
        Standard deviation of the strain components of the first history's
        points.  ``None`` (default): the largest strain component of the
        histories :func:`~returnmap.verification.check_model` uses for the
        model — the standard ones and those its family declares — so that a
        family which rescales the standard histories rescales these with
        them.
    seed : int, default 0
        Seed of the random directions.  The histories are the same for the
        same model, amplitude and seed.

    Returns
    -------
    dict[str, FieldHistory]
        Strain histories by name, each of shape ``(N, ntens)``.
    """
    if amplitude is None:
        amplitude = max(
            float(np.max(np.abs(load.data))) for load in _model_histories(model).values()
        )
    if not amplitude > 0.0:
        raise ValueError(f"amplitude must be positive, got {amplitude!r}")
    rng = np.random.default_rng(seed)
    ntens = model.ntens
    histories: dict[str, FieldHistory] = {}
    for name, size, n_per_leg in _RANDOM_SHAPES:
        corners = [np.zeros(ntens)] + [
            size * amplitude * rng.normal(size=ntens) for _ in range(_N_LEGS)
        ]
        rows = [
            a + (b - a) * t
            for a, b in zip(corners, corners[1:])
            for t in np.linspace(0.0, 1.0, n_per_leg + 1)[1:]
        ]
        histories[name] = FieldHistory(FieldType.STRAIN, np.array(rows))
    return histories


def parameter_variants(model: MaterialModel) -> dict[str, MaterialModel]:
    """The model at each end of the range every parameter documents.

    One parameter at a time is moved to the lower and to the upper bound of
    its ``param_meta`` entry, the others keeping the values ``model`` has.
    A bound the parameter already sits at gives no variant.

    Returns
    -------
    dict[str, MaterialModel]
        Models by name (``"gamma=0"``, ``"gamma=10000"``), in the order of
        ``param_names``.  Empty for a model that documents no bounds.
    """
    variants: dict[str, MaterialModel] = {}
    meta = getattr(type(model), "param_meta", None) or {}
    for name in model.param_names:
        entry = meta.get(name)
        if entry is None or entry.bounds is None:
            continue
        for bound in entry.bounds:
            if bound == getattr(model, name):
                continue
            variants[f"{name}={bound:g}"] = model.with_params(**{name: bound})
    return variants


class _GivingUp:
    """An integrator that answers a step it cannot solve with "not converged".

    The candidate of a lockstep :func:`crosscheck` is called once per step
    and is expected to return a result.  The autodiff return mapping, asked
    for an increment far beyond what it normally solves, may instead raise
    from its linear solve, or return an iterate whose tangent cannot be
    formed.  Both become a result that carries no values and no tangent, and
    says it did not converge.
    """

    def __init__(self, integrator: "StressIntegrator") -> None:
        self._integrator = integrator

    def __getattr__(self, name: str):
        return getattr(self._integrator, name)

    def stress_update(self, strain_inc, stress_n, state_n) -> StressUpdateResult:
        try:
            result = self._integrator.stress_update(strain_inc, stress_n, state_n)
        except _SOLVER_ERRORS:
            result = None
        if result is not None and result.converged:
            return result

        def blank(value) -> np.ndarray:
            return np.full(np.shape(value), np.nan)

        unsolved = ReturnMappingResult(
            stress=blank(stress_n),
            state={key: blank(value) for key, value in dict(state_n.items()).items()},
            dlambda=np.array(np.nan), n_iterations=0, residual_history=[], converged=False,
        )
        return StressUpdateResult(
            return_mapping=unsolved, ddsdde=None, is_plastic=None, increment=None,
        )


def _declared_variants(model: MaterialModel) -> dict[str, MaterialModel] | None:
    """Variants the model's family declares in its ``verification`` module, if any."""
    declare = _family_declaration(model, "variants")
    return dict(declare(model)) if declare is not None else None


def _solved_part(
    reference: "StressIntegrator", load: FieldHistory,
) -> tuple[FieldHistory | None, list["StressUpdateResult"], list[tuple]]:
    """Drive ``reference`` through ``load``, leaving out what it does not solve.

    Returns the history made of the increments it solved (``None`` if it
    solved none), their results, and the ``(strain_inc, stress_n, state_n)``
    of every increment it did not solve.
    """
    data = np.asarray(load.data, dtype=float)
    data = data.reshape(len(data), -1)
    increments = np.diff(data, axis=0, prepend=np.zeros((1, data.shape[1])))
    stress = np.zeros(data.shape[1])
    state = reference.initial_state()
    solved_rows: list[np.ndarray] = []
    results: list["StressUpdateResult"] = []
    unsolved: list[tuple] = []
    strain = np.zeros(data.shape[1])
    for strain_inc in increments:
        try:
            result = reference.stress_update(strain_inc, stress, state)
            ok = bool(result.converged) and bool(
                np.all(np.isfinite(np.asarray(result.stress, dtype=float)))
            )
        except _SOLVER_ERRORS:
            ok = False
        if not ok:
            unsolved.append((strain_inc, stress, state))
            continue
        strain = strain + strain_inc
        solved_rows.append(strain)
        results.append(result)
        stress, state = result.stress, result.state
    if not solved_rows:
        return None, results, unsolved
    return FieldHistory(FieldType.STRAIN, np.array(solved_rows)), results, unsolved


def _check_loading(
    model: MaterialModel, name: str, load: FieldHistory, *,
    tangent_tol: float, block_tol: float, solution_tol: float, stress_tol: float,
    n_tangent_points: int,
) -> RangeCase:
    """One random history: every comparison that applies, wrong answers only."""
    if load.type is not FieldType.STRAIN:
        raise ValueError(f"history {name!r} is not a strain history")
    cls = type(model)
    hand_rm = cls.user_defined_return_mapping is not MaterialModel.user_defined_return_mapping
    reference = PythonIntegrator(model, raise_on_nonconverged=False)
    solved, results, unsolved = _solved_part(reference, load)
    plastic = [r for r in results if r.is_plastic]
    n_total = len(results) + len(unsolved)
    coverage = f"{len(results)} of {n_total} increments solved, {len(plastic)} plastic"
    if solved is None or not plastic:
        return RangeCase(name, "loading", "not solved", f"{coverage} — nothing to compare")

    parts: list[str] = []
    failures: list[str] = []
    reports: dict = {}

    def record(label: str, ok: bool, text: str, report) -> None:
        reports[label] = report
        parts.append(f"{label} {text}")
        if not ok:
            failures.append(label)

    if hand_rm:
        report = check_solution(model, plastic, tol=solution_tol)
        worst = _worst([row.worst_err for row in report.rows.values()])
        record("solution", report.passed, f"{worst:.1e}", report)

        # From the same state, with the same increment: the two return
        # mappings.  The tangent is compared with differences below, where a
        # declared approximation is taken into account.
        numerical = _GivingUp(PythonNumericalIntegrator(model, raise_on_nonconverged=False))
        report = crosscheck(
            reference, numerical, solved, lockstep=True,
            stress_tol=stress_tol, state_tol=stress_tol, tangent_tol=np.inf,
        )
        both = [s for s in report.steps if s.a_converged and s.b_converged]
        wrong = [s for s in both if not s.passed]
        n_autodiff_unsolved = sum(not s.b_converged for s in report.steps)
        worst = _worst(
            [s.stress_rel_err for s in both]
            + [e for s in both for e in s.state_rel_err.values()] or [0.0]
        )
        text = f"{worst:.1e}"
        if n_autodiff_unsolved:
            text += f" (autodiff did not solve {n_autodiff_unsolved})"
        record("hand-written vs autodiff", not wrong, text, report)

    tangent_integrator: "StressIntegrator" = reference
    declared = None
    if hand_rm and cls.user_defined_tangent is not MaterialModel.user_defined_tangent:
        tangent_integrator = PythonAnalyticalIntegrator(model, raise_on_nonconverged=False)
        declared = hand_tangent_record(cls)
    status, detail, tangent_reports = _check_tangents(
        tangent_integrator, {name: results}, n_tangent_points, tangent_tol, model, declared,
    )
    record("tangent", status == "passed", detail.split(" (tol")[0], tangent_reports)

    if getattr(cls, "_hand_blocks", None):
        report = check_hand_blocks(model, plastic, tol=block_tol)
        exact = [b.worst_err for b in report.blocks.values() if b.tol == block_tol]
        record("hand blocks", report.passed, f"{(_worst(exact) if exact else 0.0):.1e}", report)

    umat = umat_of(model)
    if umat is not None and not umat.is_built():
        parts.append("UMAT not built")
    elif umat is not None:
        fortran = umat.integrator(model, raise_on_nonconverged=False)
        report = crosscheck(
            reference, fortran, solved, lockstep=True,
            stress_tol=stress_tol, state_tol=stress_tol, tangent_tol=10.0 * stress_tol,
        )
        errors = [report.max_stress_rel_err, *report.max_state_rel_err.values()]
        if report.max_tangent_rel_err is not None:
            errors.append(report.max_tangent_rel_err)
        text = f"{_worst(errors):.1e}"
        ok = report.passed
        if "parse_umat_failure" in umat.integrator_options:
            # The increments the Python model gave up on: the UMAT reads its
            # own convergence flag, and has to give up on them as well.
            n_solved_by_umat = 0
            for strain_inc, stress_n, state_n in unsolved:
                try:
                    n_solved_by_umat += bool(
                        fortran.stress_update(strain_inc, stress_n, state_n).converged
                    )
                except _SOLVER_ERRORS:
                    pass
            if n_solved_by_umat:
                ok = False
                text += (f" — converged on {n_solved_by_umat} increment(s) "
                         "the Python model did not solve")
        record("UMAT", ok, text, report)

    detail = f"{coverage}; " + ", ".join(parts)
    if failures:
        return RangeCase(
            name, "loading", "failed", f"{', '.join(failures)} failed — {detail}", reports,
        )
    return RangeCase(name, "loading", "passed", detail, reports)


def _check_variant(name: str, variant: MaterialModel, **tolerances) -> RangeCase:
    """One parameter variant through ``check_model``."""
    try:
        report: ModelReport = check_model(variant, **tolerances)
    except ValueError as exc:
        # No history can be built for it (no yield strain to scale by).
        return RangeCase(name, "parameters", "failed", f"could not be checked — {exc}")
    failed = [c for c in report.checks if c.status == "failed"]
    if failed:
        detail = "; ".join(f"{c.name}: {c.detail}" for c in failed)
        return RangeCase(name, "parameters", "failed", detail, {"model": report})
    ran = [c.name for c in report.checks if c.status == "passed"]
    return RangeCase(
        name, "parameters", "passed", f"{len(ran)} checks pass ({', '.join(ran)})",
        {"model": report},
    )


def check_range(
    model: MaterialModel,
    *,
    variants: Mapping[str, MaterialModel] | None = None,
    histories: Mapping[str, FieldHistory] | None = None,
    seed: int = 0,
    tangent_tol: float = 1e-5,
    block_tol: float = 1e-8,
    solution_tol: float = 1e-8,
    stress_tol: float = 1e-6,
    n_tangent_points: int = 4,
) -> RangeReport:
    """Check ``model`` at the ends of its parameter ranges and under random loading.

    Parameters
    ----------
    model : MaterialModel
        The model, with the parameters the variants and the random histories
        start from.
    variants : mapping of str to MaterialModel, or None
        Models to run through :func:`~returnmap.verification.check_model`.
        ``None`` (default): the variants the model's family declares (see
        the module docstring), or :func:`parameter_variants` when it
        declares none.  An empty mapping leaves the parameter cases out.
    histories : mapping of str to FieldHistory, or None
        Strain histories for the loading cases.  ``None`` (default):
        :func:`random_histories` with ``seed``.  An empty mapping leaves the
        loading cases out.
    seed : int, default 0
        Seed of the default random histories.
    tangent_tol, block_tol, solution_tol, stress_tol, n_tangent_points
        As in :func:`~returnmap.verification.check_model`, for both kinds of
        case.

    Returns
    -------
    RangeReport
        ``report.passed``, ``report.summary()``, and ``report[name]`` for one
        case with the reports behind it.

    Notes
    -----
    A parameter variant is judged exactly as ``check_model`` judges it, on
    the standard histories and those the family declares: there, the autodiff
    return mapping has to converge.  The loading cases apply the rule of the
    module docstring instead — a wrong answer fails, an increment that was
    not solved is counted.

    The UMAT is compared at every variant and on every random history when
    it is declared and built; where it is not built the cases say so and
    pass on the Python side alone.

    The bounds of ``param_meta`` are the *plausible* range a model documents.
    A model need not work at their ends — but then the range it documents is
    wider than the one it was verified on, and this report shows where.
    """
    if not isinstance(model, MaterialModel):
        raise TypeError(
            f"check_range takes a model instance, got {model!r}; "
            "for a class, check cls.at_defaults()"
        )
    tolerances = dict(
        tangent_tol=tangent_tol, block_tol=block_tol, solution_tol=solution_tol,
        stress_tol=stress_tol, n_tangent_points=n_tangent_points,
    )
    notes: list[str] = []
    if variants is None:
        variants = _declared_variants(model)
    if variants is None:
        variants = parameter_variants(model)
        if not variants:
            notes.append(
                f"{type(model).__name__} documents no parameter range "
                "(param_meta bounds): no parameter case"
            )
    if histories is None:
        histories = random_histories(model, seed=seed)
    cases = [_check_variant(name, variant, **tolerances) for name, variant in variants.items()]
    cases += [
        _check_loading(model, name, load, **tolerances) for name, load in histories.items()
    ]
    return RangeReport(type(model).__name__, cases, notes)
