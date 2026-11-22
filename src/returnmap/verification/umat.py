"""Declaring the Fortran UMAT that belongs to a model.

A model that has been ported to a UMAT says so in its family's
``verification.py``, next to the models and loading histories it is verified
on::

    from returnmap.verification import Umat

    UMATS = [
        Umat(YUKinematic3D, module="yu_kinematic_3d",
             sources=["fortran/yu_kinematic_3d.f90"], libraries=["lapack", "blas"]),
    ]

That one declaration is what the tools read:

* ``returnmap build --family <package>`` compiles every UMAT a package of
  models declares — sources, the ABAQUS stubs, the libraries to link;
* :func:`~returnmap.verification.check_model` runs the UMAT through the same
  loading histories as the Python model and compares them step by step, when
  the module has been built — and then calls the subroutine ABAQUS itself
  calls (``umat``), which unpacks ``PROPS`` and ``STATEV``, and compares it
  with the one it wraps.
"""

from __future__ import annotations

import importlib
import importlib.util
import inspect
import pkgutil
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from returnmap.core.material import MaterialModel
    from returnmap.simulation.integrator import FortranIntegrator


@dataclass(frozen=True)
class Umat:
    """The Fortran UMAT of one model class.

    Parameters
    ----------
    model : type[MaterialModel]
        The model class the UMAT is the port of — one dimension of a family
        (``YUKinematic3D``), since a UMAT is written for one.
    module : str
        Name of the f2py module the sources are compiled into.
    sources : sequence of str
        Fortran source files, relative to the directory of the model's
        module (``"fortran/yu_kinematic_3d.f90"`` for a model defined in
        ``.../yu_kinematic/model.py``), so the declaration moves with the
        family directory.
    libraries : sequence of str
        Libraries to link (``["lapack", "blas"]``).
    subroutine : str or None
        The subroutine that does one stress update.  ``None`` (default): the
        one named like the module.
    props : sequence of str or None
        The model's parameter names in the order the UMAT takes them — its
        ``PROPS``.  ``None`` (default): ``param_names``, in that order.  Give
        it when the UMAT orders them differently or takes only some of them
        (a parameter that selects a branch on the Python side, say).
    integrator_options : mapping
        Extra keyword arguments for
        :meth:`FortranIntegrator.from_model
        <returnmap.simulation.FortranIntegrator.from_model>` — typically
        ``parse_umat_failure``, to read the UMAT's own convergence flag.
    entry : str or None
        The subroutine ABAQUS calls, in the same module.  Default ``"umat"``.
        ``None``: the module has none, or it is not to be checked.
    statev : sequence of str or None
        The model's state variables in the order ``STATEV`` holds them, a
        tensor taking as many slots as it has components.  ``None``
        (default): the declaration order, without the stress.  Give it when
        the entry lays them out differently.

    Notes
    -----
    Two subroutines, two conventions.

    The **core** (``subroutine``) does one stress update with plain
    arguments, and is what the Python model is compared with.  Its inputs
    are the parameters (``props``), the stress, the state variables in
    declaration order, the strain increment; its outputs are the stress, the
    state variables in the same order, the tangent, and then whatever else
    it reports (a convergence flag).

    The **entry** (``entry``) has the argument list ABAQUS prescribes and
    only unpacks ``PROPS`` and ``STATEV`` for the core.  It is declared as
    in the shipped UMATs — ``STRESS``, ``STATEV`` and ``PNEWDT`` are
    ``intent(inout)``; ``DDSDDE`` and the energies ``intent(out)``; the rest
    ``intent(in)`` — which is what lets it be called without a description
    of its own.  It reports a step it could not solve by lowering
    ``PNEWDT``, and is then expected to leave ``STRESS`` and ``STATEV`` as
    they were.  A ``STATEV`` with slots that are not state variables of the
    model cannot be described by ``statev``; declare ``entry=None`` for
    such a module and test its entry by hand.
    """

    model: "type[MaterialModel]"
    module: str
    sources: Sequence[str]
    libraries: Sequence[str] = ()
    subroutine: str | None = None
    props: Sequence[str] | None = None
    integrator_options: Mapping[str, Any] = field(default_factory=dict)
    entry: str | None = "umat"
    statev: Sequence[str] | None = None

    def __post_init__(self) -> None:
        if self.props is not None:
            unknown = [name for name in self.props if name not in self.model.param_names]
            if unknown:
                raise ValueError(
                    f"Umat({self.model.__name__}, module={self.module!r}): props names "
                    f"{unknown}, which are not parameters of the model "
                    f"{list(self.model.param_names)}"
                )
        if self.statev is not None:
            states = [name for name in self.model.state_names if name != "stress"]
            if sorted(self.statev) != sorted(states):
                raise ValueError(
                    f"Umat({self.model.__name__}, module={self.module!r}): statev names "
                    f"{list(self.statev)}; it has to hold each state variable of the "
                    f"model once: {states}"
                )

    def source_paths(self) -> list[Path]:
        """The sources as absolute paths."""
        source_file = inspect.getsourcefile(self.model)
        if source_file is None:
            raise ValueError(f"cannot locate the module {self.model.__name__} is defined in")
        base = Path(source_file).parent
        return [base / source for source in self.sources]

    def is_built(self) -> bool:
        """Whether the compiled module can be imported."""
        _build_dir_on_path()
        try:
            return importlib.util.find_spec(self.module) is not None
        except (ImportError, ValueError):
            return False

    def integrator(self, model: "MaterialModel", **options) -> "FortranIntegrator":
        """The UMAT as an integrator, with ``model``'s parameters.

        Raises
        ------
        ModuleNotFoundError
            The module has not been built.
        """
        from returnmap.simulation.integrator import FortranIntegrator, FortranModule

        _build_dir_on_path()
        declared: dict[str, Any] = dict(self.integrator_options)
        if self.props is not None:
            names = tuple(self.props)
            declared["param_fn"] = lambda: tuple(getattr(model, name) for name in names)
        return FortranIntegrator.from_model(
            FortranModule(self.module), self.subroutine or self.module, model,
            **{**declared, **options},
        )

    def has_entry(self) -> bool:
        """Whether the built module has the ABAQUS entry this declaration names."""
        if self.entry is None or not self.is_built():
            return False
        return hasattr(importlib.import_module(self.module), self.entry)

    def entry_integrator(self, model: "MaterialModel", **options) -> "FortranIntegrator":
        """The ABAQUS entry as an integrator, with ``model``'s parameters.

        The step is passed the way ABAQUS passes it — the parameters as
        ``PROPS``, the state variables as ``STATEV`` in the declared layout —
        and read back the same way.  A step on which the entry lowers
        ``PNEWDT`` is one it did not converge on.

        Raises
        ------
        ModuleNotFoundError
            The module has not been built.
        AttributeError
            The module has no subroutine of the entry's name.
        """
        from returnmap.simulation.integrator import FortranIntegrator

        if self.entry is None:
            raise AttributeError(f"Umat({self.module!r}) declares no ABAQUS entry")
        _build_dir_on_path()
        module = importlib.import_module(self.module)
        names = tuple(self.props) if self.props is not None else tuple(model.param_names)
        return FortranIntegrator.from_model(
            _AbaqusEntry(getattr(module, self.entry), model, len(names), self.statev),
            self.entry, model,
            param_fn=lambda: tuple(getattr(model, name) for name in names),
            parse_umat_failure=_pnewdt_lowered,
            **options,
        )


def _pnewdt_lowered(ret: tuple) -> str | None:
    """The last value an ``_AbaqusEntry`` call returns is ``PNEWDT``."""
    pnewdt = float(ret[-1])
    return None if pnewdt >= 1.0 else f"the UMAT asked for a smaller increment (PNEWDT = {pnewdt:g})"


class _AbaqusEntry:
    """The ABAQUS entry of a compiled UMAT, called with the arguments of the core.

    :class:`~returnmap.simulation.FortranIntegrator` calls
    ``call(name, *parameters, stress, *state, strain_inc)`` and reads
    ``(stress, *state, tangent, ...)`` back.  This packs that call into the
    argument list ABAQUS prescribes, calls the entry, and unpacks what it
    left in ``STRESS`` and ``STATEV`` into the same shape, with ``PNEWDT`` at
    the end.

    What the entry is not given: the total strain (zero), a time, a
    temperature, a deformation gradient (identity).  A UMAT that reads those
    is outside what this compares.
    """

    def __init__(self, entry, model: "MaterialModel", n_props: int,
                 statev: Sequence[str] | None) -> None:
        self._entry = entry
        self._n_props = n_props
        self._declared = [name for name in model.state_names if name != "stress"]
        self._layout = list(statev) if statev is not None else self._declared
        initial = model.initial_state()
        self._shapes = {name: np.shape(initial[name]) for name in self._declared}
        self._ndi, self._nshr = model.dimension.ndi, model.dimension.nshr
        self._cmname = type(model).__name__[:80].ljust(80).encode()

    def call(self, name: str, *args):
        props = np.array(args[:self._n_props], dtype=np.float64)
        stress = np.array(args[self._n_props], dtype=np.float64)        # a copy: updated in place
        values = dict(zip(self._declared, args[self._n_props + 1:-1]))
        strain_inc = np.asarray(args[-1], dtype=np.float64)
        statev = np.concatenate(
            [np.ravel(np.asarray(values[key], dtype=np.float64)) for key in self._layout]
            or [np.zeros(0)]
        )
        pnewdt = np.array(1.0)
        zeros, identity = np.zeros(stress.size), np.eye(3)
        ret = self._entry(
            stress, statev, zeros, strain_inc, np.zeros(2), 1.0, 0.0, 0.0,
            np.zeros(1), np.zeros(1), self._cmname, self._ndi, self._nshr,
            props, np.zeros(3), identity, pnewdt, 1.0, identity, identity,
            1, 1, 1, 1, 1, 1,
        )
        tangent = np.asarray(ret[0] if isinstance(ret, tuple) else ret, dtype=np.float64)
        unpacked, start = {}, 0
        for key in self._layout:
            size = int(np.prod(self._shapes[key], dtype=int))
            chunk = statev[start:start + size]
            unpacked[key] = chunk.reshape(self._shapes[key]) if self._shapes[key] else chunk[0]
            start += size
        return (stress, *(unpacked[key] for key in self._declared), tangent, float(pnewdt))


def _build_dir_on_path() -> None:
    """Make modules built by ``returnmap build`` importable (<project>/fortran)."""
    from returnmap.cli import build_dir

    path = str(build_dir())
    if path not in sys.path:
        sys.path.insert(0, path)


def declared_umats(package: str) -> list[Umat]:
    """The UMATs a package of models declares, its subpackages included.

    ``package`` is a dotted name: a family (``"returnmap.models.yu_kinematic"``),
    or a package of families (``"returnmap.models"``).  A UMAT is declared in
    the ``UMATS`` list of a module named ``verification``.
    """
    root = importlib.import_module(package)
    names = [package]
    if hasattr(root, "__path__"):
        names += [
            info.name
            for info in pkgutil.walk_packages(root.__path__, prefix=f"{package}.")
            if info.ispkg and not info.name.endswith(".tests")
        ]
    found: list[Umat] = []
    for name in names:
        if importlib.util.find_spec(f"{name}.verification") is None:
            continue
        module = importlib.import_module(f"{name}.verification")
        for umat in getattr(module, "UMATS", ()):
            if not isinstance(umat, Umat):
                raise TypeError(
                    f"{name}.verification.UMATS holds {umat!r}; "
                    "declare each UMAT with returnmap.verification.Umat"
                )
            found.append(umat)
    return found


def umat_of(model: "MaterialModel") -> Umat | None:
    """The UMAT declared for exactly this model's class, or ``None``."""
    package = type(model).__module__.rpartition(".")[0]
    if not package or importlib.util.find_spec(f"{package}.verification") is None:
        return None
    module = importlib.import_module(f"{package}.verification")
    for umat in getattr(module, "UMATS", ()):
        if umat.model is type(model):
            return umat
    return None
