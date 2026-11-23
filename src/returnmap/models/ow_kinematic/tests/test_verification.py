"""The checks every model owes — :func:`returnmap.verification.check_model`.

The same few lines in every family; what is checked follows from what the
model declares, and on which histories from ``verification.py``.
"""

import pytest

from returnmap.models.ow_kinematic.verification import MODELS, UMATS
from returnmap.verification import check_model


def _model_id(model) -> str:
    return f"{type(model).__name__}, gamma={model.gamma:g}"


# Every model of MODELS a UMAT is declared for, with that UMAT.
UMAT_MODELS = [
    pytest.param(umat, model, id=f"{umat.module}, gamma={model.gamma:g}")
    for umat in UMATS
    for model in MODELS
    if type(model) is umat.model
]


@pytest.mark.parametrize("model", MODELS, ids=_model_id)
def test_passes_every_check_that_applies(model):
    report = check_model(model)
    assert report.passed, report.summary()


@pytest.mark.fortran
@pytest.mark.parametrize("umat, model", UMAT_MODELS)
def test_the_umat_agrees_with_the_python_model(umat, model):
    """The compiled UMAT, step by step over the same histories, and its
    ABAQUS entry against the subroutine it wraps (needs the module built:
    ``returnmap build --family returnmap.models.ow_kinematic``)."""
    pytest.importorskip(umat.module, reason=f"{umat.module} is not built")
    report = check_model(model)
    assert report["UMAT"].status == "passed", report.summary()
    assert report["UMAT entry"].status == "passed", report.summary()
