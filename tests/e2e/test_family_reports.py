"""The report a family's README records is what ``returnmap check`` prints.

A family written in the specification format (``docs/guide/spec-driven.md``)
keeps the report of its checks in its ``README.md``, under 検証 / 報告: the
output of ``returnmap check --family <package> --require-built --range`` in an
environment where its UMAT is built.  A change to the model, to its UMAT or to
what the checks print leaves that record behind; this compares the two, numbers
aside.

Needs the UMATs built (``returnmap build --family returnmap.models``) and runs
the range checks of every model of the family, so it is both ``fortran`` and
``slow``.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).parents[2]
MODELS = REPO / "src" / "returnmap" / "models"
_NUMBER = re.compile(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?")
_REPORT = re.compile(r"^### 報告\n.*?^```none\n(.*?)^```$", flags=re.M | re.S)


def _families_with_a_recorded_report() -> list[str]:
    return sorted(
        readme.parent.name for readme in MODELS.glob("*/README.md")
        if _REPORT.search(readme.read_text(encoding="utf-8"))
    )


def _masked(line: str) -> str:
    """A line with its numbers and its column padding taken out."""
    return " ".join(_NUMBER.sub("#", line).split())


def test_the_families_in_the_specification_format_record_a_report():
    assert _families_with_a_recorded_report() == ["ow_kinematic"]


@pytest.mark.fortran
@pytest.mark.slow
@pytest.mark.parametrize("family", _families_with_a_recorded_report())
def test_the_recorded_report_is_what_the_check_prints(family):
    from returnmap.verification.umat import declared_umats

    package = f"returnmap.models.{family}"
    unbuilt = [umat.module for umat in declared_umats(package) if not umat.is_built()]
    if unbuilt:
        pytest.skip(f"not built: {', '.join(unbuilt)} — returnmap build --family {package}")
    env = dict(os.environ, PYTHONPATH=str(REPO / "src"), PYTHONDONTWRITEBYTECODE="1",
               PYTHONIOENCODING="utf-8")
    check = subprocess.run(
        [sys.executable, "-m", "returnmap", "check", "--family", package,
         "--require-built", "--range"],
        capture_output=True, text=True, encoding="utf-8", cwd=REPO, env=env,
    )

    readme = (MODELS / family / "README.md").read_text(encoding="utf-8")
    recorded = [_masked(line) for line in _REPORT.search(readme).group(1).splitlines()
                if line.strip()]
    printed = [_masked(line) for line in check.stdout.splitlines() if line.strip()]
    assert recorded == printed
