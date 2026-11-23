# mymodels — a model family built on returnmap

The layout a repository of your own models takes, complete and runnable:
J2 plasticity with Voce hardening, from the Python model to a verified UMAT.
The guide walks through it step by step (in Japanese):
<https://toruinaba.github.io/returnmap/guide/adding-a-model.html>.

```
pyproject.toml                  returnmap pinned to a release tag
Makefile, Dockerfile            the same checks in a container with gfortran
conftest.py                     only for running from a returnmap checkout
report.py                       prints the report of the common checks
umat_crosscheck.py              compares the built UMAT step by step (needs the build below)
templates/
  family-README.md              the specification of a new family, to fill in
  brief.md                      what to hand to whoever implements it
src/mymodels/voce/
  README.md                     the specification: formulation, properties, UMAT, verification
  model.py                      J2Voce (the physics), J2Voce1D, J2Voce3D (+ hand-written solver)
  verification.py               what is checked: MODELS, UMATS
  fortran/voce_3d.f90           the UMAT: the port of J2Voce3D
  tests/test_verification.py    runs the common checks
  tests/test_properties.py      the properties the specification states, one test each
fortran/                        where `returnmap build` puts compiled modules (created on build)
```

Copy the directory out of the checkout and it is a project of its own:

```bash
uv sync --extra dev
uv run pytest                                  # 20 passed, 1 skipped (the UMAT is not built)
uv run returnmap check --family mymodels       # what was checked, row by row

uv sync --extra dev --extra fortran            # needs gfortran
uv run returnmap build --family mymodels       # compiles what verification.py declares
uv run pytest                                  # 21 passed
uv run returnmap check --family mymodels --require-built --range
uv run python umat_crosscheck.py               # the comparison with the UMAT, taken apart
```

Without gfortran on the machine, the same in a container:

```bash
make docker-image                              # once
make docker-check                              # build, tests, report; an unbuilt UMAT fails here
```

A new family starts from its specification: copy `templates/family-README.md` to
`src/mymodels/<family>/README.md` and fill it in.  The guide describes how,
and how to hand it over for implementation (in Japanese):
<https://toruinaba.github.io/returnmap/guide/spec-driven.html>.
