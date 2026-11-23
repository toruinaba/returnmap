"""The commands' own logic: what build hands to f2py, what check runs and reports.

Running an actual build needs gfortran and is covered in ``tests/e2e``.
"""

import sys
from pathlib import Path

from returnmap.cli import abaqus_stubs_path, f2py_command, main


def test_the_stubs_are_shipped_inside_the_package():
    path = abaqus_stubs_path()
    assert path.name == "abaqus_stubs.f90"
    assert path.is_file()
    import returnmap
    assert Path(returnmap.__file__).parent in path.parents
    source = path.read_text().lower()
    assert all(name in source for name in ("subroutine sinv", "subroutine sprinc", "subroutine rotsig"))


def test_plain_build_is_sources_then_module_name():
    cmd = f2py_command([Path("/src/a.f90"), Path("/src/b.f90")], "my_umat")
    assert cmd == [sys.executable, "-m", "numpy.f2py", "-c", "/src/a.f90", "/src/b.f90",
                   "-m", "my_umat"]


def test_stubs_are_compiled_in_first():
    cmd = f2py_command([Path("/src/umat.f90")], "my_umat", abaqus_stubs=True)
    assert cmd[4:6] == [str(abaqus_stubs_path()), "/src/umat.f90"]


def test_libraries_are_linked_after_the_module_name():
    cmd = f2py_command([Path("/src/umat.f90")], "my_umat", link=["lapack", "blas"])
    assert cmd[-4:] == ["-m", "my_umat", "-llapack", "-lblas"]


def test_stubs_command_prints_the_path(capsys, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["returnmap", "stubs"])
    main()
    assert capsys.readouterr().out.strip() == str(abaqus_stubs_path())


def _parsed_build_args(monkeypatch, *argv):
    """What the build command receives for a command line."""
    import returnmap.cli as cli

    seen = {}
    monkeypatch.setattr(cli, "cmd_build", lambda args: seen.update(vars(args)))
    monkeypatch.setattr(sys, "argv", ["returnmap", "build", *argv])
    cli.main()
    return seen


def test_build_takes_the_stubs_flag_and_repeated_link(monkeypatch):
    args = _parsed_build_args(
        monkeypatch, "umat.f90", "--name", "my_umat",
        "--abaqus-stubs", "--link", "lapack", "-l", "blas",
    )
    assert args["abaqus_stubs"] is True
    assert args["link"] == ["lapack", "blas"]


def test_build_defaults_to_no_stubs_and_no_libraries(monkeypatch):
    """A source that is not a UMAT (the smoke-test module) needs neither."""
    args = _parsed_build_args(monkeypatch, "test_basic.f90", "--name", "basic")
    assert args["abaqus_stubs"] is False
    assert args["link"] == []


# ---------------------------------------------------------------------------
# build --family: every UMAT a package of models declares
# ---------------------------------------------------------------------------

def _run_build(monkeypatch, *argv):
    """Run ``returnmap build`` with the compiler replaced by a recorder."""
    import returnmap.cli as cli

    built = []
    for check in ("_check_gfortran", "_check_f2py", "_check_meson"):
        monkeypatch.setattr(cli, check, lambda: None)
    monkeypatch.setattr(
        cli, "_build_module",
        lambda sources, name, **options: built.append((name, sources, options)),
    )
    monkeypatch.setattr(sys, "argv", ["returnmap", "build", *argv])
    cli.main()
    return built


def test_family_builds_what_the_models_declare(monkeypatch, tmp_path):
    built = _run_build(monkeypatch, "--family", "returnmap.models.yu_kinematic",
                       "--output-dir", str(tmp_path))
    assert [name for name, _, _ in built] == ["yu_kinematic_3d", "yu_kinematic_ps"]
    for name, sources, options in built:
        assert [p.name for p in sources] == [f"{name}.f90"]
        assert all(p.is_file() for p in sources)
        # A UMAT always needs the stubs; the libraries come from the declaration.
        assert options["abaqus_stubs"] is True
        assert options["link"] == ["lapack", "blas"]
        assert options["output_dir"] == tmp_path


def test_family_can_be_a_package_of_families(monkeypatch, tmp_path):
    built = _run_build(monkeypatch, "--family", "returnmap.models", "-o", str(tmp_path))
    assert [name for name, _, _ in built] == [
        "j2_isotropic_3d", "ow_kinematic_3d", "yu_kinematic_3d", "yu_kinematic_ps"]


def test_family_without_a_umat_has_nothing_to_build(monkeypatch, capsys):
    """Not an error: the same command serves a family before it has a UMAT."""
    built = _run_build(monkeypatch, "--family", "returnmap.models.af_kinematic")
    assert built == []
    out = capsys.readouterr().out
    assert "declares no UMAT" in out and "nothing to build" in out


def test_family_and_source_files_are_alternatives(monkeypatch, capsys):
    import pytest

    for argv in (["--family", "returnmap.models", "x.f90"], ["x.f90"], []):
        with pytest.raises(SystemExit) as exit_info:
            _run_build(monkeypatch, *argv)
        assert exit_info.value.code == 2
    assert "--family <package>" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------

def _run_check(monkeypatch, capsys, *argv):
    """Exit code and output of ``returnmap check <argv>``."""
    import pytest

    monkeypatch.setattr(sys, "argv", ["returnmap", "check", *argv])
    code = 0
    try:
        main()
    except SystemExit as exit_info:
        code = exit_info.code
    return code, capsys.readouterr().out


def test_check_prints_the_report_of_every_declared_model(monkeypatch, capsys):
    code, out = _run_check(monkeypatch, capsys, "--family", "returnmap.models.j2_isotropic")
    assert code == 0
    for name in ("J2Isotropic3D", "J2IsotropicPS", "J2Isotropic1D"):
        assert f"{name}: every check that applies passes" in out
    assert out.rstrip().endswith("3 model(s) checked: every one passes")


def test_check_exits_with_1_when_a_model_fails(monkeypatch, capsys):
    import returnmap.verification.model as module
    from returnmap.models import J2Isotropic3D

    class TangentOff(J2Isotropic3D):
        def user_defined_tangent(self, new, dlambda, inc):
            return 1.01 * super().user_defined_tangent(new, dlambda, inc)

    monkeypatch.setattr(module, "declared_models", lambda package: [TangentOff.at_defaults()])
    code, out = _run_check(monkeypatch, capsys, "--family", "anything")
    assert code == 1
    assert "FAIL hand-written tangent" in out
    # The table behind the failed line, without being asked for.
    assert "    hand-written tangent — uniaxial[0]:" in out
    assert "tangent check" in out
    assert "and 5 more" in out
    assert out.rstrip().endswith("1 model(s) checked: 1 failed — TangentOff")


def test_check_names_a_model_by_the_parameters_that_are_not_typical(monkeypatch, capsys):
    import returnmap.verification.model as module
    from returnmap.models import J2Isotropic1D

    models = [J2Isotropic1D.at_defaults(), J2Isotropic1D.at_defaults(H=0.0, sigma_y0=300.0)]
    monkeypatch.setattr(module, "declared_models", lambda package: models)
    code, out = _run_check(monkeypatch, capsys, "--family", "anything")
    assert code == 0
    assert "J2Isotropic1D (sigma_y0=300, H=0):" in out


def test_check_without_models_is_an_error(monkeypatch, capsys):
    code, out = _run_check(monkeypatch, capsys, "--family", "returnmap.verification")
    assert code == 1
    assert "declares no model to check" in out


def test_a_umat_that_is_not_built_fails_only_when_required(monkeypatch, capsys):
    import pytest

    if __import__("importlib").util.find_spec("j2_isotropic_3d") is not None:
        pytest.skip("j2_isotropic_3d is built here")
    family = ("--family", "returnmap.models.j2_isotropic")
    code, out = _run_check(monkeypatch, capsys, *family)
    assert code == 0 and "required:" not in out
    code, out = _run_check(monkeypatch, capsys, *family, "--require-built")
    assert code == 1
    assert "  required: UMAT — module j2_isotropic_3d is not built" in out
    assert "  required: UMAT entry — module j2_isotropic_3d is not built" in out
    assert out.rstrip().endswith("3 model(s) checked: 1 failed — J2Isotropic3D")


def test_check_with_range_runs_check_range_as_well(monkeypatch, capsys):
    import returnmap.verification.model as module
    from returnmap.models import J2Isotropic1D

    monkeypatch.setattr(module, "declared_models", lambda package: [J2Isotropic1D.at_defaults()])
    code, out = _run_check(monkeypatch, capsys, "--family", "anything", "--range", "--seed", "3")
    assert code == 0
    assert "J2Isotropic1D over its range: every case passes" in out
    assert "random single-step" in out
