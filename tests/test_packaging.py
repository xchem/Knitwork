"""Tests that the wheel we publish to PyPI is installable and complete."""

import ast
import subprocess
import sys
import zipfile
from email.parser import Parser
from pathlib import Path

import pytest
from packaging.requirements import Requirement

PROJECT_ROOT = Path(__file__).parent.parent
PACKAGE_DIR = PROJECT_ROOT / "knitwork"


@pytest.fixture(scope="module")
def wheel(tmp_path_factory) -> zipfile.ZipFile:
    out_dir = tmp_path_factory.mktemp("dist")
    subprocess.run(
        [sys.executable, "-m", "build", "--wheel", "--outdir", str(out_dir)],
        cwd=PROJECT_ROOT,
        check=True,
    )
    (wheel_path,) = out_dir.glob("*.whl")
    with zipfile.ZipFile(wheel_path) as zf:
        yield zf


def _wheel_metadata(wheel: zipfile.ZipFile):
    (metadata_name,) = [n for n in wheel.namelist() if n.endswith(".dist-info/METADATA")]
    return Parser().parsestr(wheel.read(metadata_name).decode())


def _third_party_imports() -> set[str]:
    """Top-level module names imported by the package, excluding the stdlib."""
    names = set()
    for source in PACKAGE_DIR.glob("*.py"):
        tree = ast.parse(source.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0:
                names.add(node.module.split(".")[0])
    return {n for n in names if n not in sys.stdlib_module_names}


def test_wheel_contains_package_modules(wheel):
    names = wheel.namelist()
    for source in PACKAGE_DIR.glob("*.py"):
        assert f"knitwork/{source.name}" in names


def test_wheel_contains_fingerprint_feature_definitions(wheel):
    # Loaded at runtime by knitwork.tools.load_sig_factory()
    assert "knitwork/FeatureswAliphaticXenon.fdef" in wheel.namelist()


def test_declared_dependencies_match_imports(wheel):
    # Every package we import must be declared (or users get ImportErrors)
    # and every declared package should be imported (or users install bloat).
    # All our import names currently match their distribution names.
    requires = _wheel_metadata(wheel).get_all("Requires-Dist") or []
    declared = {Requirement(r).name.lower() for r in requires if "extra ==" not in r}
    imported = {name.lower() for name in _third_party_imports()}
    assert declared == imported
