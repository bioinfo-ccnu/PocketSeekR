"""Installed package resources and public command entry points."""

from dataclasses import asdict
from importlib import metadata, resources
import subprocess
import sys

import yaml

from pocketseekr import PocketSettings, __version__


def test_installed_defaults_are_available_as_package_data():
    resource = resources.files("pocketseekr").joinpath("data/default.yaml")
    assert yaml.safe_load(resource.read_text()) == asdict(PocketSettings())


def test_installed_version_matches_distribution_and_cli(tmp_path):
    assert __version__ == metadata.version("pocketseekr")
    result = subprocess.run(
        [sys.executable, "-m", "pocketseekr", "--version"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip() == f"PocketSeekR {__version__}"
