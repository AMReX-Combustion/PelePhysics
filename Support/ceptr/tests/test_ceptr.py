"""Tests for cepter."""

import pathlib

import cantera as ct

import ceptr.converter as converter
from ceptr import __version__


def mechanism_path(mname):
    """Determine mechanism path."""
    this_file_dir = pathlib.Path(__file__).parent.resolve()
    model_path = "Mechanism"
    return this_file_dir.parents[1] / model_path / mname


def test_version():
    """Test the version."""
    assert __version__ == "0.1.0"


def test_air():
    """Test mechanism generation of air."""
    mech_path = mechanism_path("air")
    fname = mech_path / "mechanism.yaml"
    mechanism = ct.Solution(fname)
    conv = converter.Converter(mechanism)
    conv.writer()
    conv.formatter()


def test_lidryer():
    """Test mechanism generation of LiDryer."""
    mech_path = mechanism_path("LiDryer")
    fname = mech_path / "mechanism.yaml"
    mechanism = ct.Solution(fname)
    conv = converter.Converter(mechanism)
    conv.writer()
    conv.formatter()


def test_dodecane_lu():
    """Test mechanism generation of dodecane_lu."""
    mech_path = mechanism_path("dodecane_lu")
    fname = mech_path / "mechanism.yaml"
    mechanism = ct.Solution(fname)
    conv = converter.Converter(mechanism)
    conv.writer()
    conv.formatter()
