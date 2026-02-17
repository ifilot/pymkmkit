import zipfile
from pathlib import Path

import yaml
import pytest
from ase.io import read

from pymkmkit.vasp_freq import (
    average_mode_pairs,
    frequencies_from_partial_hessian,
    parse_vasp_frequency,
    parse_vasp_optimization,
)
from pymkmkit.yaml_writer import write_yaml


DATA_DIR = Path(__file__).parent / "data"


def _extract_outcar(zip_name, tmp_path):
    zip_path = DATA_DIR / zip_name
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(tmp_path)
    return tmp_path / zip_name.replace(".zip", "")


def _outcar_sigma0_energies(outcar_path):
    energies = []

    for line in Path(outcar_path).read_text(errors="ignore").splitlines():
        if "energy  without entropy=" in line and "energy(sigma->0)" in line:
            energies.append(float(line.split("=")[-1].strip()))

    return energies


def test_parse_outcar_zip(tmp_path):
    outcar = _extract_outcar("OUTCAR_Ni311_C.zip", tmp_path)

    data = parse_vasp_frequency(outcar)

    # basic structural checks
    assert data["structure"]["formula"].startswith("Ni")
    assert data["structure"]["n_atoms"] > 0

    # ensure lattice vectors parsed
    assert len(data["structure"]["lattice_vectors"]) == 3

    # ensure coordinates exist
    assert len(data["structure"]["coordinates_direct"]) > 0

    # energy should match the first ionic step (not the last)
    sigma0_energies = _outcar_sigma0_energies(outcar)
    assert isinstance(data["energy"]["electronic"], float)
    assert data["energy"]["electronic"] == sigma0_energies[0]
    assert data["energy"]["electronic"] != sigma0_energies[-1]

    # vibrations parsed
    assert len(data["vibrations"]["frequencies_cm-1"]) > 0
    assert "partial_hessian" in data["vibrations"]


def test_parse_ru1121_c_with_pair_averaging(tmp_path):
    outcar = _extract_outcar("OUTCAR_Ru1121_C.zip", tmp_path)

    data = parse_vasp_frequency(outcar, average_pairs=True)

    assert data["structure"]["formula"].startswith("Ru")
    assert data["vibrations"]["paired_modes_averaged"] is True
    assert data["vibrations"]["imaginary_cm-1"] is None
    assert data["vibrations"].get("pairing_note") is None


def test_parse_ru1121_c_ch_ts_with_pair_averaging_and_ts_note(tmp_path):
    outcar = _extract_outcar("OUTCAR_Ru1121_C_CH_TS.zip", tmp_path)

    data = parse_vasp_frequency(outcar, average_pairs=True)

    vibrations = data["vibrations"]
    assert vibrations["paired_modes_averaged"] is True
    assert vibrations["imaginary_cm-1"] is not None
    assert len(vibrations["imaginary_cm-1"]) == 1
    assert vibrations["imaginary_cm-1"][0] < 0
    assert "Odd number of real modes" in vibrations["pairing_note"]
    assert "Dropped unpaired real mode" in vibrations["pairing_note"]


def test_parse_outcar_as_optimization_uses_last_ionic_step(tmp_path):
    outcar = _extract_outcar("OUTCAR_Ni311_C.zip", tmp_path)

    data = parse_vasp_optimization(outcar)

    last_atoms = read(outcar, index=-1)

    assert data["calculation"]["type"] == "optimization"
    assert "vibrations" not in data

    expected_last_coords = [
        f"{atom.symbol} " + " ".join(f"{x:.8f}" for x in pos)
        for atom, pos in zip(last_atoms, last_atoms.get_scaled_positions())
    ]
    assert data["structure"]["coordinates_direct"] == expected_last_coords
    assert isinstance(data["energy"]["electronic"], float)


def test_average_mode_pairs_also_pairs_imaginary_modes():
    real_freqs = [100.0, 102.0, 150.0, 154.0]
    imag_freqs = [-400.0, -396.0]

    averaged_real, averaged_imag, note = average_mode_pairs(real_freqs, imag_freqs)

    assert averaged_real == [101.0, 152.0]
    assert averaged_imag == [-398.0]
    assert note is None


def test_parse_optimization_falls_back_when_ase_cannot_parse_positions(tmp_path):
    outcar = _extract_outcar("OUTCAR_Ru1121_empty.zip", tmp_path)

    data = parse_vasp_optimization(outcar)

    assert data["calculation"]["type"] == "optimization"
    assert data["structure"]["formula"].startswith("Ru")
    assert data["structure"]["n_atoms"] == len(data["structure"]["coordinates_direct"])

    sigma0_energies = _outcar_sigma0_energies(outcar)
    assert data["energy"]["electronic"] == sigma0_energies[-1]


def test_partial_hessian_roundtrip_recovers_vasp_frequencies(tmp_path):
    outcar = _extract_outcar("OUTCAR_Ni311_C.zip", tmp_path)

    data = parse_vasp_frequency(outcar)
    yaml_path = tmp_path / "freq.yaml"
    write_yaml(data, yaml_path)

    loaded = yaml.safe_load(yaml_path.read_text())
    partial_hessian = loaded["vibrations"]["partial_hessian"]

    recovered = frequencies_from_partial_hessian(
        partial_hessian["dof_labels"],
        partial_hessian["matrix"],
        read(outcar),
    )

    reported = sorted(loaded["vibrations"]["frequencies_cm-1"], reverse=True)
    assert len(recovered) == len(reported)

    for rec, rep in zip(recovered, reported):
        assert rec == pytest.approx(rep, abs=2.0)
