import zipfile
from pathlib import Path

from ase.io import read

from pymkmkit.vasp_freq import parse_vasp_frequency, parse_vasp_optimization


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