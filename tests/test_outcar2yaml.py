import zipfile
from pathlib import Path

from pymkmkit.vasp_freq import parse_vasp_frequency


DATA_DIR = Path(__file__).parent / "data"


def _extract_outcar(zip_name, tmp_path):
    zip_path = DATA_DIR / zip_name
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(tmp_path)
    return tmp_path / zip_name.replace(".zip", "")


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

    # energy should exist
    assert isinstance(data["energy"]["electronic"], float)

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
