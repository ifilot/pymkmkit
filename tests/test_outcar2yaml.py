import zipfile
from pathlib import Path

from pymkmkit.vasp_freq import parse_vasp_frequency


DATA_DIR = Path(__file__).parent / "data"
ZIP_PATH = DATA_DIR / "OUTCAR_Ni311_C.zip"


def test_parse_outcar_zip(tmp_path):
    # extract OUTCAR into temporary directory
    with zipfile.ZipFile(ZIP_PATH, "r") as z:
        z.extractall(tmp_path)

    outcar = tmp_path / "OUTCAR_Ni311_C"

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
