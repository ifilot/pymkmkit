import zipfile
from pathlib import Path

import yaml
from click.testing import CliRunner

from pymkmkit.cli import cli


DATA_DIR = Path(__file__).parent / "data"


def _extract_outcar(zip_name, tmp_path):
    zip_path = DATA_DIR / zip_name
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(tmp_path)
    return tmp_path / zip_name.replace(".zip", "")


def test_freq2yaml_cli_writes_frequency_yaml(tmp_path):
    outcar = _extract_outcar("OUTCAR_Ru1121_C.zip", tmp_path)
    output = tmp_path / "freq.yaml"

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["freq2yaml", str(outcar), "--average-pairs", "-o", str(output)],
    )

    assert result.exit_code == 0
    assert output.exists()

    parsed = yaml.safe_load(output.read_text())
    assert parsed["calculation"]["type"] == "frequency"
    assert parsed["vibrations"]["paired_modes_averaged"] is True


def test_opt2yaml_cli_writes_optimization_yaml(tmp_path):
    outcar = _extract_outcar("OUTCAR_Ni311_C.zip", tmp_path)
    output = tmp_path / "opt.yaml"

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["opt2yaml", str(outcar), "-o", str(output)],
    )

    assert result.exit_code == 0
    assert output.exists()

    parsed = yaml.safe_load(output.read_text())
    assert parsed["calculation"]["type"] == "optimization"
    assert "vibrations" not in parsed
