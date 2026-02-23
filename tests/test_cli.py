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


def test_read_network_cli_prints_reaction_and_barriers(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    is_file = states_dir / "is.yaml"
    ts_file = states_dir / "ts.yaml"
    fs_file = states_dir / "fs.yaml"

    is_file.write_text(
        """
energy:
  electronic: -1.0
vibrations:
  frequencies_cm-1: [500.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    ts_file.write_text(
        """
energy:
  electronic: 0.0
vibrations:
  frequencies_cm-1: [1000.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    fs_file.write_text(
        """
energy:
  electronic: -2.0
vibrations:
  frequencies_cm-1: [200.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: IS
    file: states/is.yaml
  - name: FS
    file: states/fs.yaml
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: test step
    reaction: IS* => FS*
    forward:
      ts:
        - name: TS
          stoichiometry: 1
      is:
        - name: IS
          stoichiometry: 1
      normalization: 1
    backward:
      ts:
        - name: TS
          stoichiometry: 1
      is:
        - name: FS
          stoichiometry: 1
      normalization: 2
""".strip()
        + "\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["read_network", str(network_file)])

    assert result.exit_code == 0
    assert "Reaction: IS* => FS*" in result.output
    assert "Forward electronic barrier: 1.000000" in result.output
    assert "Forward ZPE correction: 0.030996" in result.output
    assert "Forward total barrier: 1.030996" in result.output
    assert "Reverse electronic barrier: 1.000000" in result.output
    assert "Reverse ZPE correction: 0.049594" in result.output
    assert "Reverse total barrier: 1.049594" in result.output


def test_read_network_cli_normalizes_zpe_when_not_paired(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    is_file = states_dir / "is.yaml"
    ts_file = states_dir / "ts.yaml"

    is_file.write_text(
        """
energy:
  electronic: -1.0
vibrations:
  frequencies_cm-1: [200.0]
  paired_modes_averaged: false
""".strip()
        + "\n",
        encoding="utf-8",
    )
    ts_file.write_text(
        """
energy:
  electronic: 0.0
vibrations:
  frequencies_cm-1: [1000.0]
  paired_modes_averaged: false
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: IS
    file: states/is.yaml
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: test step
    reaction: IS* => TS*
    forward:
      ts:
        - name: TS
          stoichiometry: 1
      is:
        - name: IS
          stoichiometry: 1
      normalization: 2
    backward:
      ts:
        - name: TS
          stoichiometry: 1
      is:
        - name: IS
          stoichiometry: 1
      normalization: 2
""".strip()
        + "\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["read_network", str(network_file)])

    assert result.exit_code == 0
    assert "Forward ZPE correction: 0.024797" in result.output


def test_read_network_cli_allows_missing_vibrations_as_zero_zpe(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    is_file = states_dir / "is.yaml"
    ts_file = states_dir / "ts.yaml"

    is_file.write_text(
        """
energy:
  electronic: -1.0
""".strip()
        + "\n",
        encoding="utf-8",
    )
    ts_file.write_text(
        """
energy:
  electronic: 0.0
vibrations:
  frequencies_cm-1: [1000.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: IS
    file: states/is.yaml
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: test step
    reaction: IS* => TS*
    forward:
      ts:
        - name: TS
          stoichiometry: 1
      is:
        - name: IS
          stoichiometry: 1
      normalization: 1
    backward:
      ts:
        - name: TS
          stoichiometry: 1
      is:
        - name: IS
          stoichiometry: 1
      normalization: 1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["read_network", str(network_file)])

    assert result.exit_code == 0
    assert "Forward electronic barrier: 1.000000" in result.output
    assert "Forward ZPE correction: 0.061992" in result.output
    assert "Forward total barrier: 1.061992" in result.output
