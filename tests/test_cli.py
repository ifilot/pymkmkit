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
    assert "Forward barrier: 1.030996 eV (inc. ZPE-corr: 0.030996)" in result.output
    assert "Reverse barrier: 1.049594 eV (inc. ZPE-corr: 0.049594)" in result.output


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
    assert "Forward barrier: 0.524797 eV (inc. ZPE-corr: 0.024797)" in result.output


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
    assert "Forward barrier: 1.061992 eV (inc. ZPE-corr: 0.061992)" in result.output


def test_read_network_cli_adsorption_heat_with_zpe(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    gas_file = states_dir / "gas.yaml"
    empty_file = states_dir / "empty.yaml"
    ads_file = states_dir / "ads.yaml"

    gas_file.write_text(
        """
energy:
  electronic: -0.2
vibrations:
  frequencies_cm-1: [100.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    empty_file.write_text(
        """
energy:
  electronic: -1.0
vibrations:
  frequencies_cm-1: [200.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    ads_file.write_text(
        """
energy:
  electronic: -1.5
vibrations:
  frequencies_cm-1: [600.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: GAS
    file: states/gas.yaml
  - name: EMPTY
    file: states/empty.yaml
  - name: ADS
    file: states/ads.yaml
network:
  - name: adsorption
    type: ads
    reaction: GAS + * => ADS*
    is:
      - name: GAS
        stoichiometry: 1
      - name: EMPTY
        stoichiometry: 1
    fs:
      - name: ADS
        stoichiometry: 1
    normalization: 1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["read_network", str(network_file)])

    assert result.exit_code == 0
    assert "Reaction: GAS + * => ADS*" in result.output
    assert "Adsorption heat: -0.281402 eV (inc. ZPE-corr: 0.018598)" in result.output


def test_evaluate_paths_cli_sums_reaction_heats(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    is_file = states_dir / "is.yaml"
    ts_file = states_dir / "ts.yaml"
    fs_file = states_dir / "fs.yaml"
    gas_file = states_dir / "gas.yaml"

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
    gas_file.write_text(
        """
energy:
  electronic: -0.2
vibrations:
  frequencies_cm-1: [100.0]
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
  - name: GAS
    file: states/gas.yaml
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: surface_step
    reaction: IS* => FS*
    forward:
      ts:
        - name: TS
      is:
        - name: IS
    backward:
      ts:
        - name: TS
      is:
        - name: FS
  - name: adsorption_step
    type: ads
    reaction: GAS + * => IS*
    is:
      - name: GAS
      - name: FS
    fs:
      - name: IS
paths:
  - name: combined
    steps:
      - name: surface_step
      - name: adsorption_step
        factor: 2
""".strip()
        + "\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["evaluate_paths", str(network_file)])

    assert result.exit_code == 0
    assert "Path combined: 1.406199 eV" in result.output


def test_build_ped_cli_writes_output_file(tmp_path):
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
  - name: surface_step
    reaction: IS* => FS*
    forward:
      ts:
        - name: TS
      is:
        - name: IS
    backward:
      ts:
        - name: TS
      is:
        - name: FS
paths:
  - name: methanation
    steps:
      - name: surface_step
        factor: 2.5
        label: state_1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output_file = tmp_path / "ped.png"

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["build_ped", str(network_file), "methanation", str(output_file)],
    )

    assert result.exit_code == 0
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    assert f"Potential energy diagram written to: {output_file}" in result.output


def test_build_ped_cli_uses_placeholder_labels_when_missing(tmp_path):
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
  - name: surface_step
    reaction: IS* => FS*
    forward:
      ts:
        - name: TS
      is:
        - name: IS
    backward:
      ts:
        - name: TS
      is:
        - name: FS
paths:
  - name: methanation
    steps:
      - name: surface_step
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output_file = tmp_path / "ped_placeholder.png"

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["build_ped", str(network_file), "methanation", str(output_file)],
    )

    assert result.exit_code == 0
    assert output_file.exists()
    assert output_file.stat().st_size > 0
