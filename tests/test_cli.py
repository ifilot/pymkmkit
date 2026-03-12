import zipfile
from pathlib import Path

import yaml
import pytest
from click.testing import CliRunner

from pymkmkit.cli import cli
from pymkmkit.network_reader import _format_chemical_subscripts, read_network


DATA_DIR = Path(__file__).parent / "data"


def _extract_outcar(zip_name, tmp_path):
    zip_path = DATA_DIR / zip_name
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(tmp_path)
    default = tmp_path / zip_name.replace(".zip", "")
    if default.exists():
        return default
    return tmp_path


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
    assert parsed["calculation"]["version"] == "6.5.1"
    assert parsed["calculation"]["executed_at"] == "2025-11-02T10:24:49Z"
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
    assert parsed["calculation"]["version"] == "5.3.5"
    assert parsed["calculation"]["executed_at"] == "2019-05-22T13:55:02Z"
    assert "vibrations" not in parsed


def test_asevib2yaml_cli_writes_frequency_yaml_from_vib_folders(tmp_path):
    root = _extract_outcar("CeO2_Pd4_CO.zip", tmp_path)
    outcar = root / "OUTCAR"
    output = tmp_path / "asevib.yaml"

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["asevib2yaml", str(outcar), "-o", str(output)],
    )

    assert result.exit_code == 0
    assert output.exists()

    parsed = yaml.safe_load(output.read_text())
    assert parsed["calculation"]["type"] == "frequency"
    assert len(parsed["vibrations"]["partial_hessian"]["dof_labels"]) == 12


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
    assert "Reaction: IS* => FS* (surf)" in result.output
    assert "Forward barrier:" in result.output
    assert "1.0310 eV" in result.output
    assert "0.0310 eV" in result.output
    assert "Reverse barrier:" in result.output
    assert "1.0496 eV" in result.output
    assert "0.0496 eV" in result.output


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
    assert "Forward barrier:" in result.output
    assert "0.5248 eV" in result.output
    assert "0.0248 eV" in result.output


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
    assert "Forward barrier:" in result.output
    assert "1.0620 eV" in result.output
    assert "0.0620 eV" in result.output


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
    assert "Reaction: GAS + * => ADS* (ads)" in result.output
    assert "Adsorption heat:" in result.output
    assert "-0.2814 eV" in result.output
    assert "0.0186 eV" in result.output



def test_read_network_cli_displays_kj_per_mol_option(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    is_file = states_dir / "is.yaml"
    ts_file = states_dir / "ts.yaml"

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
      is:
        - name: IS
    backward:
      ts:
        - name: TS
      is:
        - name: IS
""".strip()
        + "\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["read_network", str(network_file), "--unit", "kj/mol"])

    assert result.exit_code == 0
    assert "Reaction: IS* => TS* (surf)" in result.output
    assert "Forward barrier:" in result.output
    assert "99.5 kJ/mol" in result.output
    assert "3.0 kJ/mol" in result.output

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


def test_build_ped_cli_uses_path_startlabel_for_state_0(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    for state_path in ("is.yaml", "ts.yaml", "fs.yaml"):
        (states_dir / state_path).write_text(
            """
energy:
  electronic: 0.0
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
    startlabel: CO + 3H2
    steps:
      - name: surface_step
        label: CO* + 3H2
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output_file = tmp_path / "ped_startlabel.png"

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["build_ped", str(network_file), "methanation", str(output_file)],
    )

    assert result.exit_code == 0
    assert output_file.exists()


def test_format_chemical_subscripts_only_subscripts_attached_digits():
    assert _format_chemical_subscripts("CO* + 3H2") == "CO* + 3H$_{2}$"


def test_network2fnf_cli_writes_nodes_and_edges_ev(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    (states_dir / "a.yaml").write_text(
        """
energy:
  electronic: -1.0
vibrations:
  frequencies_cm-1: [100.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_dir / "b.yaml").write_text(
        """
energy:
  electronic: -0.5
vibrations:
  frequencies_cm-1: [200.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_dir / "ts.yaml").write_text(
        """
energy:
  electronic: 0.3
vibrations:
  frequencies_cm-1: [300.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A*
    file: states/a.yaml
  - name: B*
    file: states/b.yaml
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: surf step
    type: surf
    reaction: A* => B*
    forward:
      ts:
        - name: TS
      is:
        - name: A*
    backward:
      ts:
        - name: TS
      is:
        - name: B*
  - name: ads step
    type: ads
    reaction: A + * => A*
    is:
      - name: A*
    fs:
      - name: B*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    runner = CliRunner()
    result = runner.invoke(cli, ["network2fnf", str(network_file), "-o", str(output)])

    assert result.exit_code == 0
    payload = yaml.safe_load(output.read_text())

    assert payload["pymkmkit"]["units"] == "eV"
    assert payload["pymkmkit"]["energy_type"] == "elec+zpe"
    assert payload["nodes"] == [{"label": "A*"}, {"label": "B*"}]

    surf_edge = payload["edges"][0]
    assert surf_edge["type"] == "surf"
    assert set(surf_edge) >= {"forward", "backward", "nodes"}

    ads_edge = payload["edges"][1]
    assert ads_edge["type"] == "ads"
    assert "ads" in ads_edge
    assert ads_edge["nodes"] == ["A*", "B*"]


def test_network2fnf_cli_structures_copies_yaml_and_sets_node_structure_paths(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    (states_dir / "a.yaml").write_text(
        """
energy:
  electronic: -1.0
vibrations:
  frequencies_cm-1: [100.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_dir / "b.yaml").write_text(
        """
energy:
  electronic: -0.5
vibrations:
  frequencies_cm-1: [200.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_dir / "ts.yaml").write_text(
        """
energy:
  electronic: 0.3
vibrations:
  frequencies_cm-1: [300.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A*
    file: states/a.yaml
  - name: B*
    file: states/b.yaml
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: surf step
    type: surf
    reaction: A* => B*
    forward:
      ts:
        - name: TS
      is:
        - name: A*
    backward:
      ts:
        - name: TS
      is:
        - name: B*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["network2fnf", str(network_file), "-o", str(output), "--structures", "structures"],
    )

    assert result.exit_code == 0
    payload = yaml.safe_load(output.read_text())

    assert payload["nodes"] == [
        {"label": "A*", "structure": "structures/a.yaml"},
        {"label": "B*", "structure": "structures/b.yaml"},
    ]
    assert (tmp_path / "structures" / "a.yaml").exists()
    assert (tmp_path / "structures" / "b.yaml").exists()


def test_network2fnf_cli_structures_renames_duplicate_filenames(tmp_path):
    states_1 = tmp_path / "states_1"
    states_2 = tmp_path / "states_2"
    states_1.mkdir()
    states_2.mkdir()

    (states_1 / "a.yaml").write_text(
        """
energy:
  electronic: -1.0
vibrations:
  frequencies_cm-1: [100.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_2 / "a.yaml").write_text(
        """
energy:
  electronic: -0.5
vibrations:
  frequencies_cm-1: [200.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "ts.yaml").write_text(
        """
energy:
  electronic: 0.3
vibrations:
  frequencies_cm-1: [300.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A1*
    file: states_1/a.yaml
  - name: A2*
    file: states_2/a.yaml
transition_states:
  - name: TS
    file: ts.yaml
network:
  - name: surf step
    type: surf
    reaction: A1* => A2*
    forward:
      ts:
        - name: TS
      is:
        - name: A1*
    backward:
      ts:
        - name: TS
      is:
        - name: A2*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["network2fnf", str(network_file), "-o", str(output), "--structures", "structures"],
    )

    assert result.exit_code == 0
    payload = yaml.safe_load(output.read_text())

    assert payload["nodes"] == [
        {"label": "A1*", "structure": "structures/a__1.yaml"},
        {"label": "A2*", "structure": "structures/a__2.yaml"},
    ]
    assert (tmp_path / "structures" / "a__1.yaml").exists()
    assert (tmp_path / "structures" / "a__2.yaml").exists()


def test_network2fnf_cli_prune_removes_matching_two_node_steps(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    for name, energy in (("a.yaml", -1.0), ("b.yaml", -0.5), ("c.yaml", -0.2), ("ts.yaml", 0.4)):
        (states_dir / name).write_text(
            f"""
energy:
  electronic: {energy}
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
            + "\n",
            encoding="utf-8",
        )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A*
    file: states/a.yaml
  - name: B*
    file: states/b.yaml
  - name: C*
    file: states/c.yaml
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: step ab
    type: surf
    reaction: A* => B*
    forward:
      ts:
        - name: TS
      is:
        - name: A*
    backward:
      ts:
        - name: TS
      is:
        - name: B*
  - name: step bc
    type: surf
    reaction: B* => C*
    forward:
      ts:
        - name: TS
      is:
        - name: B*
    backward:
      ts:
        - name: TS
      is:
        - name: C*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["network2fnf", str(network_file), "-o", str(output), "--prune", "B*"],
    )

    assert result.exit_code == 0
    assert "Removed elementary reaction steps via --prune" in result.output
    assert "step ab (A* <-> B*)" in result.output
    assert "step bc (B* <-> C*)" in result.output

    payload = yaml.safe_load(output.read_text())
    assert payload["edges"] == []


def test_network2fnf_cli_supports_kj_mol(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    (states_dir / "is.yaml").write_text(
        """
energy:
  electronic: -1.0
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_dir / "fs.yaml").write_text(
        """
energy:
  electronic: 0.0
vibrations:
  frequencies_cm-1: []
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
network:
  - name: ads step
    type: ads
    reaction: IS => FS
    is:
      - name: IS
    fs:
      - name: FS
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf_kj.yaml"
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["network2fnf", str(network_file), "-o", str(output), "--unit", "kj/mol"],
    )

    assert result.exit_code == 0
    payload = yaml.safe_load(output.read_text())
    assert payload["pymkmkit"]["units"] == "kJ/mol"
    assert payload["edges"][0]["ads"] == 96.48533212


def test_network2fnf_cli_ignores_gas_and_warns_for_multinode_steps(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    for name, energy in (("a.yaml", -1.0), ("b.yaml", -0.5), ("c.yaml", -0.2), ("g.yaml", -0.1), ("ts.yaml", 0.4)):
        (states_dir / name).write_text(
            f"""
energy:
  electronic: {energy}
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
            + "\n",
            encoding="utf-8",
        )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A*
    file: states/a.yaml
    type: surf
  - name: B*
    file: states/b.yaml
    type: surf
  - name: C*
    file: states/c.yaml
    type: surf
  - name: G
    file: states/g.yaml
    type: gas
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: ads step
    type: ads
    reaction: G + A* => B*
    is:
      - name: G
      - name: A*
    fs:
      - name: B*
  - name: bimolecular surf
    type: surf
    reaction: A* + C* => B*
    forward:
      ts:
        - name: TS
      is:
        - name: A*
        - name: C*
    backward:
      ts:
        - name: TS
      is:
        - name: B*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    runner = CliRunner()
    result = runner.invoke(cli, ["network2fnf", str(network_file), "-o", str(output)])

    assert result.exit_code == 0
    assert "WARNING" in result.output
    assert "bimolecular surf" in result.output
    payload = yaml.safe_load(output.read_text())
    assert payload["edges"][0]["nodes"] == ["A*", "B*"]
    output_lines = output.read_text().splitlines()
    warning_lines = [line.strip() for line in output_lines if "# WARNING:" in line]
    assert warning_lines == [
        "# WARNING: Elementary reaction step 'bimolecular surf' has more than two nodes. "
        "Remove this step or prune the number of nodes."
    ]
    warning_idx = next(i for i, line in enumerate(output_lines) if "# WARNING:" in line)
    step_idx = next(i for i, line in enumerate(output_lines) if "- name: bimolecular surf" in line)
    assert warning_idx == step_idx - 1


def test_network2fnf_cli_split_splits_supported_multinode_steps(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    for name, energy in (("a.yaml", -1.0), ("b.yaml", -0.5), ("c.yaml", -0.2), ("ts.yaml", 0.4)):
        (states_dir / name).write_text(
            f"""
energy:
  electronic: {energy}
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
            + "\n",
            encoding="utf-8",
        )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A*
    file: states/a.yaml
    type: surf
  - name: B*
    file: states/b.yaml
    type: surf
  - name: C*
    file: states/c.yaml
    type: surf
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: bimolecular surf
    type: surf
    reaction: A* + C* => B*
    forward:
      ts:
        - name: TS
      is:
        - name: A*
        - name: C*
    backward:
      ts:
        - name: TS
      is:
        - name: B*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    runner = CliRunner()
    result = runner.invoke(cli, ["network2fnf", str(network_file), "-o", str(output), "--split"])

    assert result.exit_code == 0
    assert "WARNING" not in result.output
    payload = yaml.safe_load(output.read_text())
    split_edges = [edge for edge in payload["edges"] if edge["name"] == "bimolecular surf"]
    assert len(split_edges) == 2
    assert split_edges[0]["nodes"] == ["A*", "B*"]
    assert split_edges[1]["nodes"] == ["C*", "B*"]
    assert "# WARNING:" not in output.read_text()


def test_network2fnf_cli_split_still_warns_for_unsupported_multinode_steps(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    for name, energy in (("a.yaml", -1.0), ("b.yaml", -0.5), ("c.yaml", -0.2), ("d.yaml", -0.1), ("ts.yaml", 0.4)):
        (states_dir / name).write_text(
            f"""
energy:
  electronic: {energy}
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
            + "\n",
            encoding="utf-8",
        )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A*
    file: states/a.yaml
    type: surf
  - name: B*
    file: states/b.yaml
    type: surf
  - name: C*
    file: states/c.yaml
    type: surf
  - name: D*
    file: states/d.yaml
    type: surf
transition_states:
  - name: TS
    file: states/ts.yaml
network:
  - name: unsupported multinode surf
    type: surf
    reaction: A* + C* => B* + D*
    forward:
      ts:
        - name: TS
      is:
        - name: A*
        - name: C*
    backward:
      ts:
        - name: TS
      is:
        - name: B*
        - name: D*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    runner = CliRunner()
    result = runner.invoke(cli, ["network2fnf", str(network_file), "-o", str(output), "--split"])

    assert result.exit_code == 0
    assert "WARNING" in result.output
    assert "unsupported multinode surf" in result.output
    payload = yaml.safe_load(output.read_text())
    assert payload["edges"][0]["nodes"] == ["A*", "C*", "B*", "D*"]
    assert "# WARNING:" in output.read_text()

def test_read_network_rejects_mismatched_single_transition_states(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    for name, energy in (("is.yaml", -1.0), ("fs.yaml", -2.0), ("ts_a.yaml", 0.0), ("ts_b.yaml", 0.2)):
        (states_dir / name).write_text(
            f"""
energy:
  electronic: {energy}
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
            + "\n",
            encoding="utf-8",
        )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: IS*
    file: states/is.yaml
  - name: FS*
    file: states/fs.yaml
transition_states:
  - name: TS_A
    file: states/ts_a.yaml
  - name: TS_B
    file: states/ts_b.yaml
network:
  - name: mismatched ts step
    type: surf
    reaction: IS* => FS*
    forward:
      ts:
        - name: TS_A
      is:
        - name: IS*
    backward:
      ts:
        - name: TS_B
      is:
        - name: FS*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="must use the same TS state"):
        read_network(network_file)


def test_read_network_supports_rearrangement_steps(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    (states_dir / "is.yaml").write_text(
        """
energy:
  electronic: -10.0
vibrations:
  frequencies_cm-1: [1000.0, 2000.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_dir / "fs.yaml").write_text(
        """
energy:
  electronic: -9.5
vibrations:
  frequencies_cm-1: [1200.0, 2100.0]
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: CO_2O*
    file: states/is.yaml
  - name: CO2_O*
    file: states/fs.yaml
network:
  - name: CO oxidation via LH1
    type: rearrangement
    reaction: CO_2O* => CO2_O*
    is:
      - name: CO_2O*
    fs:
      - name: CO2_O*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    steps = read_network(network_file)

    assert len(steps) == 1
    step = steps[0]
    assert step.step_type == "rearrangement"
    assert step.forward_barrier_electronic == pytest.approx(0.5)
    assert step.reverse_barrier_electronic == pytest.approx(-0.5)
    assert step.forward_zpe_correction == pytest.approx(0.5 * 300.0 * 1.239841984e-4)
    assert step.reverse_zpe_correction == pytest.approx(-0.5 * 300.0 * 1.239841984e-4)
    assert step.forward_total_barrier == pytest.approx(
        0.5 + 0.5 * 300.0 * 1.239841984e-4
    )
    assert step.reverse_total_barrier == pytest.approx(
        -(0.5 + 0.5 * 300.0 * 1.239841984e-4)
    )
    assert step.reaction_heat_total == pytest.approx(step.forward_total_barrier)




def test_read_network_cli_displays_rearrangement_type(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    (states_dir / "is.yaml").write_text(
        """
energy:
  electronic: -10.0
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_dir / "fs.yaml").write_text(
        """
energy:
  electronic: -9.0
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: CO_2O*
    file: states/is.yaml
  - name: CO2_O*
    file: states/fs.yaml
network:
  - name: CO oxidation via LH1
    type: rearrangement
    reaction: CO_2O* => CO2_O*
    is:
      - name: CO_2O*
    fs:
      - name: CO2_O*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["read_network", str(network_file)])

    assert result.exit_code == 0
    assert "Reaction: CO_2O* => CO2_O* (rearrangement)" in result.output
    assert "Forward barrier:" in result.output
    assert "Reverse barrier:" in result.output

def test_network2fnf_cli_supports_rearrangement_edges(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    (states_dir / "is.yaml").write_text(
        """
energy:
  electronic: -2.0
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_dir / "fs.yaml").write_text(
        """
energy:
  electronic: -1.7
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A*
    file: states/is.yaml
    type: surf
  - name: B*
    file: states/fs.yaml
    type: surf
network:
  - name: rearrange
    type: rearrangement
    reaction: A* => B*
    is:
      - name: A*
    fs:
      - name: B*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    runner = CliRunner()
    result = runner.invoke(cli, ["network2fnf", str(network_file), "-o", str(output)])

    assert result.exit_code == 0
    payload = yaml.safe_load(output.read_text())
    edge = payload["edges"][0]
    assert edge["type"] == "rearrangement"
    assert edge["nodes"] == ["A*", "B*"]
    assert edge["forward"] == pytest.approx(0.3)
    assert edge["backward"] == pytest.approx(-0.3)


def test_network2fnf_cli_fails_if_output_exists_without_merge(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    (states_dir / "a.yaml").write_text(
        """
energy:
  electronic: -1.0
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (states_dir / "b.yaml").write_text(
        """
energy:
  electronic: -0.8
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
        + "\n",
        encoding="utf-8",
    )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A*
    file: states/a.yaml
  - name: B*
    file: states/b.yaml
network:
  - name: ads step
    type: ads
    reaction: A* => B*
    is:
      - name: A*
    fs:
      - name: B*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    output.write_text("nodes: []\nedges: []\n", encoding="utf-8")

    runner = CliRunner()
    result = runner.invoke(cli, ["network2fnf", str(network_file), "-o", str(output)])

    assert result.exit_code != 0
    assert "Output file already exists" in result.output


def test_network2fnf_cli_merge_adds_and_skips_entries(tmp_path):
    states_dir = tmp_path / "states"
    states_dir.mkdir()

    for name, energy in (("a.yaml", -1.0), ("b.yaml", -0.6), ("c.yaml", -0.3)):
        (states_dir / name).write_text(
            f"""
energy:
  electronic: {energy}
vibrations:
  frequencies_cm-1: []
  paired_modes_averaged: true
""".strip()
            + "\n",
            encoding="utf-8",
        )

    network_file = tmp_path / "network.yaml"
    network_file.write_text(
        """
stable_states:
  - name: A*
    file: states/a.yaml
  - name: B*
    file: states/b.yaml
  - name: C*
    file: states/c.yaml
network:
  - name: skip me
    type: ads
    reaction: A* => B*
    is:
      - name: A*
    fs:
      - name: B*
  - name: add me
    type: ads
    reaction: B* => C*
    is:
      - name: B*
    fs:
      - name: C*
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "fnf.yaml"
    output.write_text(
        """
pymkmkit:
  version: 0.1.0
  units: eV
  energy_type: elec+zpe
nodes:
  - label: A*
  - label: B*
edges:
  - name: existing AB
    type: ads
    nodes: [B*, A*]
    ads: 0.1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["network2fnf", str(network_file), "-o", str(output), "--merge"],
    )

    assert result.exit_code == 0
    assert "Added nodes:" in result.output
    assert "  + C*" in result.output
    assert "Skipped nodes:" in result.output
    assert "  - A*" in result.output
    assert "  - B*" in result.output
    assert "Added elementary reaction steps:" in result.output
    assert "  + add me (B* <-> C*)" in result.output
    assert "Skipped elementary reaction steps:" in result.output
    assert "  - skip me (A* <-> B*)" in result.output

    payload = yaml.safe_load(output.read_text())
    assert payload["nodes"] == [{"label": "A*"}, {"label": "B*"}, {"label": "C*"}]
    assert [edge["name"] for edge in payload["edges"]] == ["existing AB", "add me"]
