from pathlib import Path

import click

from pymkmkit.network_reader import build_ped, evaluate_paths, read_network
from pymkmkit.vasp_freq import parse_vasp_frequency, parse_vasp_optimization
from pymkmkit.yaml_writer import write_yaml


@click.group()
def cli():
    """pymkmkit: DFT → microkinetic data tools"""
    pass


# ---------------------------------------------------------------------------
# Shared CLI utilities
# ---------------------------------------------------------------------------
def _ensure_output_dir(output):
    """Create parent directories for an output file path.

    Parameters
    ----------
    output : str | pathlib.Path
        Destination file path provided by the CLI command.
    """
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# OUTCAR parsers
# ---------------------------------------------------------------------------
@cli.command("freq2yaml")
@click.argument(
    "outcar",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(dir_okay=False),
    help="Output YAML file (required)."
)
@click.option(
    "--average-pairs",
    is_flag=True,
    help=(
        "Average vibrational modes in sequential pairs. "
        "Use when two identical adsorbates are present."
    ),
)
def freq2yaml(outcar, output, average_pairs):
    """Convert a VASP frequency-calculation ``OUTCAR`` file to YAML.

    Parameters
    ----------
    outcar : str
        Input OUTCAR path.
    output : str
        Output YAML path.
    average_pairs : bool
        Whether sequential vibrational modes should be pair-averaged.
    """

    _ensure_output_dir(output)

    data = parse_vasp_frequency(
        outcar,
        average_pairs=average_pairs,
    )

    write_yaml(data, output)

    click.echo(f"YAML written to: {output}")

    if average_pairs:
        click.echo("Sequential mode pairs averaged.")


@cli.command("opt2yaml")
@click.argument(
    "outcar",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(dir_okay=False),
    help="Output YAML file (required)."
)
def opt2yaml(outcar, output):
    """Convert a VASP optimization ``OUTCAR`` file to YAML.

    Parameters
    ----------
    outcar : str
        Input OUTCAR path.
    output : str
        Output YAML path.
    """

    _ensure_output_dir(output)

    data = parse_vasp_optimization(outcar)

    write_yaml(data, output)

    click.echo(f"YAML written to: {output}")


# ---------------------------------------------------------------------------
# Network analysis
# ---------------------------------------------------------------------------
@cli.command("read_network")
@click.argument(
    "network_file",
    type=click.Path(exists=True, dir_okay=False)
)
def read_network_command(network_file):
    """Load a network definition and print per-step barrier information.

    Parameters
    ----------
    network_file : str
        Path to the network YAML file.
    """

    steps = read_network(network_file)

    for step in steps:
        click.echo(f"Reaction: {step.reaction}")
        if step.step_type == "ads":
            click.echo(
                "  Adsorption heat: "
                f"{step.forward_total_barrier:.6f} eV "
                f"(inc. ZPE-corr: {step.forward_zpe_correction:.6f})"
            )
            continue

        click.echo(
            "  Forward barrier: "
            f"{step.forward_total_barrier:.6f} eV "
            f"(inc. ZPE-corr: {step.forward_zpe_correction:.6f})"
        )
        click.echo(
            "  Reverse barrier: "
            f"{step.reverse_total_barrier:.6f} eV "
            f"(inc. ZPE-corr: {step.reverse_zpe_correction:.6f})"
        )


@cli.command("evaluate_paths")
@click.argument(
    "network_file",
    type=click.Path(exists=True, dir_okay=False)
)
def evaluate_paths_command(network_file):
    """Evaluate and print total reaction energy for each named path.

    Parameters
    ----------
    network_file : str
        Path to the network YAML file.
    """

    paths = evaluate_paths(network_file)

    for path in paths:
        click.echo(f"Path {path.name}: {path.total_reaction_energy:.6f} eV")


@cli.command("build_ped")
@click.argument(
    "network_file",
    type=click.Path(exists=True, dir_okay=False)
)
@click.argument("path_name")
@click.argument(
    "output_file",
    required=False,
    type=click.Path(dir_okay=False),
)
def build_ped_command(network_file, path_name, output_file):
    """Generate a potential energy diagram for a selected reaction path.

    Parameters
    ----------
    network_file : str
        Path to the network YAML file.
    path_name : str
        Name of the path to render.
    output_file : str | None
        Optional output image path. If omitted, a GUI window is shown.
    """

    build_ped(network_file, path_name, output_file)

    if output_file:
        click.echo(f"Potential energy diagram written to: {output_file}")
