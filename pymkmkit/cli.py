from pathlib import Path

import click

from pymkmkit.vasp_freq import parse_vasp_frequency, parse_vasp_optimization
from pymkmkit.yaml_writer import write_yaml


@click.group()
def cli():
    """pymkmkit: DFT → microkinetic data tools"""
    pass


def _ensure_output_dir(output):
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)


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
    """Parse a VASP OUTCAR frequency calculation into YAML."""

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
    """Parse a VASP OUTCAR geometry optimization into YAML."""

    _ensure_output_dir(output)

    data = parse_vasp_optimization(outcar)

    write_yaml(data, output)

    click.echo(f"YAML written to: {output}")
