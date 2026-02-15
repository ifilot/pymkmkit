from pathlib import Path
import click

from pymkmkit.vasp_freq import parse_vasp_frequency
from pymkmkit.yaml_writer import write_yaml


@click.group()
def cli():
    """pymkmkit: DFT → microkinetic data tools"""
    pass


@cli.command()
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
def outcar2yaml(outcar, output, average_pairs):
    """
    Parse a VASP OUTCAR frequency calculation into YAML.
    """

    output_path = Path(output)

    # create directory if needed
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = parse_vasp_frequency(
        outcar,
        average_pairs=average_pairs,
    )

    write_yaml(data, output)

    click.echo(f"YAML written to: {output}")

    if average_pairs:
        click.echo("Sequential mode pairs averaged.")
