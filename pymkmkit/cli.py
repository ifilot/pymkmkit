from pathlib import Path

import click
import yaml
from rich.console import Console

from pymkmkit.network_reader import build_fnf, build_ped, evaluate_paths, read_network
from pymkmkit.vasp_freq import (
    parse_ase_vibrations,
    parse_vasp_frequency,
    parse_vasp_optimization,
)
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


def _format_energy(value: float, *, unit: str) -> str:
    """Format an energy value with unit-specific precision and alignment."""
    if unit == "kj/mol":
        converted = value * 96.48533212
        return f"{converted:>8.1f} kJ/mol"
    return f"{value:>8.4f} eV"


def _dump_yaml(data: dict, output_file: str) -> None:
    """Write dictionary content to a YAML file."""
    with Path(output_file).open("w", encoding="utf-8") as stream:
        yaml.safe_dump(data, stream, sort_keys=False)


def _edge_node_key(edge: dict) -> tuple[str, ...]:
    """Return an order-insensitive edge key from its node labels."""
    nodes = edge.get("nodes", [])
    if not isinstance(nodes, list) or len(nodes) != 2:
        raise ValueError(f"Invalid edge nodes definition: {nodes}")
    return tuple(sorted(str(node) for node in nodes))


def _merge_fnf_payloads(existing_payload: dict, new_payload: dict) -> tuple[dict, dict]:
    """Merge newly generated FNF data into an existing FNF payload.

    Returns
    -------
    tuple[dict, dict]
        Merged payload and a report dictionary with added/skipped nodes and edges.
    """
    merged = dict(existing_payload)
    merged.setdefault("pymkmkit", new_payload.get("pymkmkit", {}))
    merged_nodes = list(existing_payload.get("nodes", []))
    merged_edges = list(existing_payload.get("edges", []))

    existing_node_labels = {
        str(node.get("label"))
        for node in merged_nodes
        if isinstance(node, dict) and node.get("label") is not None
    }
    added_nodes: list[str] = []
    skipped_nodes: list[str] = []
    for node in new_payload.get("nodes", []):
        label = str(node.get("label"))
        if label in existing_node_labels:
            skipped_nodes.append(label)
            continue
        merged_nodes.append(node)
        existing_node_labels.add(label)
        added_nodes.append(label)

    existing_edge_keys = {_edge_node_key(edge) for edge in merged_edges}
    added_edges: list[str] = []
    skipped_edges: list[str] = []
    for edge in new_payload.get("edges", []):
        edge_key = _edge_node_key(edge)
        edge_name = str(edge.get("name", "unnamed_step"))
        descriptor = f"{edge_name} ({edge_key[0]} <-> {edge_key[1]})"
        if edge_key in existing_edge_keys:
            skipped_edges.append(descriptor)
            continue
        merged_edges.append(edge)
        existing_edge_keys.add(edge_key)
        added_edges.append(descriptor)

    merged["nodes"] = merged_nodes
    merged["edges"] = merged_edges
    return merged, {
        "added_nodes": added_nodes,
        "skipped_nodes": skipped_nodes,
        "added_edges": added_edges,
        "skipped_edges": skipped_edges,
    }


def _print_merge_report(report: dict) -> None:
    """Print a colored merge report for nodes and elementary steps."""
    console = Console()

    console.print("[green]Added nodes:[/green]")
    if report["added_nodes"]:
        for node in report["added_nodes"]:
            console.print(f"  [green]+[/green] {node}")
    else:
        console.print("  [dim]none[/dim]")

    console.print("[yellow]Skipped nodes:[/yellow]")
    if report["skipped_nodes"]:
        for node in report["skipped_nodes"]:
            console.print(f"  [yellow]-[/yellow] {node}")
    else:
        console.print("  [dim]none[/dim]")

    console.print("[green]Added elementary reaction steps:[/green]")
    if report["added_edges"]:
        for edge in report["added_edges"]:
            console.print(f"  [green]+[/green] {edge}")
    else:
        console.print("  [dim]none[/dim]")

    console.print("[yellow]Skipped elementary reaction steps:[/yellow]")
    if report["skipped_edges"]:
        for edge in report["skipped_edges"]:
            console.print(f"  [yellow]-[/yellow] {edge}")
    else:
        console.print("  [dim]none[/dim]")


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


@cli.command("asevib2yaml")
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
def asevib2yaml(outcar, output):
    """Convert OUTCAR + sibling ``vibX`` ASE caches to a single YAML file."""

    _ensure_output_dir(output)

    data = parse_ase_vibrations(outcar)

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
@click.option(
    "--unit",
    type=click.Choice(["ev", "kj/mol"], case_sensitive=False),
    default="ev",
    show_default=True,
    help="Energy unit for displayed values.",
)
def read_network_command(network_file, unit):
    """Load a network definition and print per-step barrier information.

    Parameters
    ----------
    network_file : str
        Path to the network YAML file.
    """

    steps = read_network(network_file)
    console = Console()

    for step in steps:
        console.print(f"Reaction: [cyan]{step.reaction}[/cyan] [plum1]({step.step_type})[/plum1]")
        if step.step_type == "ads":
            console.print(
                "  Adsorption heat: "
                f"[green]{_format_energy(step.forward_total_barrier, unit=unit)}[/green] "
                "(inc. ZPE-corr: "
                f"[green]{_format_energy(step.forward_zpe_correction, unit=unit)}[/green])"
            )
            continue

        console.print(
            "  Forward barrier: "
            f"[green]{_format_energy(step.forward_total_barrier, unit=unit)}[/green] "
            "(inc. ZPE-corr: "
            f"[green]{_format_energy(step.forward_zpe_correction, unit=unit)}[/green])"
        )
        console.print(
            "  Reverse barrier: "
            f"[green]{_format_energy(step.reverse_total_barrier, unit=unit)}[/green] "
            "(inc. ZPE-corr: "
            f"[green]{_format_energy(step.reverse_zpe_correction, unit=unit)}[/green])"
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


@cli.command("network2fnf")
@click.argument(
    "network_file",
    type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(dir_okay=False),
    help="Output FNF YAML file (required)."
)
@click.option(
    "--unit",
    type=click.Choice(["ev", "kj/mol"], case_sensitive=False),
    default="ev",
    show_default=True,
    help="Energy unit used for edge values.",
)
@click.option(
    "--merge",
    is_flag=True,
    help=(
        "Merge into an existing FNF file. Existing nodes are matched by label and "
        "existing elementary steps are matched by the combination of nodes (order-insensitive)."
    ),
)
def network2fnf_command(network_file, output, unit, merge):
    """Convert a network YAML definition to a formatted network file (FNF)."""
    _ensure_output_dir(output)
    output_path = Path(output)

    if output_path.exists() and not merge:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --merge to append new entries."
        )

    fnf_payload = build_fnf(network_file, unit=unit)

    if merge and output_path.exists():
        with output_path.open("r", encoding="utf-8") as stream:
            existing_payload = yaml.safe_load(stream) or {}
        merged_payload, report = _merge_fnf_payloads(existing_payload, fnf_payload)
        _dump_yaml(merged_payload, output)
        _print_merge_report(report)
    else:
        _dump_yaml(fnf_payload, output)

    click.echo(f"YAML written to: {output}")
