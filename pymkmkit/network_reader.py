from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

EV_PER_CM1 = 1.239841984e-4


@dataclass(frozen=True)
class State:
    name: str
    file: Path
    electronic_energy: float
    zpe_energy: float
    paired_modes_averaged: bool


@dataclass(frozen=True)
class ElementaryStep:
    name: str
    step_type: str
    reaction: str
    forward_equation: str
    reverse_equation: str
    forward_barrier_electronic: float
    reverse_barrier_electronic: float
    forward_zpe_correction: float
    reverse_zpe_correction: float
    forward_total_barrier: float
    reverse_total_barrier: float


def _resolve_state_file(base_dir: Path, state_file: str) -> Path:
    path = (base_dir / state_file).resolve()
    if path.exists():
        return path

    yaml_path = path.with_suffix(".yaml")
    if yaml_path.exists():
        return yaml_path

    yml_path = path.with_suffix(".yml")
    if yml_path.exists():
        return yml_path

    raise FileNotFoundError(f"State file not found: {state_file}")


def _read_state_data(state_path: Path) -> tuple[float, float, bool]:
    with state_path.open("r", encoding="utf-8") as stream:
        state_data = yaml.safe_load(stream) or {}

    try:
        electronic_energy = float(state_data["energy"]["electronic"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            f"State file '{state_path}' is missing a valid energy.electronic value"
        ) from exc

    vibrations = state_data.get("vibrations", {})
    frequencies = vibrations.get("frequencies_cm-1", [])

    if not isinstance(frequencies, list):
        raise ValueError(
            f"State file '{state_path}' has invalid vibrations.frequencies_cm-1 format"
        )

    try:
        frequencies_cm1 = [float(freq) for freq in frequencies]
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"State file '{state_path}' has non-numeric frequency values"
        ) from exc

    zpe_energy = 0.5 * sum(frequencies_cm1) * EV_PER_CM1

    paired_modes_averaged = bool(vibrations.get("paired_modes_averaged", False))

    return electronic_energy, zpe_energy, paired_modes_averaged


def _load_states(network_data: dict, base_dir: Path) -> dict[str, State]:
    states: dict[str, State] = {}

    for section in ("stable_states", "transition_states"):
        for state in network_data.get(section, []):
            name = state.get("name")
            state_file = state.get("file")
            if not name or not state_file:
                raise ValueError(f"Invalid state entry in '{section}': {state}")

            resolved_file = _resolve_state_file(base_dir, state_file)
            electronic_energy, zpe_energy, paired = _read_state_data(resolved_file)
            states[name] = State(
                name=name,
                file=resolved_file,
                electronic_energy=electronic_energy,
                zpe_energy=zpe_energy,
                paired_modes_averaged=paired,
            )

    return states


def _sum_energy(terms: list[dict], states: dict[str, State]) -> tuple[float, str]:
    if not terms:
        return 0.0, "0"

    total = 0.0
    pieces: list[str] = []

    for term in terms:
        name = term.get("name")
        stoich = term.get("stoichiometry", 1)

        if name not in states:
            raise ValueError(f"Unknown state '{name}' referenced in network step")

        try:
            stoich_value = float(stoich)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid stoichiometry for state '{name}': {stoich}") from exc

        total += stoich_value * states[name].electronic_energy
        pieces.append(f"{stoich_value:g}*E({name})")

    return total, " + ".join(pieces)


def _sum_zpe(
    terms: list[dict],
    states: dict[str, State],
    normalization_value: float,
) -> tuple[float, str]:
    if not terms:
        return 0.0, "0"

    total = 0.0
    pieces: list[str] = []

    for term in terms:
        name = term.get("name")
        stoich = term.get("stoichiometry", 1)

        if name not in states:
            raise ValueError(f"Unknown state '{name}' referenced in network step")

        try:
            stoich_value = float(stoich)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid stoichiometry for state '{name}': {stoich}") from exc

        zpe_term = stoich_value * states[name].zpe_energy
        if states[name].paired_modes_averaged:
            pieces.append(f"{stoich_value:g}*ZPE({name})")
        else:
            zpe_term /= normalization_value
            pieces.append(f"({stoich_value:g}*ZPE({name})/{normalization_value:g})")

        total += zpe_term

    return total, " + ".join(pieces)


def _compute_barrier(
    direction_data: dict,
    states: dict[str, State],
) -> tuple[float, float, float, str]:
    ts_energy, ts_energy_expr = _sum_energy(direction_data.get("ts", []), states)
    is_energy, is_energy_expr = _sum_energy(direction_data.get("is", []), states)

    normalization = direction_data.get("normalization", 1)
    try:
        normalization_value = float(normalization)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid normalization value: {normalization}") from exc

    if normalization_value == 0:
        raise ValueError("Normalization cannot be zero")

    electronic_barrier = (ts_energy - is_energy) / normalization_value
    electronic_equation = (
        f"({ts_energy_expr} - ({is_energy_expr})) / {normalization_value:g}"
    )

    ts_zpe, ts_zpe_expr = _sum_zpe(direction_data.get("ts", []), states, normalization_value)
    is_zpe, is_zpe_expr = _sum_zpe(direction_data.get("is", []), states, normalization_value)

    zpe_correction = ts_zpe - is_zpe
    zpe_equation = f"({ts_zpe_expr} - ({is_zpe_expr}))"

    total_barrier = electronic_barrier + zpe_correction
    full_equation = (
        f"E_el: {electronic_equation}; "
        f"ZPE corr: {zpe_equation}; "
        f"Total: ({electronic_equation}) + ({zpe_equation})"
    )

    return electronic_barrier, zpe_correction, total_barrier, full_equation


def _compute_adsorption_heat(
    step_data: dict,
    states: dict[str, State],
) -> tuple[float, float, float, str]:
    fs_energy, fs_energy_expr = _sum_energy(step_data.get("fs", []), states)
    is_energy, is_energy_expr = _sum_energy(step_data.get("is", []), states)

    normalization = step_data.get("normalization", 1)
    try:
        normalization_value = float(normalization)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid normalization value: {normalization}") from exc

    if normalization_value == 0:
        raise ValueError("Normalization cannot be zero")

    electronic_heat = (fs_energy - is_energy) / normalization_value
    electronic_equation = (
        f"({fs_energy_expr} - ({is_energy_expr})) / {normalization_value:g}"
    )

    fs_zpe, fs_zpe_expr = _sum_zpe(step_data.get("fs", []), states, normalization_value)
    is_zpe, is_zpe_expr = _sum_zpe(step_data.get("is", []), states, normalization_value)

    zpe_correction = fs_zpe - is_zpe
    zpe_equation = f"({fs_zpe_expr} - ({is_zpe_expr}))"

    total_heat = electronic_heat + zpe_correction
    full_equation = (
        f"E_el: {electronic_equation}; "
        f"ZPE corr: {zpe_equation}; "
        f"Total: ({electronic_equation}) + ({zpe_equation})"
    )

    return electronic_heat, zpe_correction, total_heat, full_equation


def read_network(network_file: str | Path) -> list[ElementaryStep]:
    network_path = Path(network_file).resolve()

    with network_path.open("r", encoding="utf-8") as stream:
        network_data = yaml.safe_load(stream) or {}

    states = _load_states(network_data, network_path.parent)

    steps: list[ElementaryStep] = []
    for step in network_data.get("network", []):
        step_type = step.get("type", "surf")
        if step_type not in {"surf", "ads"}:
            raise ValueError(
                f"Invalid network step type '{step_type}'. Expected 'surf' or 'ads'."
            )

        if step_type == "ads":
            (
                forward_el,
                forward_zpe,
                forward_total,
                forward_eq,
            ) = _compute_adsorption_heat(step, states)
            reverse_el = 0.0
            reverse_zpe = 0.0
            reverse_total = 0.0
            reverse_eq = "N/A"
        else:
            forward_data = step.get("forward", {})
            backward_data = step.get("backward", {})

            (
                forward_el,
                forward_zpe,
                forward_total,
                forward_eq,
            ) = _compute_barrier(forward_data, states)
            (
                reverse_el,
                reverse_zpe,
                reverse_total,
                reverse_eq,
            ) = _compute_barrier(backward_data, states)

        steps.append(
            ElementaryStep(
                name=step.get("name", "unnamed_step"),
                step_type=step_type,
                reaction=step.get("reaction", ""),
                forward_equation=forward_eq,
                reverse_equation=reverse_eq,
                forward_barrier_electronic=forward_el,
                reverse_barrier_electronic=reverse_el,
                forward_zpe_correction=forward_zpe,
                reverse_zpe_correction=reverse_zpe,
                forward_total_barrier=forward_total,
                reverse_total_barrier=reverse_total,
            )
        )

    return steps
