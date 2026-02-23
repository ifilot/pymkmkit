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
    reaction_heat_electronic: float
    reaction_heat_zpe_correction: float
    reaction_heat_total: float


@dataclass(frozen=True)
class ReactionPath:
    name: str
    total_reaction_energy: float


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
            reaction_heat_el = forward_el
            reaction_heat_zpe = forward_zpe
            reaction_heat_total = forward_total
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
            reaction_heat_el = forward_el - reverse_el
            reaction_heat_zpe = forward_zpe - reverse_zpe
            reaction_heat_total = forward_total - reverse_total

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
                reaction_heat_electronic=reaction_heat_el,
                reaction_heat_zpe_correction=reaction_heat_zpe,
                reaction_heat_total=reaction_heat_total,
            )
        )

    return steps


def evaluate_paths(network_file: str | Path) -> list[ReactionPath]:
    network_path = Path(network_file).resolve()
    with network_path.open("r", encoding="utf-8") as stream:
        network_data = yaml.safe_load(stream) or {}

    steps = read_network(network_file)
    steps_by_name = {step.name: step for step in steps}

    paths: list[ReactionPath] = []
    for path in network_data.get("paths", []):
        path_name = path.get("name")
        if not path_name:
            raise ValueError(f"Invalid path entry without a name: {path}")

        total_reaction_energy = 0.0
        for path_step in path.get("steps", []):
            step_name = path_step.get("name")
            if step_name not in steps_by_name:
                raise ValueError(f"Unknown step '{step_name}' referenced in path '{path_name}'")

            factor = path_step.get("factor", 1)
            try:
                factor_value = float(factor)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid factor for step '{step_name}' in path '{path_name}': {factor}"
                ) from exc

            step_heat = steps_by_name[step_name].reaction_heat_total
            total_reaction_energy += factor_value * step_heat

        paths.append(
            ReactionPath(
                name=path_name,
                total_reaction_energy=total_reaction_energy,
            )
        )

    return paths


def build_ped(
    network_file: str | Path,
    path_name: str,
    output_file: str | Path | None = None,
) -> None:
    """Build a potential energy diagram (PED) for a named path."""

    import numpy as np
    import matplotlib.pyplot as plt

    network_path = Path(network_file).resolve()
    with network_path.open("r", encoding="utf-8") as stream:
        network_data = yaml.safe_load(stream) or {}

    steps = read_network(network_file)
    steps_by_name = {step.name: step for step in steps}

    selected_path = None
    for path in network_data.get("paths", []):
        if path.get("name") == path_name:
            selected_path = path
            break

    if selected_path is None:
        raise ValueError(f"Path '{path_name}' not found in network file")

    condensed_steps: list[tuple[ElementaryStep, float, str]] = []
    for index, path_step in enumerate(selected_path.get("steps", []), start=1):
        step_name = path_step.get("name")
        if step_name not in steps_by_name:
            raise ValueError(f"Unknown step '{step_name}' referenced in path '{path_name}'")

        factor = path_step.get("factor", 1)
        try:
            factor_value = float(factor)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid factor for step '{step_name}' in path '{path_name}': {factor}"
            ) from exc

        label = str(path_step.get("label", f"state_{index}"))
        condensed_steps.append((steps_by_name[step_name], factor_value, label))

    fig, ax = plt.subplots(figsize=(10, 7))

    stable_lw = 3.0
    connector_lw = 1.5
    label_offset = 0.0

    plateau_width = 0.35
    connector_width = 0.30
    x_cursor = 0.0
    current_energy = 0.0

    state_centers: list[float] = []
    state_energies: list[float] = []
    state_labels: list[str] = ["state_0"]

    # Always start from a horizontal line at the zero reference energy.
    ax.plot(
        [x_cursor, x_cursor + plateau_width],
        [current_energy, current_energy],
        color="black",
        linewidth=stable_lw,
    )
    state_centers.append(x_cursor + 0.5 * plateau_width)
    state_energies.append(current_energy)
    x_cursor += plateau_width

    for step, factor_value, label in condensed_steps:
        magnitude = abs(factor_value)
        direction = 1 if factor_value >= 0 else -1

        if direction > 0:
            final_energy = current_energy + magnitude * step.reaction_heat_total
            barrier = magnitude * step.forward_total_barrier
        else:
            final_energy = current_energy - magnitude * step.reaction_heat_total
            barrier = magnitude * step.reverse_total_barrier

        if step.step_type == "surf":
            ts_energy = current_energy + barrier
            x_start = x_cursor
            x_end = x_cursor + connector_width
            x_peak = 0.5 * (x_start + x_end)

            x_parabola = np.linspace(x_start, x_end, 120)
            a = np.array([
                [x_start**2, x_start, 1.0],
                [x_peak**2, x_peak, 1.0],
                [x_end**2, x_end, 1.0],
            ])
            b = np.array([current_energy, ts_energy, final_energy])
            coef = np.linalg.solve(a, b)
            y_parabola = coef[0] * x_parabola**2 + coef[1] * x_parabola + coef[2]
            ax.plot(
                x_parabola,
                y_parabola,
                color="tab:blue",
                linewidth=connector_lw,
                linestyle="--",
            )
        else:
            ax.plot(
                [x_cursor, x_cursor + connector_width],
                [current_energy, final_energy],
                color="tab:blue",
                linewidth=connector_lw,
                linestyle="--",
            )

        x_cursor += connector_width
        current_energy = final_energy

        # Draw stable-state plateau thicker than connecting lines.
        ax.plot(
            [x_cursor, x_cursor + plateau_width],
            [current_energy, current_energy],
            color="black",
            linewidth=stable_lw,
        )
        state_centers.append(x_cursor + 0.5 * plateau_width)
        state_energies.append(current_energy)
        state_labels.append(label)
        x_cursor += plateau_width

    # Draw state labels (-90°) with their top aligned to the state line.
    y_min = min(state_energies) if state_energies else 0.0
    y_max = max(state_energies) if state_energies else 0.0
    y_span = max(y_max - y_min, 1e-6)
    dy = label_offset * y_span
    for x_state, y_state, state_label in zip(state_centers, state_energies, state_labels):
        ax.text(
            x_state,
            y_state - dy,
            state_label,
            rotation=-90,
            rotation_mode="anchor",
            ha="left",
            va="top",
            fontsize=8,
            clip_on=False,
        )

    # Final reaction heat annotation with vertical arrow and vertical text.
    final_energy = state_energies[-1] if state_energies else 0.0
    arrow_x = x_cursor + 0.12
    ax.annotate(
        "",
        xy=(arrow_x, final_energy),
        xytext=(arrow_x, 0.0),
        arrowprops={"arrowstyle": "->", "linewidth": 1.5, "color": "black"},
    )

    delta_text = f"ΔE = {final_energy:.3f} eV"
    text_y = 0.5 * final_energy
    if abs(final_energy) < 1e-12:
        text_y = 0.0
    ax.text(
        arrow_x + 0.04,
        text_y,
        delta_text,
        rotation=-90,
        va="center",
        ha="left",
        fontsize=9,
    )

    ax.set_xlabel("Reaction coordinate")
    ax.set_ylabel("Energy (eV)")
    ax.set_title(f"Potential energy diagram: {path_name}")
    ax.set_xticks([])
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    y_min_plot = min(state_energies + [0.0])
    y_max_plot = max(state_energies + [0.0])
    y_span_plot = max(y_max_plot - y_min_plot, 1e-6)
    ax.set_ylim(y_min_plot - 0.45 * y_span_plot, y_max_plot + 0.20 * y_span_plot)
    fig.tight_layout()

    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300)
        plt.close(fig)
    else:
        plt.show()

