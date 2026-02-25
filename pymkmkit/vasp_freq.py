from pathlib import Path
import re
import numpy as np
from ase import Atoms
from ase.io import read
from ase.io.formats import UnknownFileTypeError
from ase.io import ParseError
from pymkmkit.yaml_writer import InlineList
from pymkmkit._version import get_version
from datetime import datetime, timezone

# ============================================================
# INCAR extraction (physics-relevant parameters only)
# ============================================================

IMPORTANT_KEYS = {
    "ENCUT": float,
    "PREC": str,
    "EDIFF": float,
    "EDIFFG": float,
    "ISMEAR": int,
    "SIGMA": float,
    "ISPIN": int,
    "IBRION": int,
    "POTIM": float,
    "NFREE": int,
    "ISIF": int,
    "NSW": int,
    "LREAL": str,
    "LASPH": str,
    "GGA": str,
    "METAGGA": str,
    "IVDW": int,
    "ALGO": str,
}


def _clean_value(value, cast):
    """Strip units/comments and cast to desired type."""
    value = value.split()[0].replace(";", "")
    try:
        return cast(value)
    except Exception:
        return value


def extract_incar_settings(text):
    """Extract selected INCAR parameters from OUTCAR text.

    Parameters
    ----------
    text : str
        Full OUTCAR content.

    Returns
    -------
    dict
        Mapping of key INCAR tags to parsed values.
    """
    incar = {}

    for line in text.splitlines():
        stripped = line.strip()

        for key, cast in IMPORTANT_KEYS.items():
            if stripped.startswith(key):
                parts = stripped.split("=")
                if len(parts) > 1:
                    incar[key] = _clean_value(parts[1], cast)

    return incar


# ============================================================
# POTCAR extraction (unique entries, preserve order)
# ============================================================

def extract_potcar_info(text):
    """Extract unique POTCAR descriptors from OUTCAR text."""
    pots = []

    for line in text.splitlines():
        if "POTCAR:" in line:
            entry = line.split("POTCAR:")[1].strip()
            if entry not in pots:
                pots.append(entry)

    return pots


# ============================================================
# Vibrational frequencies
# ============================================================



def extract_ionic_energies(text):
    """Extract ionic-step electronic energies from OUTCAR text."""
    energies = []

    for line in text.splitlines():
        if "energy  without entropy=" in line and "energy(sigma->0)" in line:
            try:
                energies.append(float(line.split("=")[-1].strip()))
            except ValueError:
                continue

    return energies


def extract_total_energies(text):
    """Extract free energies (TOTEN) from OUTCAR text."""
    energies = []

    for line in text.splitlines():
        if "free  energy   TOTEN" in line:
            try:
                energies.append(float(line.split("=")[-1].split()[0]))
            except ValueError:
                continue

    return energies


def _parse_atomic_symbols(text):
    """Infer atom symbols for each site from OUTCAR species metadata."""
    counts = None
    symbols = []

    for line in text.splitlines():
        if "ions per type" in line:
            counts = [int(x) for x in re.findall(r"\d+", line)]
        elif "TITEL" in line:
            parts = line.split()
            for part in reversed(parts):
                cleaned = part.strip("_.,")
                if cleaned and cleaned[0].isalpha():
                    symbol = ''.join(ch for ch in cleaned if ch.isalpha())
                    if symbol:
                        symbols.append(symbol.capitalize())
                    break

    if not counts or not symbols:
        raise ValueError("Could not parse species metadata from OUTCAR")

    unique_symbols = []
    for symbol in symbols:
        if symbol not in unique_symbols:
            unique_symbols.append(symbol)

    if len(unique_symbols) < len(counts):
        raise ValueError("Insufficient POTCAR species entries for ions per type")

    expanded = []
    for symbol, count in zip(unique_symbols, counts):
        expanded.extend([symbol] * count)

    return expanded


def _parse_last_lattice_vectors(lines):
    """Parse the final reported lattice vectors from OUTCAR lines."""
    for i in range(len(lines) - 1, -1, -1):
        if "direct lattice vectors" in lines[i]:
            vectors = []
            for j in range(i + 1, i + 4):
                parts = lines[j].split()
                vectors.append([float(parts[0]), float(parts[1]), float(parts[2])])
            return vectors

    raise ValueError("Could not parse lattice vectors from OUTCAR")


def _parse_last_direct_positions(lines, n_atoms):
    """Parse final direct (fractional) coordinates for all atoms."""
    marker = "position of ions in fractional coordinates (direct lattice)"

    for i in range(len(lines) - 1, -1, -1):
        if marker in lines[i]:
            positions = []
            for j in range(i + 1, i + 1 + n_atoms):
                parts = lines[j].split()
                positions.append([float(parts[0]), float(parts[1]), float(parts[2])])
            return positions

    raise ValueError("Could not parse direct coordinates from OUTCAR")


def _parse_atoms_from_outcar_text(text):
    """Construct an ASE ``Atoms`` object directly from OUTCAR text."""
    lines = text.splitlines()
    symbols = _parse_atomic_symbols(text)
    cell = _parse_last_lattice_vectors(lines)
    scaled_positions = _parse_last_direct_positions(lines, len(symbols))

    atoms = Atoms(symbols=symbols, cell=cell, pbc=True)
    atoms.set_scaled_positions(scaled_positions)
    return atoms


def _read_last_optimization_atoms(outcar_path, text):
    """Read final optimization geometry, with text-parsing fallback."""
    try:
        return read(outcar_path, index=-1)
    except (ParseError, UnknownFileTypeError, KeyError, ValueError):
        return _parse_atoms_from_outcar_text(text)

def extract_frequencies(text):
    """
    Extract vibrational frequencies from OUTCAR.

    VASP prints the frequency list twice. We detect the repetition by checking
    when the mode index restarts at 1 and stop reading.
    """
    real = []
    imaginary = []

    last_index = 0

    for line in text.splitlines():
        if " f  =" in line or " f/i=" in line:
            parts = line.split()

            try:
                index = int(parts[0])
            except ValueError:
                continue

            # detect repeated block (index resets)
            if index <= last_index:
                break

            last_index = index

            match = re.search(r"([-+]?\d*\.?\d+)\s+cm-1", line)
            if not match:
                continue

            value = float(match.group(1))

            if "f/i=" in line:
                imaginary.append(-abs(value))
            else:
                real.append(value)

    return real, imaginary


_DOF_LABEL_RE = re.compile(r"^(\d+)([XYZ])$")


def extract_perturbed_hessian(text):
    """Extract the perturbed Hessian block from OUTCAR.

    Returns
    -------
    tuple[list[str], list[list[float]]] | (None, None)
        DOF labels and Hessian matrix values for the perturbed coordinates.
    """
    lines = text.splitlines()
    start = None

    for idx, line in enumerate(lines):
        if "SECOND DERIVATIVES (NOT SYMMETRIZED)" in line:
            start = idx

    if start is None:
        return None, None

    row_order = []
    values = {}
    i = start + 1

    while i < len(lines):
        line = lines[i]

        if "Eigenvectors and eigenvalues of the dynamical matrix" in line:
            break

        labels = [tok for tok in line.split() if _DOF_LABEL_RE.match(tok)]
        if not labels:
            i += 1
            continue

        i += 1
        while i < len(lines):
            row_line = lines[i]

            if not row_line.strip():
                break

            parts = row_line.split()
            if not parts or not _DOF_LABEL_RE.match(parts[0]):
                break

            row_label = parts[0]
            numbers = parts[1:1 + len(labels)]

            if len(numbers) != len(labels):
                raise ValueError("Malformed Hessian block in OUTCAR")

            if row_label not in row_order:
                row_order.append(row_label)
                values[row_label] = {}

            for col_label, value in zip(labels, numbers):
                values[row_label][col_label] = float(value)

            i += 1

    if not row_order:
        return None, None

    matrix = [
        [values[row][col] for col in row_order]
        for row in row_order
    ]
    return row_order, matrix


def frequencies_from_partial_hessian(dof_labels, hessian_matrix, atoms):
    """Recover vibrational frequencies from a partial Hessian matrix.

    Parameters
    ----------
    dof_labels : list[str]
        Degree-of-freedom labels (for example ``1X``).
    hessian_matrix : list[list[float]]
        Hessian values in eV/Å² units for selected DOFs.
    atoms : ase.Atoms
        Atomic structure used to assign atomic masses.

    Returns
    -------
    list[float]
        Frequencies in cm⁻¹ sorted in descending order.
    """
    dof_masses = []
    for label in dof_labels:
        match = _DOF_LABEL_RE.match(label)
        if not match:
            raise ValueError(f"Invalid DOF label: {label}")
        atom_index = int(match.group(1)) - 1
        dof_masses.append(atoms[atom_index].mass)

    hessian = np.array(hessian_matrix, dtype=float)
    mass = np.sqrt(np.outer(dof_masses, dof_masses))
    dynamical = -hessian / mass
    eigenvals = np.linalg.eigvalsh(dynamical)

    # Conversion: sqrt(eV/amu)/Ang -> cm-1
    factor_cm = 521.4708983725064
    frequencies = factor_cm * np.sqrt(np.maximum(eigenvals, 0.0))

    return sorted((float(v) for v in frequencies), reverse=True)

# ============================================================
# Structure utilities
# ============================================================

def formula_from_atom_order(atoms):
    """Build a chemical formula preserving first-seen element order."""
    counts = {}
    order = []

    for atom in atoms:
        sym = atom.symbol
        if sym not in counts:
            counts[sym] = 0
            order.append(sym)
        counts[sym] += 1

    formula = ""
    for sym in order:
        n = counts[sym]
        formula += f"{sym}{n if n > 1 else ''}"

    return formula


def lattice_vectors(atoms, precision=8):
    """Return lattice vectors rounded for stable YAML serialization."""
    fmt = f"{{:.{precision}f}}"
    return [
        InlineList([float(fmt.format(x)) for x in vec])
        for vec in atoms.cell.array
    ]

def geometry_direct_strings(atoms, precision=8):
    """
    Return fractional coordinates as list of strings.

    Using strings avoids YAML ambiguity and keeps diffs clean.
    """
    scaled = atoms.get_scaled_positions()

    fmt = f"{{:.{precision}f}}"
    lines = []

    for atom, pos in zip(atoms, scaled):
        x, y, z = pos
        lines.append(
            f"{atom.symbol} {fmt.format(x)} {fmt.format(y)} {fmt.format(z)}"
        )

    return lines


# ============================================================
# Main parser
# ============================================================

def parse_vasp_frequency(outcar_path, average_pairs=False):
    """Parse a VASP frequency OUTCAR file into pymkmkit YAML schema.

    Parameters
    ----------
    outcar_path : str | pathlib.Path
        Path to the input OUTCAR file.
    average_pairs : bool, default=False
        If ``True``, average vibrational modes in sequential pairs.

    Returns
    -------
    dict
        Structured dictionary ready to be serialized to YAML.
    """
    path = Path(outcar_path)

    if not path.exists():
        raise FileNotFoundError(outcar_path)

    text = path.read_text(errors="ignore")

    # ASE reads final structure & energy
    atoms = read(outcar_path)

    incar = extract_incar_settings(text)
    potcar = extract_potcar_info(text)
    real_freqs, imag_freqs = extract_frequencies(text)
    hessian_dofs, hessian_matrix = extract_perturbed_hessian(text)
    ionic_energies = extract_ionic_energies(text)
    if ionic_energies:
        electronic_energy = ionic_energies[0]
    else:
        electronic_energy = float(atoms.get_potential_energy())

    # ---- vibrational processing ----

    pairing_note = None

    if average_pairs:
        real_freqs, imag_freqs, pairing_note = average_mode_pairs(real_freqs, imag_freqs)

    vibration_block = {
        "frequencies_cm-1": real_freqs,
        "imaginary_cm-1": imag_freqs if imag_freqs else None,
    }

    if hessian_dofs and hessian_matrix:
        vibration_block["partial_hessian"] = {
            "dof_labels": hessian_dofs,
            "matrix": hessian_matrix,
        }

    if average_pairs:
        vibration_block["paired_modes_averaged"] = True

        if pairing_note is not None:
            vibration_block["pairing_note"] = pairing_note

    # ---- assemble output ----

    return {
        "pymkmkit": {
            "version": get_version(),
            "generated": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        },
        "structure": {
            "formula": formula_from_atom_order(atoms),
            "n_atoms": len(atoms),
            "lattice_vectors": lattice_vectors(atoms),
            "coordinates_direct": geometry_direct_strings(atoms),
            "pbc": [bool(x) for x in atoms.pbc],
        },
        "calculation": {
            "code": "VASP",
            "type": "frequency",
            "incar": incar,
            "potcar": potcar,
        },
        "energy": {
            "electronic": electronic_energy,
        },
        "vibrations": vibration_block,
    }


def parse_vasp_optimization(outcar_path):
    """Parse a VASP optimization OUTCAR file into pymkmkit YAML schema."""
    path = Path(outcar_path)

    if not path.exists():
        raise FileNotFoundError(outcar_path)

    text = path.read_text(errors="ignore")

    # geometry optimizations should store the last ionic step
    atoms = _read_last_optimization_atoms(outcar_path, text)

    incar = extract_incar_settings(text)
    potcar = extract_potcar_info(text)

    ionic_energies = extract_ionic_energies(text)
    total_energies = extract_total_energies(text)

    if ionic_energies:
        electronic_energy = ionic_energies[-1]
    elif total_energies:
        electronic_energy = total_energies[-1]
    else:
        electronic_energy = float(atoms.get_potential_energy())

    return {
        "pymkmkit": {
            "version": get_version(),
            "generated": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        },
        "structure": {
            "formula": formula_from_atom_order(atoms),
            "n_atoms": len(atoms),
            "lattice_vectors": lattice_vectors(atoms),
            "coordinates_direct": geometry_direct_strings(atoms),
            "pbc": [bool(x) for x in atoms.pbc],
        },
        "calculation": {
            "code": "VASP",
            "type": "optimization",
            "incar": incar,
            "potcar": potcar,
        },
        "energy": {
            "electronic": electronic_energy,
        },
    }

def average_mode_pairs(real_freqs, imag_freqs):
    """
    Average real frequencies in sequential pairs (1+2, 3+4, ...).

    TS edge case:
    If number of real modes is odd and exactly one imaginary mode exists,
    drop the last unpaired real mode and proceed.
    """

    note = None

    if len(real_freqs) % 2 == 1:
        if len(imag_freqs) == 1 and len(real_freqs) > 0:
            dropped = real_freqs[-1]
            real_freqs = real_freqs[:-1]

            note = (
                "Odd number of real modes with one imaginary mode (TS case). "
                f"Dropped unpaired real mode ({dropped:.6f} cm-1) "
                "before pair averaging."
            )
        else:
            raise ValueError(
                "Cannot pair-average: odd number of real modes and "
                "imaginary mode count != 1."
            )

    averaged_real = [
        (real_freqs[i] + real_freqs[i + 1]) / 2.0
        for i in range(0, len(real_freqs), 2)
    ]

    if len(imag_freqs) % 2 == 1 and len(imag_freqs) != 1:
        raise ValueError(
            "Cannot pair-average: odd number of imaginary modes "
            "that cannot be paired sequentially."
        )

    if len(imag_freqs) == 1:
        averaged_imag = imag_freqs
    else:
        averaged_imag = [
            (imag_freqs[i] + imag_freqs[i + 1]) / 2.0
            for i in range(0, len(imag_freqs), 2)
        ]

    return averaged_real, averaged_imag, note
