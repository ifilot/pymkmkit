from pathlib import Path
import re
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
    energies = []

    for line in text.splitlines():
        if "energy  without entropy=" in line and "energy(sigma->0)" in line:
            try:
                energies.append(float(line.split("=")[-1].strip()))
            except ValueError:
                continue

    return energies


def extract_total_energies(text):
    energies = []

    for line in text.splitlines():
        if "free  energy   TOTEN" in line:
            try:
                energies.append(float(line.split("=")[-1].split()[0]))
            except ValueError:
                continue

    return energies


def _parse_atomic_symbols(text):
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
    for i in range(len(lines) - 1, -1, -1):
        if "direct lattice vectors" in lines[i]:
            vectors = []
            for j in range(i + 1, i + 4):
                parts = lines[j].split()
                vectors.append([float(parts[0]), float(parts[1]), float(parts[2])])
            return vectors

    raise ValueError("Could not parse lattice vectors from OUTCAR")


def _parse_last_direct_positions(lines, n_atoms):
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
    lines = text.splitlines()
    symbols = _parse_atomic_symbols(text)
    cell = _parse_last_lattice_vectors(lines)
    scaled_positions = _parse_last_direct_positions(lines, len(symbols))

    atoms = Atoms(symbols=symbols, cell=cell, pbc=True)
    atoms.set_scaled_positions(scaled_positions)
    return atoms


def _read_last_optimization_atoms(outcar_path, text):
    try:
        return read(outcar_path, index=-1)
    except (ParseError, UnknownFileTypeError, KeyError, ValueError):
        return _parse_atoms_from_outcar_text(text)

def extract_frequencies(text):
    """
    Extract vibrational frequencies from OUTCAR.

    VASP prints the frequency list twice.
    We detect the repetition by checking when the
    mode index restarts at 1 and stop reading.
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

            try:
                value = float(parts[-2])  # cm-1 column
            except Exception:
                continue

            if "f/i=" in line:
                imaginary.append(-abs(value))
            else:
                real.append(value)

    return real, imaginary

# ============================================================
# Structure utilities
# ============================================================

def formula_from_atom_order(atoms):
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
    path = Path(outcar_path)

    if not path.exists():
        raise FileNotFoundError(outcar_path)

    text = path.read_text(errors="ignore")

    # ASE reads final structure & energy
    atoms = read(outcar_path)

    incar = extract_incar_settings(text)
    potcar = extract_potcar_info(text)
    real_freqs, imag_freqs = extract_frequencies(text)
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
