from pathlib import Path
from ase.io import read
from pymkmkit.yaml_writer import InlineList
from pymkmkit._version import get_version
from datetime import datetime

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

    # ---- vibrational processing ----

    pairing_note = None

    if average_pairs:
        real_freqs, pairing_note = average_mode_pairs(real_freqs, imag_freqs)

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
            "generated": datetime.utcnow().isoformat() + "Z",
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
            "electronic": float(atoms.get_potential_energy()),
        },
        "vibrations": vibration_block,
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

    averaged = [
        (real_freqs[i] + real_freqs[i + 1]) / 2.0
        for i in range(0, len(real_freqs), 2)
    ]

    return averaged, note
