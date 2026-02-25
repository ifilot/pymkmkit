import yaml


class InlineList(list):
    """Force YAML to print list in flow style."""
    pass


def represent_inline_list(dumper, data):
    """Serialize :class:`InlineList` values in YAML flow style.

    Parameters
    ----------
    dumper : yaml.Dumper
        Active PyYAML dumper instance.
    data : InlineList
        Sequence to render.

    Returns
    -------
    yaml.nodes.SequenceNode
        YAML node representation of the list.
    """
    return dumper.represent_sequence(
        "tag:yaml.org,2002:seq",
        data,
        flow_style=True
    )


yaml.SafeDumper.add_representer(InlineList, represent_inline_list)


def clean_none(d):
    """Recursively remove keys/items whose value is ``None``.

    Parameters
    ----------
    d : Any
        Arbitrary Python data structure destined for YAML output.

    Returns
    -------
    Any
        Cleaned structure with ``None`` entries removed while preserving
        :class:`InlineList` instances.
    """
    # Preserve dicts
    if isinstance(d, dict):
        return {k: clean_none(v) for k, v in d.items() if v is not None}

    # Preserve InlineList type (CRUCIAL)
    if isinstance(d, InlineList):
        return InlineList(clean_none(v) for v in d)

    # Normal lists
    if isinstance(d, list):
        return [clean_none(v) for v in d]

    return d


def write_yaml(data, filename):
    """Write parser output to YAML while preserving intended formatting.

    Parameters
    ----------
    data : dict
        Data structure to serialize.
    filename : str | pathlib.Path
        Target YAML file path.
    """
    with open(filename, "w") as f:
        yaml.dump(
            clean_none(data),
            f,
            sort_keys=False,
            Dumper=yaml.SafeDumper,
        )
