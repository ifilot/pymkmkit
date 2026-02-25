from importlib.metadata import version, PackageNotFoundError


def get_version():
    """Return the installed package version.

    Returns
    -------
    str
        Installed ``pymkmkit`` version string, or ``"unknown"`` when the
        distribution metadata cannot be resolved (for example in editable or
        source-only environments).
    """
    try:
        return version("pymkmkit")
    except PackageNotFoundError:
        return "unknown"
