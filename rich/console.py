import re


class Console:
    """Minimal Console-compatible wrapper used for local output rendering."""

    def print(self, message: str) -> None:
        plain = re.sub(r"\[/?[a-zA-Z0-9_#=-]+\]", "", str(message))
        print(plain)
