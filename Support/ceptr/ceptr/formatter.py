"""Helper functions for CPP formatting."""


def format_species(species):
    """Remove characters not allowed in preprocessor defines."""
    s = species.strip()
    # Ionic species
    if s[-1] == "-":
        s = s[:-1] + "n"
    if s[-1] == "+":
        s = s[:-1] + "p"
    # Excited species
    s = s.replace("*", "D")
    # Parenthesis
    s = s.replace("(", "_").replace(")", "")
    # Dash
    s = s.replace("-", "")
    return s
