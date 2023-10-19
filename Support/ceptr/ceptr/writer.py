"""Helper functions for printing."""

def comment(string):
    """Comment a string."""
    return f"""// {string}"""


def writer(fstream, string=""):
    """Write string to file followed by newline."""
    fstream.write(string)
    fstream.write("\n")
