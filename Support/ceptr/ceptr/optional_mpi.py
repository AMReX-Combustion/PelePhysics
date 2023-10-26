"""Optional import of mpi4py, only if available in the environment."""

import sys

try:
    from mpi4py import MPI

except ModuleNotFoundError:
    pass


def use_mpi():
    """Detect if mpi4py was imported."""
    return "mpi4py" in sys.modules


def world_comm():
    """Return the world communicator."""
    return MPI.COMM_WORLD
