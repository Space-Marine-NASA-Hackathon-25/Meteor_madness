from pathlib import Path
from services.api_calls import *

DIR_KERNELS = "kernels"
KERNELS_PATHS = ["spk/de432s.bsp", "lsk/naif0012.tls", "pck/gm_de440.tpc"]

def get_paths_to_kernels(dir_kernels=DIR_KERNELS, paths=KERNELS_PATHS):
    """ Returns absolute path for each kernel file """
    # https://naif.jpl.nasa.gov/pub/naif/generic_kernels/
    BASE = Path(__file__).resolve().parent
    KERNELS = BASE / dir_kernels
    paths_to_kernels = []
    for p in paths:
        paths_to_kernels.append(str(KERNELS / p))
    return paths_to_kernels