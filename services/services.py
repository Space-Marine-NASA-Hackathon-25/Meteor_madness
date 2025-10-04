from pathlib import Path
import spiceypy as spice
from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Sun
from poliastro.twobody import Orbit
import numpy as np
from numpy import linalg, arccos, degrees


DIR_KERNELS = "kernels"
KERNELS_PATHS = ["spk/de432s.bsp", "lsk/naif0012.tls", "pck/gm_de440.tpc"]

def get_paths_to_kernels(dir_kernels=DIR_KERNELS, paths=KERNELS_PATHS):
    BASE = Path(__file__).resolve().parent.parent
    KERNELS = BASE / dir_kernels
    return [str(KERNELS / p) for p in paths]




sbdb_data = {
    "data": [
        {
            "full_name": "(164207) Cardea",
            "diameter": 200,
            "mass": 5e9,
            "GM": 0,
            "spec": "C-type"
        }
    ]
}

horizons_data = {
    "vectors": [
        {
            "datetime": "2025-10-04 00:00",
            "x": 1.234,
            "y": 0.987,
            "z": 0.123,
            "vx": -0.012,
            "vy": 0.023,
            "vz": 0.001
        }
    ]
}

phys = sbdb_data.get("data", [{}])[0]
vec = horizons_data.get("vectors", [{}])[0]

AU_in_km = 149597870.7
DAY_in_s = 86400

x = vec.get('x', 0) * AU_in_km
y = vec.get('y', 0) * AU_in_km
z = vec.get('z', 0) * AU_in_km
vx = vec.get('vx', 0) * AU_in_km / DAY_in_s
vy = vec.get('vy', 0) * AU_in_km / DAY_in_s
vz = vec.get('vz', 0) * AU_in_km / DAY_in_s
datetime_str = vec.get('datetime', '2000-01-01 00:00')

for k in get_paths_to_kernels():
    spice.furnsh(k)

et = spice.utc2et(datetime_str)


state_vec, lt = spice.spkezr("399", et, "J2000", "LT+S", "0")  

r_rel = np.array([x, y, z]) - np.array(state_vec[:3])
v_rel = np.array([vx, vy, vz]) - np.array(state_vec[3:])

r_vec = r_rel * u.km
v_vec = v_rel * u.km / u.s
epoch = Time(datetime_str)

orbit = Orbit.from_vectors(Sun, r_vec, v_vec, epoch)

r_hat = r_vec.value / linalg.norm(r_vec.value)
v_hat = v_vec.value / linalg.norm(v_vec.value)
entry_angle = degrees(arccos(-r_hat @ v_hat))

mass = phys.get("mass", 1e10)
v_mag = linalg.norm(v_vec.value)
kinetic_energy = 0.5 * mass * (v_mag * 1000)**2

print(f"Entry coordinates (km): {r_vec}")
print(f"Entry speed (km/s): {v_mag:.3f}")
print(f"Entry angle (°): {entry_angle:.2f}")
print(f"Meteor mass (kg): {mass}")
print(f"Entry kinetic energy (J): {kinetic_energy:.3e}")
print("\n=== Poliastro Orbit Object ===")
print(orbit)
