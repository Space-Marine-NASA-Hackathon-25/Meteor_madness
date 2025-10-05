from pathlib import Path
import spiceypy as spice
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import GCRS, ITRS, CartesianRepresentation, EarthLocation
from poliastro.bodies import Sun
from poliastro.twobody import Orbit
from numpy import linalg, arccos, degrees
from service_helpers.helpers import get_paths_to_kernels, au_in_km
import folium


AU_in_km = au_in_km()
DAY_in_s = 86400



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
            "datetime": "2025-10-05 00:00",
            "x": 1.234,
            "y": 0.987,
            "z": 0.123,
            "vx": -0.012,
            "vy": 0.023,
            "vz": 0.001
        }
    ]
}


for k in get_paths_to_kernels():
    spice.furnsh(k)


phys = sbdb_data["data"][0]
vec = horizons_data["vectors"][0]

x = vec["x"] * AU_in_km
y = vec["y"] * AU_in_km
z = vec["z"] * AU_in_km
vx = vec["vx"] * AU_in_km / DAY_in_s
vy = vec["vy"] * AU_in_km / DAY_in_s
vz = vec["vz"] * AU_in_km / DAY_in_s
datetime_str = vec["datetime"]

epoch = Time(datetime_str)
et = spice.utc2et(datetime_str)


state_earth, _ = spice.spkezr("399", et, "J2000", "NONE", "0")
earth_pos = np.array(state_earth[:3])
earth_vel = np.array(state_earth[3:])

meteor_pos = np.array([x, y, z])
meteor_vel = np.array([vx, vy, vz])

r_rel = meteor_pos - earth_pos
v_rel = meteor_vel - earth_vel

r_vec = r_rel * u.km
v_vec = v_rel * u.km / u.s


orbit = Orbit.from_vectors(Sun, r_vec, v_vec, epoch)


r_hat = r_vec.value / linalg.norm(r_vec.value)
v_hat = v_vec.value / linalg.norm(v_vec.value)
entry_angle = degrees(arccos(-r_hat @ v_hat))

mass = phys.get("mass", 1e10)
v_mag = linalg.norm(v_vec.value)
kinetic_energy = 0.5 * mass * (v_mag * 1000)**2


gcrs = GCRS(CartesianRepresentation(r_vec), obstime=epoch)
itrs = gcrs.transform_to(ITRS(obstime=epoch))
earth_location = EarthLocation(itrs.x, itrs.y, itrs.z)

lat = earth_location.lat.deg
lon = earth_location.lon.deg
alt = earth_location.height.to(u.km).value


spice.kclear()


print(f"Entry coordinates (km): {r_vec}")
print(f"Entry speed (km/s): {v_mag:.3f}")
print(f"Entry angle (°): {entry_angle:.2f}")
print(f"Meteor mass (kg): {mass}")
print(f"Entry kinetic energy (J): {kinetic_energy:.3e}")

print(f"{lat:.4f}°")
print(f"{lon:.4f}°")
print(f"{alt:.2f} km")

print("\n=== Poliastro Orbit Object ===")
print(orbit)

m = folium.Map(location=[lat, lon], zoom_start=5)
folium.Marker([lat, lon]).add_to(m)
m.save("meteor_entry_map.html")
