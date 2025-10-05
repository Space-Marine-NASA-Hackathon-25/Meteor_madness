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
      "des": "(29075) 1950 DA",
      "diameter": 1300,
      "mass": 1.2e12,
      "GM": 1.2e-3,
      "spec": "S-type"
    },
    {
      "des": "(385343) 2002 LV",
      "diameter": 1420,
      "mass": 1.0e12,
      "GM": 1.0e-3,
      "spec": "S-type"
    },
    {
      "des": "(456938) 2007 YV56",
      "diameter": 1000,
      "mass": 5.0e11,
      "GM": 5.0e-4,
      "spec": "S-type"
    },
    {
      "des": "(454101) 2013 BP73",
      "diameter": 660,
      "mass": 2.0e10,
      "GM": 2.0e-5,
      "spec": "S-type"
    },
    {
      "des": "(2012) UE34",
      "diameter": 130,
      "mass": 1.0e6,
      "GM": 1.0e-9,
      "spec": "S-type"
    }
  ]
}

horizons_data = {
    "vectors": [
    {
      "datetime": "2025-10-05 00:00",
      "x": 1.234, "y": 0.987, "z": 0.123,
      "vx": -0.012, "vy": 0.023, "vz": 0.001
    },
    {
      "datetime": "2025-10-05 00:00",
      "x": 2.543, "y": -0.876, "z": 0.321,
      "vx": 0.015, "vy": -0.019, "vz": 0.002
    },
    {
      "datetime": "2025-10-05 00:00",
      "x": -1.876, "y": 1.234, "z": -0.145,
      "vx": -0.020, "vy": 0.014, "vz": -0.001
    },
    {
      "datetime": "2025-10-05 00:00",
      "x": 0.543, "y": -1.098, "z": 0.234,
      "vx": 0.010, "vy": 0.012, "vz": 0.003
    },
    {
      "datetime": "2025-10-05 00:00",
      "x": -0.321, "y": 0.654, "z": -0.098,
      "vx": -0.005, "vy": -0.008, "vz": 0.001
    }
  ]
}


nameMeteor = "(456938) 2007 YV56"
indexMeteor = 0

for x in sbdb_data["data"]:
    if nameMeteor == x["des"]:
        break
    
    indexMeteor += 1


if indexMeteor == 5:
    print("It's good meteor")
    
for k in get_paths_to_kernels():
    spice.furnsh(k)


phys = sbdb_data["data"][indexMeteor]
vec = horizons_data["vectors"][indexMeteor]

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
kinetic_energy = 0.5 * mass * (v_mag * 1000) ** 2   # ← додав **


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


def draw_crater(map_obj, lat, lon, kinetic_energy):
    crater_radius = (kinetic_energy / (4.3e14)) ** (1/4)  # ← замінив пробіли на **
    crater_radius_km = crater_radius / 1000
    folium.Circle(
        location=[lat, lon],
        radius=crater_radius,
        color="red",
        fill=True,
        fill_opacity=0.4,
        popup=f"Crater radius: {crater_radius_km:.2f} km"
    ).add_to(map_obj)
    return map_obj


m = folium.Map(location=[lat, lon], zoom_start=5)
folium.Marker([lat, lon]).add_to(m)
m = draw_crater(m, lat, lon, kinetic_energy)
m.save("templates/meteor_entry_map.html")
