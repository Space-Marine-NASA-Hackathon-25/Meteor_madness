from pathlib import Path
import spiceypy as spice
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import GCRS, ITRS, CartesianRepresentation, EarthLocation
from poliastro.bodies import Sun
from poliastro.twobody import Orbit
from numpy import linalg, arccos, degrees
from service_helpers.helpers import get_paths_to_kernels, au_in_km, get_datetimes
from service_helpers.api_calls import call_horizons, call_sentry, call_smdb
import folium


AU_in_km = au_in_km()
DAY_in_s = 86400
DES = "2008 JL3"
EPOCH = None


def get_phys_meteor(des=DES):
    threat = True
    try:
        mass, diameter, energy = call_sentry(des)
    except ValueError:
        threat = False
        mass = diameter = energy = None
    return threat, mass, diameter, energy


def get_earth_position():
    for k in get_paths_to_kernels():
        spice.furnsh(k)
    _, et, _ = get_datetimes()
    state_earth, _ = spice.spkezr("399", et, "J2000", "NONE", "0")
    earth_pos = np.array(state_earth[:3])
    earth_vel = np.array(state_earth[3:])
    spice.kclear()
    # print("from get_earth_position", earth_pos, earth_vel)
    return earth_pos, earth_vel
    

def get_relative_vecs(des=DES):
    global EPOCH
    x, y, z, vx, vy, vz, EPOCH = call_horizons(des)
    meteor_pos = [x, y, z]
    meteor_vel = [vx, vy, vz]
    earth_pos, earth_vel = get_earth_position()
    r_rel = meteor_pos - earth_pos
    v_rel = meteor_vel - earth_vel
    r_vec = r_rel * u.km
    v_vec = v_rel * u.km / u.s
    # print("from get_relative_vecs", x, y, z, vx, vy, vz, EPOCH, r_vec, v_vec)
    return r_vec, v_vec


def sim_orbit_n_conseq(des=DES):
    print("from last")
    r_vec, v_vec = get_relative_vecs()
    orbit = Orbit.from_vectors(Sun, r_vec, v_vec, EPOCH)
    r_hat = r_vec.value / linalg.norm(r_vec.value)
    v_hat = v_vec.value / linalg.norm(v_vec.value)
    entry_angle = degrees(arccos(-r_hat @ v_hat))
    
    _, mass, _, _ = get_phys_meteor(des)
    v_mag = linalg.norm(v_vec.value)
    kinetic_energy = 0.5 * mass * (v_mag * 1000)**2
    gcrs = GCRS(CartesianRepresentation(r_vec), obstime=EPOCH)
    itrs = gcrs.transform_to(ITRS(obstime=EPOCH))
    earth_location = EarthLocation(itrs.x, itrs.y, itrs.z)
    lat = earth_location.lat.deg
    lon = earth_location.lon.deg
    alt = earth_location.height.to(u.km).value
    
    # map
    m = get_map_obj(lat, lon)
    draw_crater(m, lat, lon, kinetic_energy)
    m.save("templates/meteor_entry_map.html")
    
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
    
    return orbit, r_vec, mass, entry_angle, v_mag, kinetic_energy, lat, lon


def get_map_obj(lat, lon):
    m = folium.Map(location=[lat, lon], zoom_start=5)
    folium.Marker([lat, lon]).add_to(m)
    return m
   

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

# sim_orbit_n_conseq()
