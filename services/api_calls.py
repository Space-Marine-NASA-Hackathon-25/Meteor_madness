import requests
import json
from astroquery.jplhorizons import Horizons
from helpers import get_datetimes, get_paths_to_kernels, au_in_km
import spiceypy as spice
import numpy as np

SUN_ID   = 10
EARTH_ID = 399
AU_in_km = au_in_km()

FIELDS = [
    "full_name", "H", "G", "density", "albedo",
    "e", "a", "q", "i", "om", "w", "ma", "tp", "per", "n", "ad", "moid_ld"]


def call_sentry(designation):
    pass


def call_horizons(designation):
    """ Get data from JPL Horizons for the asteroid """
    spice.furnsh(get_paths_to_kernels())
    _, _, jd = get_datetimes()
    obj = Horizons(
        id=f"{designation}",
        id_type="smallbody",
        location="@sun",
        epochs=jd
    )
    orbital_el = get_orbital_elements(obj)
    vector_el = get_vector_elements(obj)
    spice.kclear()
    return orbital_el, vector_el


def get_orbital_elements(obj):
    els_dict = dict(obj.elements()[0]) # obj.elements() return astropy.table by default
    # print(els_dict)
    _, gm_sun_pre = spice.bodvcd(bodyid=SUN_ID, item="GM", maxn=1)
    gm_sun = gm_sun_pre[0]
    orbital_elements = [
         spice.convrt(els_dict["q"],"AU", "km"),
         els_dict["e"],
         np.radians(els_dict["incl"]),
         np.radians(els_dict["Omega"]), 
         np.radians(els_dict["w"]),
         np.radians(els_dict["M"]),
         spice.utc2et(els_dict["datetime_str"]),
         gm_sun # for calculations in SPICE
     ]
    return orbital_elements


def get_vector_elements(obj):
    els = dict(obj.vectors()[0])
    physical_elements = [
        els["H"],
        els["G"],
        els["x"],
        els["y"],
        els["z"],
        els["vx"],
        els["vy"],
        els["vz"],
        spice.convrt(els["range"], "AU", "km"),
        (els["range_rate"] * AU_in_km), # AU/day -> km/day
    ]
    return physical_elements

print(call_horizons("2008 JL3"))


def call_smdb(designation):
    """ Get data from JPL SBDB for the asteroid """
    response = requests.get(
        "https://ssd-api.jpl.nasa.gov/sbdb.api", 
        params={"des": designation})
    return json.dumps(response.json(), indent = 4)

# print(call_smdb("2008 JL3"))
# {
#     "object": {
#         "des": "2008 JL3",
#         "kind": "au",
#         "orbit_id": "12",
#         "fullname": "(2008 JL3)",
#         "prefix": null,
#         "orbit_class": {
#             "code": "APO",
#             "name": "Apollo"
#         },
#         "neo": true,
#         "spkid": "3410213",
#         "pha": false
#     },
#     "signature": {
#         "version": "1.3",
#         "source": "NASA/JPL Small-Body Database (SBDB) API"
#     },
#     "orbit": {
#         "first_obs": "2008-05-02",
#         "n_del_obs_used": null,
#         "condition_code": "7",
#         "n_dop_obs_used": null,
#         "not_valid_before": null,
#         "source": "JPL",
#         "soln_date": "2021-04-15 01:41:18",
#         "producer": "Otto Matic",
#         "t_jup": "3.493",
#         "pe_used": "DE441",
#         "model_pars": [],
#         "n_obs_used": 35,
#         "data_arc": "7",
#         "orbit_id": "12",
#         "elements": [
#             {
#                 "sigma": "0.00062",
#                 "label": "e",
#                 "name": "e",
#                 "units": null,
#                 "title": "eccentricity",
#                 "value": "0.547"
#             },
#             {
#                 "sigma": "0.0025",
#                 "label": "a",
#                 "name": "a",
#                 "title": "semi-major axis",
#                 "value": "2.15",
#                 "units": "au"
#             },
#             {
#                 "title": "perihelion distance",
#                 "value": "0.975",
#                 "units": "au",
#                 "sigma": "0.0002",
#                 "name": "q",
#                 "label": "q"
#             },
#             {
#                 "sigma": "0.00023",
#                 "name": "i",
#                 "label": "i",
#                 "title": "inclination; angle with respect to x-y ecliptic plane",
#                 "value": "0.892",
#                 "units": "deg"
#             },
#             {
#                 "units": "deg",
#                 "title": "longitude of the ascending node",
#                 "value": "40.4",
#                 "sigma": "0.0092",
#                 "name": "om",
#                 "label": "node"
#             },
#             {
#                 "title": "argument of perihelion",
#                 "units": "deg",
#                 "value": "156",
#                 "label": "peri",
#                 "name": "w",
#                 "sigma": "0.0041"
#             },
#             {
#                 "name": "ma",
#                 "label": "M",
#                 "sigma": "3.6",
#                 "units": "deg",
#                 "title": "mean anomaly",
#                 "value": "201"
#             },
#             {
#                 "title": "time of perihelion passage",
#                 "value": "2461508.719",
#                 "units": "TDB",
#                 "sigma": "12",
#                 "name": "tp",
#                 "label": "tp"
#             },
#             {
#                 "label": "period",
#                 "name": "per",
#                 "sigma": "2",
#                 "value": "1150",
#                 "title": "sidereal orbital period",
#                 "units": "d"
#             },
#             {
#                 "units": "deg/d",
#                 "title": "mean motion",
#                 "value": "0.312",
#                 "sigma": "0.00054",
#                 "label": "n",
#                 "name": "n"
#             },
#             {
#                 "sigma": "0.0039",
#                 "name": "ad",
#                 "label": "Q",
#                 "title": "aphelion distance",
#                 "value": "3.33",
#                 "units": "au"
#             }
#         ],
#         "two_body": null,
#         "not_valid_after": null,
#         "equinox": "J2000",
#         "rms": "0.37",
#         "sb_used": "SB441-N16",
#         "epoch": "2461000.5",
#         "moid_jup": "1.63",
#         "last_obs": "2008-05-09",
#         "cov_epoch": "2454590.5",
#         "moid": "0.000249",
#         "comment": null
#     }
# }


def convert_smdb2spice():
    pass

