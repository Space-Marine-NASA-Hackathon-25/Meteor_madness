import requests


FIELDS = [
    "full_name", "H", "G", "density", "albedo",
    "e", "a", "q", "i", "om", "w", "ma", "tp", "per", "n", "ad", "moid_ld"]


def call_sentry(designation):
    pass


def call_horizons():
    pass


def call_smdb(designation, fields=FIELDS):
    """ Get data from JPL SBDB for PHA """
    request_dict = {
        "fields": ",".join(fields),
        # "AND": ["spkid|EQ|54543542"], # for query a specific object (not working); for close approach 2025 QB21
        "sb-group":"pha", 
        "limit":"15" # testing purposes only!
    }
    response = requests.get('https://ssd-api.jpl.nasa.gov/sbdb_query.api', params=request_dict)
    return response.json()


def convert_horizons2spice():
    pass


def convert_smdb2spice():
    pass

