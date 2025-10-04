from flask import Blueprint, render_template, request
from services.services import *

bp = Blueprint("routes", __name__, url_prefix="")

@bp.route("/")
def main_page():
    """ Base route """
    return render_template("index.html")


@bp.route("/search", methods=["GET", "POST"])
def search_asteroid_daesignation():
    """ Search APIs for info + parameters? POST return text data """
    if request.method == "POST":
        pass
    else:
        return render_template("search.html")
    
    
@bp.route("/result")
def result_page(threat=False):
    """ Provide current most-like trajectory for PHA """
    if threat:
        return render_template("result.html")
    else:
        render_template("no_threat.html")


@bp.route("/simulation")
def simulation():
    """ Run simulation of hitting Earth """
    return render_template("simulation.html")