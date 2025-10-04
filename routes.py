from flask import Blueprint, render_template
from services import *

bp = Blueprint("routes", __name__, url_prefix="")

@bp.route("/")
def main_page():
    return render_template("index.html")