delta_mod = {"14.01": "me", "28.03": "me2", "42.04": "me3", "42.01": "ac", "79.96": "p"}

from pyteomics import mass

ion_formulas = {
    "a": {"a": mass.Composition(formula="H-2O-1" + "C-1O-1")},
    "b": {"b": mass.Composition(formula="H-2O-1")},
    "x": {"x": mass.Composition(formula="H-2O-1" + "CO2")},
    "y": {"y": mass.Composition(formula="")},
    "cdot": {"cdot": mass.Composition(formula="H-2O-1" + "NH3")},
    "c": {"c": mass.Composition(formula="H-2O-1" + "NH4")},
    "c-1": {"c-1": mass.Composition(formula="H-2O-1" + "NH2")},
    "c+1": {"c+1": mass.Composition(formula="H-2O-1" + "NH5")},
    "zdot": {"zdot": mass.Composition(formula="H-2O-1" + "N-1" + "OH")},
    "z+1": {"z+1": mass.Composition(formula="H-2O-1" + "N-1" + "OH2")},  # z+M(H)
    "z+2": {"z+2": mass.Composition(formula="H-2O-1" + "N-1" + "OH3")},  # z+M(2H)
    "z+3": {"z+3": mass.Composition(formula="H-2O-1" + "N-1" + "OH4")},
    "c-zdot": {
        "c-zdot": mass.Composition(formula="H-2O-1" + "H-2O-1" + "OH5")
    },  # -O   = c + z (from msnbase issue 82)
    "c-z+1": {"c-z+1": mass.Composition(formula="H-2O-1" + "H-2O-1" + "OH6")},
    "cdot-zdot": {"cdot-zdot": mass.Composition(formula="H-2O-1" + "H-2O-1" + "OH4")},
    "cdot-z+1": {"cdot-z+1": mass.Composition(formula="H-2O-1" + "H-2O-1" + "OH5")},
    "n-n": {"n-n": mass.Composition(formula="P-1")},
    "b-y": {"b-y": mass.Composition(formula="H-2O-1" + "H-2O-1" + "")},
    "a-x": {"a-x": mass.Composition(formula="H-2O-1" + "H-2O-1" + "C-1O-1" + "CO2")},
}

ion_direction = {
    "a": "n-term",
    "b": "n-term",
    "x": "c-term",
    "y": "c-term",
    "cdot": "n-term",
    "c": "n-term",
    "c-1": "n-term",
    "c+1": "n-term",
    "zdot": "c-term",
    "z+1": "c-term",
    "z+2": "c-term",
    "z+3": "c-term",
    "c-zdot": "intern",
    "c-z+1": "intern",
    "cdot-zdot": "intern",
    "cdot-z+1": "intern",
    "n-n": "intern",
    "b-y": "intern",
    "a-x": "intern",
}


# ---------------------------------------------------------------------------- #
#                               For visualization                              #
# ---------------------------------------------------------------------------- #


colors = [
    "#9b2226",
    "#005f73",
    "#ee9b00",
    "#0a9396",
    "#94d2bd",
    "#ca6702",
    "#e9d8a6",
    "#bb3e03",
    "#001219",
    "#006BA6",
    "#35A7FF",
    "#EFA8B8",
    "#BFACC8",
    "#476A6F",
    "#7067CF",
    "#364156",
    "#98FB98",
    "#8A2BE2",
    "#35682D",
    "#252850",
    "#7E7B52",
]

colors = colors + (["#808080"] * 1000)
