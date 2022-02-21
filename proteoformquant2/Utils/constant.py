delta_mod = { 
        "14.01": "me",
        "28.03": "me2",
        "42.04": "me3",
        "42.01": "ac",
        "79.96": "p"
    }

from pyteomics import mass
ion_formulas = {
        'a':        {'a': mass.Composition(formula='H-2O-1' + 'C-1O-1')},
        'b':        {'b': mass.Composition(formula='H-2O-1')},
        'x':        {'x': mass.Composition(formula='H-2O-1' + 'CO2')},
        'y':        {'y': mass.Composition(formula='')},

        'cdot':     {'cdot': mass.Composition(formula='H-2O-1' +'NH3')},
        'c':        {'c': mass.Composition(formula='H-2O-1' +'NH4')},
        'c-1':      {'c-1': mass.Composition(formula='H-2O-1' + 'NH2')},
        'c+1':      {'c+1': mass.Composition(formula='H-2O-1' + 'NH5')},

        'zdot':     {'zdot': mass.Composition(formula='H-2O-1' + 'N-1' + 'OH')},
        'z+1':      {'z+1': mass.Composition(formula='H-2O-1' + 'N-1' + 'OH2')}, #z+M(H)
        'z+2':      {'z+2': mass.Composition(formula='H-2O-1' + 'N-1' + 'OH3')}, #z+M(2H)
        'z+3':      {'z+3': mass.Composition(formula='H-2O-1' + 'N-1' + 'OH4')},

        "c-zdot":      {'cIzdot': mass.Composition(formula='H-2O-1' +'H-2O-1'+'OH5')},  # -O   = c + z (from msnbase issue 82)
        "c-z+1":      {'cIz+1': mass.Composition(formula='H-2O-1' +'H-2O-1'+'OH6')},
        "cdot-zdot":      {'cdotIzdot': mass.Composition(formula='H-2O-1' +'H-2O-1' + 'OH4')},
        "cdot-z+1":      {'cdotIz+1': mass.Composition(formula='H-2O-1' +'H-2O-1' +'OH5')},
        "n-n":     {'nIn': mass.Composition(formula='P-1')},
        "b-y":      {'bIy': mass.Composition(formula='H-2O-1' + 'H-2O-1' + '')},
        "a-x":      {'aIx': mass.Composition(formula='H-2O-1' + 'H-2O-1' + 'C-1O-1'+ 'CO2')}
    }


# ---------------------------------------------------------------------------- #
#                               For visualization                              #
# ---------------------------------------------------------------------------- #


colors = [
    
    "#9b2226",
    "#005f73",
    "#bb3e03",
    "#0a9396",
    "#94d2bd",
    "#ca6702",
    "#e9d8a6",
    "#ee9b00",
    "#001219"
    
]