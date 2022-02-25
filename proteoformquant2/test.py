import dash
from dash import html
from dash import dcc
from matplotlib.pyplot import figure
import plotly.graph_objects as go
from dash.dependencies import Input, Output
import pickle
from matplotlib.pyplot import cm
from Classes.elution_profile import ElutionProfile
import numpy as np
import dash_bootstrap_components as dbc
import plotly.io as pio
import json
from Utils import misc
import plotly.express as px
import math
from dash.exceptions import PreventUpdate
from scipy import stats
from kneed import KneeLocator
from statistics import mean
import operator
from itertools import compress

with open('pfq_out_obj_test_1b.pkl', 'rb') as inp:
    exp = pickle.load(inp) 


#S = exp.spectra["scan=2779"]
S = exp.spectra["scan=2991"]
#S = exp.spectra["scan=2497"]

psms = S.psms[0:4]

print("psms list")
print([psm.get_modification_brno() for psm in psms])


mod_pos_list = [] #position (1 based) of modified residues across multiple psm (e.g K14acK24ac K9me3 -> [9,14,24])
for psm in psms:
    for mod in psm.proteoform.get_modification_dict():
        if mod["location"] not in mod_pos_list:
            mod_pos_list.append(mod["location"])
mod_pos_list = sorted(mod_pos_list)
print("modification positions")
print(mod_pos_list)

mod_pos_list_n = sorted([*mod_pos_list, 1, len(psm.proteoform.get_peptide_sequence())+1])
mod_pos_list_c = sorted([*mod_pos_list, 0, len(psm.proteoform.get_peptide_sequence())])
print("modification positions direction dependnet N")
print(mod_pos_list_n)
print("modification positions direction dependnet C")
print(mod_pos_list_c)

intervals_n_term = [[mod_pos_list_n[p],mod_pos_list_n[p+1]-1] for p in range(len(mod_pos_list_n)-1)]
intervals_c_term = [[mod_pos_list_c[p]+1,mod_pos_list_c[p+1]] for p in range(len(mod_pos_list_c)-1)]
print("positions intervals N")
print(intervals_n_term)
print("positions intervals C")
print(intervals_c_term)


mod_mass_matrix = np.zeros( (len(psms), len(mod_pos_list)) )
for row, psm in enumerate(psms):
    print(psm.get_modification_brno())
    for mod in psm.proteoform.get_modification_dict():
        print(mod)
        col = mod_pos_list.index(mod["location"])
        mod_mass_matrix[row,col] = mod["monoisotopicMassDelta"]
print("modification mass matrix")
print(mod_mass_matrix)


#cterm
additive_mass_c = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1]+1))
for col in range(1, additive_mass_c.shape[1]):
    additive_mass_c[:,col] = np.add( additive_mass_c[:,col-1], mod_mass_matrix[:,col-1])

#nterm
additive_mass_n = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1]+1))
for col in range(additive_mass_n.shape[1]-1, 0, -1):
    #print(additive_mass_n[:,col-1])
    #print(mod_mass_matrix[:,col-1])
    additive_mass_n[:,col-1] = np.add( additive_mass_n[:,col], mod_mass_matrix[:,col-1])

print("Additive modification masses matrix direction dependent")
print(additive_mass_c)
print(additive_mass_n)

unique_matrix_c =  np.zeros((additive_mass_c.shape[0], additive_mass_c.shape[1]))
for col in range(unique_matrix_c.shape[1]):
    grp_mass = []
    for row, m in enumerate(additive_mass_c[:,col]):
        if m not in grp_mass:
            grp_mass.append(m)
        grp_index = grp_mass.index(m)
        unique_matrix_c[row,col] = grp_index + 1

unique_matrix_n =  np.zeros((additive_mass_n.shape[0], additive_mass_n.shape[1]))
for col in range(unique_matrix_n.shape[1]):
    grp_mass = []
    for row, m in enumerate(additive_mass_n[:,col]):
        if m not in grp_mass:
            grp_mass.append(m)
        grp_index = grp_mass.index(m)
        unique_matrix_n[row,col] = grp_index + 1

print("grouping matrix n and c term")
print(intervals_c_term)
print(unique_matrix_c)
print(intervals_n_term)
print(unique_matrix_n)


def reduce_unique_matrix(unique_matrix, intervals):

    n_col = unique_matrix.shape[1]
    n_row = unique_matrix.shape[0]

    #first pass to remove interval without any unique ions acroos proteoforms
    col_rm = []
    for col in range(unique_matrix.shape[1]):
        if (unique_matrix[:,col] == np.ones(n_row)).all(): #if no unique ions in interval
            #print("no unique ions in interval{0}".format( intervals[col] ))
            col_rm.append(col)
    
    unique_matrix = np.delete(unique_matrix, col_rm, 1)
    intervals = np.delete(intervals, col_rm, 0)


    pair_found = True
    while pair_found:
        pair_found = False
        prev_col = np.zeros(n_row)
        for col in range(unique_matrix.shape[1]):

            if (unique_matrix[:,col] == prev_col).all():
                pair_found = True
                
                unique_matrix = np.delete(unique_matrix, col-1, 1)
                intervals[col,:] = [intervals[col-1,0],intervals[col,1]]
                intervals = np.delete(intervals, col-1, 0)
                break
            prev_col = unique_matrix[:,col]

    return unique_matrix, intervals.tolist()
        
    


red_unique_matrix_n, red_intervals_n_term = reduce_unique_matrix(unique_matrix_n, intervals_n_term)
red_unique_matrix_c, red_intervals_c_term = reduce_unique_matrix(unique_matrix_c, intervals_c_term)



def get_intensities_matrix(unique_matrix, intervals, psms, direction):
    intensity_matrix =  np.ones((unique_matrix.shape[0], unique_matrix.shape[1]))
    for col in range(unique_matrix.shape[1]):
        print("new col")
        for row in range(unique_matrix.shape[0]):
            

            fragments_at_range = psms[row].get_fragments_at_range(intervals[col], direction= direction)
            print("dir = {4} GRP = {1} / {2}, range used: {3} Number of Possible fragments: {0} frag found:".format(len(fragments_at_range), 
                                                                                                                    unique_matrix[row,col],
                                                                                                                    psms[row].get_modification_brno(),
                                                                                                                    intervals[col],
                                                                                                                    direction))
            intensity_matrix[row,col] = psms[row].spectrum.get_sum_intens_fragment_list(psms[row], fragments_at_range)

    return intensity_matrix

intensity_matrix_c = get_intensities_matrix(red_unique_matrix_c, red_intervals_c_term,  psms, "c-term") 
intensity_matrix_n = get_intensities_matrix(red_unique_matrix_n, red_intervals_n_term,  psms, "n-term") 

t_red_unique_matrix_c = red_unique_matrix_c.transpose()
t_red_unique_matrix_n = red_unique_matrix_n.transpose()
t_intensity_matrix_c = intensity_matrix_c.transpose()
t_intensity_matrix_n = intensity_matrix_n.transpose()

print("Cterm")
print(red_intervals_c_term)
print(red_unique_matrix_c)
print(intensity_matrix_c)
print("Nterm")
print(red_intervals_n_term)
print(red_unique_matrix_n)
print(intensity_matrix_n)

print("transposed, and intensities")

print(t_red_unique_matrix_c)
print(t_intensity_matrix_c )
print(t_red_unique_matrix_n )
print(t_intensity_matrix_n)

def get_equation_system(unique_matrix_t, intensity_matrix_t):

    var = 0
    equations = []
    variables = []
    for row in range(unique_matrix_t.shape[0]):
        eq = unique_matrix_t[row,:]
        unique_grps = np.unique(eq).tolist()

        for u in unique_grps[:-1]:
            retaining_ratios = [v==u for v in eq]

            print("............")

            print(retaining_ratios)

            #compute X
            sum_intens_at_loc = intensity_matrix_t[row,:]
            print(sum_intens_at_loc)
            sum_intens_at_loc = np.unique(sum_intens_at_loc).tolist()
            print(sum_intens_at_loc)
            sum_intens_at_loc = sum(np.unique(sum_intens_at_loc).tolist())


            intens_proteos_at_loc = list(compress(intensity_matrix_t[row,:], retaining_ratios))
            print(intens_proteos_at_loc)
            intens_proteos_at_loc = np.unique(intens_proteos_at_loc).tolist()
            print(intens_proteos_at_loc)
            intens_proteos_at_loc  = sum(intens_proteos_at_loc)


            var = intens_proteos_at_loc/sum_intens_at_loc

            print("{0} / {1} = {2}".format(intens_proteos_at_loc,intens_proteos_at_loc,var))

            variables.append(var) 
            equations.append(  1*np.array(retaining_ratios)   )

    return equations, variables


equations_c, variables_c = get_equation_system(t_red_unique_matrix_c, t_intensity_matrix_c)

print(equations_c)
print(variables_c)

#  np.linalg.lstsq(a, b)

# mod_sumix__matrc = mod_sum_matrix_n = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1]))
# unique_matrix_c = unique_matrix_n = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1]))


# #fill mod_sum_matrix from left (cterm)
# mod_mass_sum_c = np.zeros(mod_sum_matrix_c.shape[0])
# for c in range(unique_matrix_c.shape[1]):
#     print(mod_mass_sum_c)
#     mod_mass_sum_c = np.add(mod_mass_sum_c, mod_sum_matrix_c[:, c])
#     print(mod_mass_sum_c)
#     

#         fragments_at_range = psms[r].get_fragments_at_range( intervals_c_term[c], direction="c-term" )
#         mod_sum_matrix_c[r,c] = S.get_sum_intens_fragment_list(psms[r], fragments_at_range)

#         print("{0} has {1} annotated fragments between pos {2} andd {3}".format(psms[r].get_modification_brno(), len(fragments_at_range),intervals_c_term[c][0], intervals_c_term[c][1] ))
        

# print(mod_sum_matrix_c)
# print(unique_matrix_c)





# print("Transposed")
# unique_matrix_t = unique_matrix.transpose()
# unique_frag_sum_matrix_t = unique_frag_sum_matrix.transpose()
# print(unique_matrix_t)
# print(unique_frag_sum_matrix_t)




