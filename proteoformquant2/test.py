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
import itertools

with open('pfq_out_obj_test_1b.pkl', 'rb') as inp:
    exp = pickle.load(inp) 


#S = exp.spectra["scan=2779"]
#S = exp.spectra["index=1853"]
#S = exp.spectra["index=989"]
S = exp.spectra["index=3151"]


psms = S.psms[0:5]

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

intervals_n = [[mod_pos_list_n[p],mod_pos_list_n[p+1]-1] for p in range(len(mod_pos_list_n)-1)]
intervals_c = [[mod_pos_list_c[p]+1,mod_pos_list_c[p+1]] for p in range(len(mod_pos_list_c)-1)]

print("positions intervals N")
print(intervals_n)
print("positions intervals C")
print(intervals_c)


mod_mass_matrix = np.zeros( (len(psms), len(mod_pos_list)) )
for row, psm in enumerate(psms):
    #print(psm.get_modification_brno())
    for mod in psm.proteoform.get_modification_dict():
        #print(mod)
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
print("Nterm")
print(additive_mass_n)
print("Cterm")
print(additive_mass_c)


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


print("Grouping matrix N and C term")
print("Nterm")
print(intervals_n)
print(unique_matrix_n)
print("Cterm")
print(intervals_c)
print(unique_matrix_c)


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
        
    


red_unique_matrix_n, red_intervals_n = reduce_unique_matrix(unique_matrix_n, intervals_n)
red_unique_matrix_c, red_intervals_c = reduce_unique_matrix(unique_matrix_c, intervals_c)



def get_intensities_matrix(unique_matrix, intervals, psms, direction):
    intensity_matrix =  np.ones((unique_matrix.shape[0], unique_matrix.shape[1]))
    for col in range(unique_matrix.shape[1]):
        #print("new col")
        for row in range(unique_matrix.shape[0]):
            

            fragments_at_range = psms[row].get_fragments_at_range(intervals[col], direction= direction)
            # print("dir = {4} GRP = {1} / {2}, range used: {3} Number of Possible fragments: {0} frag found:".format(len(fragments_at_range), 
            #                                                                                                         unique_matrix[row,col],
            #                                                                                                         psms[row].get_modification_brno(),
            #                                                                                                         intervals[col],
            #                                                                                                         direction))
            intensity_matrix[row,col] = psms[row].spectrum.get_sum_intens_fragment_list(psms[row], fragments_at_range)

    return intensity_matrix

intensity_matrix_c = get_intensities_matrix(red_unique_matrix_c, red_intervals_c,  psms, "c-term") 
intensity_matrix_n = get_intensities_matrix(red_unique_matrix_n, red_intervals_n,  psms, "n-term") 

t_red_unique_matrix_c = red_unique_matrix_c.transpose()
t_red_unique_matrix_n = red_unique_matrix_n.transpose()
t_intensity_matrix_c = intensity_matrix_c.transpose()
t_intensity_matrix_n = intensity_matrix_n.transpose()


print("### reduced matrices and intensity  ###")
print("Nterm")
print(red_intervals_n)
print(red_unique_matrix_n)
print(intensity_matrix_n)
print("Cterm")
print(red_intervals_c)
print(red_unique_matrix_c)
print(intensity_matrix_c)


print("### Transposed intensities C and N term ###")
print("Nterm")
print(t_red_unique_matrix_n )
print(t_intensity_matrix_n)
print("Cterm")
print(t_red_unique_matrix_c)
print(t_intensity_matrix_c )



print("### Combined Matrices ###")

unique_matrix = np.vstack((t_red_unique_matrix_c, t_red_unique_matrix_n))
intensity_matrix = np.vstack((t_intensity_matrix_c, t_intensity_matrix_n))

print(unique_matrix)
print(intensity_matrix)

def merge_unique_intensity_matrix(A,B):
    #print("### In the merging function  ###")

    y, inverse = np.unique(A, axis=0, return_inverse=True)

    A_out = np.empty((0, A.shape[1]), int)
    B_out = np.empty((0, A.shape[1]), int)

    # print(inverse)

    for g in np.unique(inverse):
        indices_identical = [i for i, x in enumerate(inverse) if x == g]
        #print("group: ", g , "indices", indices_identical)

        A_out = np.append(A_out, np.array( [A[indices_identical[0],:]] ), axis=0 ) #add first group row from unique matrix
        B_out = np.append(B_out, np.array( [B[indices_identical[0],:]] ), axis=0 ) #add first group row from intensity matrix

        # print("f b out")
        # print(B_out)

        for i in range(1,len(indices_identical)): #sum up intensities of identical rows
            # print(B_out[-1,:])
            # print(B[indices_identical[i],:])
            # print(np.sum( [ B_out[-1,:], B[indices_identical[i],:] ], axis=0))
            B_out[-1,:] = np.sum( [ B_out[-1,:], B[indices_identical[i],:] ], axis=0)

    return A_out, B_out

unique_matrix, intensity_matrix = merge_unique_intensity_matrix(unique_matrix, intensity_matrix)

print("### Removed redundancy in equations ###")
print(unique_matrix)
print(intensity_matrix)

def get_equation_system(unique_matrix_t, intensity_matrix_t):

    var = 0
    equations = []
    variables = []
    for row in range(unique_matrix_t.shape[0]):
        eq = unique_matrix_t[row,:]
        unique_grps = np.unique(eq).tolist()

        for u in unique_grps[:-1]:
            retaining_ratios = [v==u for v in eq]

            #print("............")

            #print(retaining_ratios)

            #compute X
            sum_intens_at_loc = intensity_matrix_t[row,:]
            #print(sum_intens_at_loc)
            sum_intens_at_loc = np.unique(sum_intens_at_loc).tolist()
            #print(sum_intens_at_loc)
            sum_intens_at_loc = sum(np.unique(sum_intens_at_loc).tolist())


            intens_proteos_at_loc = list(compress(intensity_matrix_t[row,:], retaining_ratios))
            #print(intens_proteos_at_loc)
            intens_proteos_at_loc = np.unique(intens_proteos_at_loc).tolist()
            # print(intens_proteos_at_loc)
            intens_proteos_at_loc  = sum(intens_proteos_at_loc)

            if sum_intens_at_loc != 0:
                var = intens_proteos_at_loc/sum_intens_at_loc
            else:
                var = 0

            #print("{0} / {1} = {2}".format(intens_proteos_at_loc,intens_proteos_at_loc,var))

            variables.append(var) 
            equations.append(  1*np.array(retaining_ratios)   )


    #append sum all proteo equal to 1
    variables.append(1)
    equations.append(np.ones(unique_matrix_t.shape[1]))


    return equations, variables


equations, variables = get_equation_system(unique_matrix, intensity_matrix)


print("### Final Equation system ###")
print(np.array(equations))
print(variables)

from scipy.optimize import nnls

print("### Result Least Square (non-negative) ###")
print(nnls(equations, variables)[0])
print("sum ratios")
print(sum(nnls(equations, variables)[0]))
print("Residuals")
print(nnls(equations, variables)[1])






