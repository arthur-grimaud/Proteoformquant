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
from scipy.optimize import nnls
import itertools

with open('pfq_out_obj_test_4b.pkl', 'rb') as inp:
    exp = pickle.load(inp) 


#S = exp.spectra["scan=2779"]
#S = exp.spectra["index=1853"]
#S = exp.spectra["index=989"]
S = exp.spectra["index=1937"]





def group_matrix(A):
    "From a matrix A return a matrix M with indices corresponding to identical groups column wise"
    M =  np.zeros((A.shape[0], A.shape[1]))
    for col in range(A.shape[1]):
        grp_mass = []
        for row, m in enumerate(A[:,col]):
            if m not in grp_mass:
                grp_mass.append(m)
            grp_index = grp_mass.index(m)
            M[row,col] = grp_index + 1
    return M


def reduce_intervals(  unique_matrix, intervals):
    """From a unique/group matrix where columns correspond to intervals in "intervals", 
    merge or delete intervals that are uninformative on the ratios of proteoforms """
    n_col = unique_matrix.shape[1]
    n_row = unique_matrix.shape[0]

    #first pass to remove interval without any unique ions acroos proteoforms
    col_rm = []
    for col in range(unique_matrix.shape[1]):
        if (unique_matrix[:,col] == np.ones(n_row)).all(): #if no unique ions in interval
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


def get_intensity_matrix( unique_matrix, intervals, psms, direction):
    intensity_matrix =  np.ones((unique_matrix.shape[0], unique_matrix.shape[1]))
    for col in range(unique_matrix.shape[1]):
        for row in range(unique_matrix.shape[0]):
            fragments_at_range = psms[row].get_fragments_at_range(intervals[col], direction= direction)
            intensity_matrix[row,col] = psms[row].spectrum.get_sum_intens_fragment_list(psms[row], fragments_at_range)

    return intensity_matrix


def merge_equations(  A,B):
    """ from matrices A and B of the same dimenssion return A_out and B_out 
    where identical rows in A_out have been merged and corresponding rows in B_out have been added"""

    y, inverse = np.unique(A, axis=0, return_inverse=True)

    A_out = np.empty((0, A.shape[1]), int)
    B_out = np.empty((0, A.shape[1]), int)

    for g in np.unique(inverse):
        indices_identical = [i for i, x in enumerate(inverse) if x == g]
        A_out = np.append(A_out, np.array( [A[indices_identical[0],:]] ), axis=0 ) #add first group row from unique matrix
        B_out = np.append(B_out, np.array( [B[indices_identical[0],:]] ), axis=0 ) #add first group row from intensity matrix

        for i in range(1,len(indices_identical)): #sum up intensities of identical rows
            B_out[-1,:] = np.sum( [ B_out[-1,:], B[indices_identical[i],:] ], axis=0)

    return A_out, B_out

def equation_system(  unique_matrix_t, intensity_matrix_t):

    var = 0
    equations = []
    variables = []
    for row in range(unique_matrix_t.shape[0]):
        eq = unique_matrix_t[row,:]
        unique_grps = np.unique(eq).tolist()

        for u in unique_grps[:-1]:
            retaining_ratios = [v==u for v in eq]

            sum_intens_at_loc = intensity_matrix_t[row,:]
            sum_intens_at_loc = np.unique(sum_intens_at_loc).tolist()
            sum_intens_at_loc = sum(np.unique(sum_intens_at_loc).tolist())

            intens_proteos_at_loc = list(compress(intensity_matrix_t[row,:], retaining_ratios))
            intens_proteos_at_loc = np.unique(intens_proteos_at_loc).tolist()
            intens_proteos_at_loc  = sum(intens_proteos_at_loc)

            if sum_intens_at_loc != 0:
                var = intens_proteos_at_loc/sum_intens_at_loc
            else:
                var = 0

            variables.append(var) 
            equations.append(  1*np.array(retaining_ratios)   )


    variables.append(1)
    equations.append(np.ones(unique_matrix_t.shape[1]))

    return equations, variables



def ratios_multiple(psms):
    print("****START****")
    #Position (1 based) of modified residues across multiple psm (e.g K14acK24ac K9me3 -> [9,14,24])
    mod_pos_list = [] 
    for psm in psms:
        print(psm.get_modification_brno())
        for mod in psm.proteoform.get_modification_dict():
            #print(mod)
            if mod["location"] not in mod_pos_list:
                mod_pos_list.append(mod["location"])

    mod_pos_list = sorted(mod_pos_list)
    #print("::::")
    #print(mod_pos_list)
   
    proteforms_multiquant = [psm.get_modification_brno() for psm in psms]

    #Add start and end of the sequence, sort the postions and create "interval of interest" for both n and c term frags 
    mod_pos_list_n = sorted([*mod_pos_list, 1, len(psms[0].proteoform.get_peptide_sequence())+1])
    mod_pos_list_c = sorted([*mod_pos_list, 0, len(psms[0].proteoform.get_peptide_sequence())])
    intervals_n = [[mod_pos_list_n[p],mod_pos_list_n[p+1]-1] for p in range(len(mod_pos_list_n)-1)]
    intervals_c = [[mod_pos_list_c[p]+1,mod_pos_list_c[p+1]] for p in range(len(mod_pos_list_c)-1)]

    # print("::::")
    # print(intervals_n)
    # print("::::")
    # print(intervals_c)
    #Create matrix with modification induced mass shift (row = proteoform, col position)
    mod_mass_matrix = np.zeros( (len(psms), len(mod_pos_list)) )
    for row, psm in enumerate(psms):
        #print(psm.get_modification_brno())
        for mod in psm.proteoform.get_modification_dict():
            #print(mod)
            col = mod_pos_list.index(mod["location"])
            mod_mass_matrix[row,col] = mod["monoisotopicMassDelta"]


    # print("***********")
    # print(mod_mass_matrix)
    # print("***********")
    #Compute matrix with cumulative mass shift for n and c term fragments ( row: proteoform, column:corresponding position interval)
    
    additive_mass_n = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1]+1))
    for col in range(1, additive_mass_n.shape[1]):
        additive_mass_n[:,col] = np.add( additive_mass_n[:,col-1], mod_mass_matrix[:,col-1])

    additive_mass_c = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1]+1))
    for col in range(additive_mass_c.shape[1]-1, 0, -1):
        additive_mass_c[:,col-1] = np.add( additive_mass_c[:,col], mod_mass_matrix[:,col-1])
    

    # print("-----------------------------")

    # print(additive_mass_n)
    # print(additive_mass_c)

    # print("-----------------------------")


    #Compute grouping matrix, give a common index to intervals with same mass across proteoforms (e.g [[54,54][54,59]] = [[1,1][1,2]])
    unique_matrix_n =  group_matrix(additive_mass_n)  
    unique_matrix_c =  group_matrix(additive_mass_c)


    # print("----------||||-----------")

    # print(unique_matrix_n)
    # print(unique_matrix_c)

    # print("------------|||||------------")

    #reduce intervales
    red_unique_matrix_n, red_intervals_n = reduce_intervals(unique_matrix_n, intervals_n)
    red_unique_matrix_c, red_intervals_c = reduce_intervals(unique_matrix_c, intervals_c)

    #compute intensity matrices
    intensity_matrix_n = get_intensity_matrix(red_unique_matrix_n, red_intervals_n,  psms, "n-term") 
    intensity_matrix_c = get_intensity_matrix(red_unique_matrix_c, red_intervals_c,  psms, "c-term") 


    # print("..................")
    # print(red_unique_matrix_n)
    # print(red_intervals_n)
    # print(red_unique_matrix_c)
    # print(red_intervals_c)
    # print("..................")
    
    #transpose
    t_red_unique_matrix_n = red_unique_matrix_n.transpose()
    t_red_unique_matrix_c = red_unique_matrix_c.transpose()
    t_intensity_matrix_n = intensity_matrix_n.transpose()
    t_intensity_matrix_c = intensity_matrix_c.transpose()
    


    #Combine matrices
    unique_matrix = np.vstack((t_red_unique_matrix_c, t_red_unique_matrix_n))
    intensity_matrix = np.vstack((t_intensity_matrix_c, t_intensity_matrix_n))
    #save matrices: 
    unique_matrix = unique_matrix
    intensity_matrix = intensity_matrix

    #print(unique_matrix)
    #print(intensity_matrix)

    #reduce equations that are redundant
    unique_matrix, intensity_matrix = merge_equations(unique_matrix, intensity_matrix)

    #print(unique_matrix)
    #print(intensity_matrix)

    #get eq system
    equations, variables = equation_system(unique_matrix, intensity_matrix)

    #print(equations)
    #print(variables)

    #Non linear least square solving:
    results = nnls(equations, variables)
    ratios_psms = [ratio/sum(results[0]) for ratio in results[0]] #normalized ratios
    quant_residuals = results[1]

    return quant_residuals
    # #assign rations: 
    # for idx, psm in enumerate(psms):
    #     psm.ratio = ratios_psms[idx]





psms = S.psms

combinations = list(itertools.combinations(range(0,len(psms)-1), 3))

for comb in combinations:

    print(comb)
    psms = [psm for i,psm in enumerate(S.psms) if i in comb]
    print(ratios_multiple(psms))