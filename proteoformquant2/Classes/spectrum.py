from array import array
from Classes.psm import Psm

from itertools import compress
from scipy.optimize import nnls
import numpy as np
import pprint

class Spectrum():

    def __init__(self, spectrumID, identMzid = None):

     
        self.id = spectrumID  #Unique ID for the spectrum

        #print(identMzid)

        try:
            self.spectrum_title = identMzid['spectrum title']
        except(KeyError):
            self.spectrum_title = identMzid['spectrumID']

        self.psms : list(Psm) = []  #list of PSM for that spectrum
        
        if(identMzid != None):
            self.set_ident_data_mzid(identMzid)

        self.sumIntensAnnotFrag = 0

        self.quant_residuals = 0

        #MultiQuant Attr (will be removed)

        self.unique_matrix = None
        self.intensity_matrix = None
        self.unique_matrix_r = None
        self.intensity_matrix_r = None
        self.equations = None
        self.variables = None
        self.proteforms_multiquant = None
        self.intervals = None
            
        pass
    

    #Getters:

    def get_id(self):
        return self.id

    def getPrecMz(self):
        return self.precMz
    
    def getPrecIntens(self):
        return self.precIntens

    def get_rt(self):
        return self.rt

    def getFragIntens(self):
        return self.fragIntens

    def getFragMz(self):
        return self.fragMz

    def getSumIntensFrag(self):
        return sum(self.fragIntens)

    def get_sum_intens_frag_at_mz(self, mz_range):
        #return the sum of intensity at specified mz range

        intensities = [pair[1] for pair in zip(self.fragMz,self.fragIntens) if pair[0] > mz_range[0] and pair[0] < mz_range[1]]
        print(intensities)
        return sum(intensities)
    
    def getSumIntensAnnotFrag(self):
        return self.sumIntensAnnotFrag

    def get_validated_psm(self):
        return [psm for psm in self.psms if psm.isValidated == True]

    def get_number_validated_psm(self):
        return len([psm for psm in self.psms if psm.isValidated == True])


    #Setters


    def set_ident_data_mzid(self, identMzid):
        """fill PSM list from identification result dict from pyteomics"""
        for identItem in identMzid["SpectrumIdentificationItem"]: #Iterate over identification item and create an instance of the object psm for each
            #print("*************PSM************")
            #print(identItem)
            self.psms.append(Psm(rank = len(self.psms)+1, spectrum=self, identificationItem = identItem))


    def set_spec_data_mgf(self, specMgf):
        "add spetrum information from a pyteomics mgf object"
        self.fragIntens: array = specMgf["intensity array"]
        self.fragMz: array = specMgf["m/z array"]
        self.precIntens: float = specMgf["params"]["pepmass"][1]
        self.precMz: float = specMgf["params"]["pepmass"][0]
        self.rt: float = specMgf["params"]["rtinseconds"]
        
    def set_spec_data_mzml(self, spec_mzml):
        "add spetrum information from a pyteomics mzml object"
        print(spec_mzml['scanList']['scan'][0]['scan start time'])

        self.fragIntens: array = spec_mzml["intensity array"]
        self.fragMz: array = spec_mzml["m/z array"]
        self.precIntens: float = spec_mzml['precursorList']['precursor'][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
        self.precMz: float = spec_mzml['precursorList']['precursor'][0]["selectedIonList"]["selectedIon"][0]["peak intensity"]
        self.rt: float = spec_mzml['scanList']['scan'][0]['scan start time']

    def setSumIntensAnnotFrag(self):
        """For proteoforms in self.psms where isvalidated is true, set self.sumIntensAnnotFrag as the sum of annotated peaks"""
        self.sumIntensAnnotFrag = 0
        seenPeaks=[]
        for psm in self.psms:
            if psm.isValidated:
                ann = psm.get_annotation()
                for type in ann.values():
                    for i in range(len(type["intens"])):
                        if type["index"] not in seenPeaks:
                            seenPeaks.append(i)
                            self.sumIntensAnnotFrag += type["intens"][i]
        #print("sum intensities annotated: " + str(self.sumIntensAnnotFrag))



    #other methods:

    def annotateFragPsm(self, frag_mz_tol):
        for psm in self.psms:
            psm.setAnnotatedFragments(frag_mz_tol = frag_mz_tol)
            
    def updateValidation(self):
        pass


    # --------------------------- PSMs Ratio in Spectra -------------------------- #

    def update_psms_ratio(self, psms_rank = []):
        """compute the ratios of psms in "psms"
        if psms_rank is empty proteoforms that are validated (isvalidated == T) are considered
        else the proteforms of rank in psms_rank are considered"""

        if len(psms_rank) == 0:
            psms = self.get_validated_psm()
        else:
            psms = [psm for idx, psm in self.psms if idx in psms_rank]
        #print(psms)

        if len(psms) > 0:
            self.__ratios_multiple(psms)
        else:
            print("No valid PSM")
        #self.__ratios_pair(get_validated_psm())


    def __ratios_multiple(self, psms):
       
        #Position (1 based) of modified residues across multiple psm (e.g K14acK24ac K9me3 -> [9,14,24])
        mod_pos_list = [] 
        for psm in psms:
            for mod in psm.proteoform.get_modification_dict():
                if mod["location"] not in mod_pos_list:
                    mod_pos_list.append(mod["location"])
        mod_pos_list = sorted(mod_pos_list)
        self.proteforms_multiquant = [psm.get_modification_brno() for psm in psms]

        #Add start and end of the sequence, sort the postions and create "interval of interest" for both n and c term frags 
        mod_pos_list_n = sorted([*mod_pos_list, 1, len(psms[0].proteoform.get_peptide_sequence())+1])
        mod_pos_list_c = sorted([*mod_pos_list, 0, len(psms[0].proteoform.get_peptide_sequence())])
        intervals_n = [[mod_pos_list_n[p],mod_pos_list_n[p+1]-1] for p in range(len(mod_pos_list_n)-1)]
        intervals_c = [[mod_pos_list_c[p]+1,mod_pos_list_c[p+1]] for p in range(len(mod_pos_list_c)-1)]


        #Create matrix with modification induced mass shift (row = proteoform, col position)
        mod_mass_matrix = np.zeros( (len(psms), len(mod_pos_list)) )
        for row, psm in enumerate(psms):
            #print(psm.get_modification_brno())
            for mod in psm.proteoform.get_modification_dict():
                #print(mod)
                col = mod_pos_list.index(mod["location"])
                mod_mass_matrix[row,col] = mod["monoisotopicMassDelta"]

        #Compute matrix with cumulative mass shift for n and c term fragments ( row: proteoform, column:corresponding position interval)
        additive_mass_n = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1]+1))
        for col in range(additive_mass_n.shape[1]-1, 0, -1):
            additive_mass_n[:,col-1] = np.add( additive_mass_n[:,col], mod_mass_matrix[:,col-1])
        
        additive_mass_c = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1]+1))
        for col in range(1, additive_mass_c.shape[1]):
            additive_mass_c[:,col] = np.add( additive_mass_c[:,col-1], mod_mass_matrix[:,col-1])

        
        #Compute grouping matrix, give a common index to intervals with same mass across proteoforms (e.g [[54,54][54,59]] = [[1,1][1,2]])
        unique_matrix_n =  self.__group_matrix(additive_mass_n)  
        unique_matrix_c =  self.__group_matrix(additive_mass_c)

        #reduce intervales
        red_unique_matrix_n, red_intervals_n = self.__reduce_intervals(unique_matrix_n, intervals_n)
        red_unique_matrix_c, red_intervals_c = self.__reduce_intervals(unique_matrix_c, intervals_c)

        #compute intensity matrices
        intensity_matrix_n = self.__intensity_matrix(red_unique_matrix_n, red_intervals_n,  psms, "n-term") 
        intensity_matrix_c = self.__intensity_matrix(red_unique_matrix_c, red_intervals_c,  psms, "c-term") 
        
        #transpose
        t_red_unique_matrix_n = red_unique_matrix_n.transpose()
        t_red_unique_matrix_c = red_unique_matrix_c.transpose()
        t_intensity_matrix_n = intensity_matrix_n.transpose()
        t_intensity_matrix_c = intensity_matrix_c.transpose()
        

        #Combine matrices
        unique_matrix = np.vstack((t_red_unique_matrix_c, t_red_unique_matrix_n))
        intensity_matrix = np.vstack((t_intensity_matrix_c, t_intensity_matrix_n))
        #save matrices: 
        self.unique_matrix = unique_matrix
        self.intensity_matrix = intensity_matrix
    
        #reduce equations that are redundant
        unique_matrix, intensity_matrix = self.__merge_equations(unique_matrix, intensity_matrix)
        #save matrices
        self.unique_matrix_r = unique_matrix
        self.intensity_matrix_r = intensity_matrix

        #get eq system
        equations, variables = self.__equation_system(unique_matrix, intensity_matrix)

        self.equations = unique_matrix
        self.variables = intensity_matrix

        #Non linear least square solving:
        results = nnls(equations, variables)
        ratios_psms = [ratio/sum(results[0]) for ratio in results[0]] #normalized ratios
        self.quant_residuals = results[1]

        #assign rations: 
        for idx, psm in enumerate(psms):
            psm.ratio = ratios_psms[idx]


        


    def __group_matrix(self, A):
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


    def __reduce_intervals(self, unique_matrix, intervals):
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


    def __intensity_matrix(self, unique_matrix, intervals, psms, direction):
        intensity_matrix =  np.ones((unique_matrix.shape[0], unique_matrix.shape[1]))
        for col in range(unique_matrix.shape[1]):
            for row in range(unique_matrix.shape[0]):
                fragments_at_range = psms[row].get_fragments_at_range(intervals[col], direction= direction)
                intensity_matrix[row,col] = psms[row].spectrum.get_sum_intens_fragment_list(psms[row], fragments_at_range)

        return intensity_matrix


    def __merge_equations(self, A,B):
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

    def __equation_system(self, unique_matrix_t, intensity_matrix_t):

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





    # ------------------ quantification of two proteoforms only ------------------ #


    def __ratios_pair(self, psms):
        #TODO works only for 2 validated psms

        validatedPsms = psms

        uniquePairwise = [] # list of list of unique ions names formatted as such [[List Unique validatedPsms[0 and 1]], [List Unique validatedPsms[1 and 2]], .... ]
        ratioPairwise = []

        if len(validatedPsms) == 1:
            validatedPsms[0].ratio = 1.0
        else:
            for p in range(0, len(validatedPsms)-1):
                uniqueFragments = []
                for fragType in validatedPsms[p].proteoform.theoFrag.keys():
                    for fragment in validatedPsms[p].proteoform.theoFrag[fragType]:
                        if validatedPsms[p].proteoform.theoFrag[fragType][fragment] != validatedPsms[p+1].proteoform.theoFrag[fragType][fragment]:
                            uniqueFragments.append(fragment)
                uniquePairwise.append(uniqueFragments) 

            for p in range(0, len(validatedPsms)-1):

                A = self.get_sum_intens_fragment_list(validatedPsms[p],uniquePairwise[p])
                B = self.get_sum_intens_fragment_list(validatedPsms[p+1],uniquePairwise[p])

                validatedPsms[p].ratio = A/(A+B)
                validatedPsms[p+1].ratio = B/(A+B)  


    # def update_psms_ratio(self):

    #     validatedPsms = self.get_validated_psm()


        
    #     uniquePairwise = [] # list of list of unique ions names formatted as such [[List Unique validatedPsms[0 and 1]], [List Unique validatedPsms[1 and 2]], .... ]
    #     ratioPairwise = []

    #     if len(validatedPsms) == 1: #Only one validated proteoform
    #         validatedPsms[0].ratio = 1.0

    #     elif len(validatedPsms) == 2: #Pair of validated proteoform
    #         uniques = self.__get_unique_fragments(validatedPsms[0], validatedPsms[1])
    #         A = self.get_sum_intens_fragment_list(validatedPsms[0],uniques)
    #         B = self.get_sum_intens_fragment_list(validatedPsms[1],uniques)
    #         validatedPsms[0].ratio = A/(A+B)
    #         validatedPsms[1].ratio = B/(A+B)  

    #     else: #More than 2 validated proteoforms
            
    #         for p_A in range(0, len(validatedPsms)):
    #             p_B = (p_A + 1) % len(validatedPsms)
                
    #             uniques = self.__get_unique_fragments(validatedPsms[p_A], validatedPsms[p_B])

    #             A = self.get_sum_intens_fragment_list(validatedPsms[p_A],uniques) 
    #             B = self.get_sum_intens_fragment_list(validatedPsms[p_B],uniques) 
                
                

    #             self.__get_unique_fragments()
    #             uniquePairwise.append(uniqueFragments) 

    #         for p in range(0, len(validatedPsms)-1): 

    #             A = self.get_sum_intens_fragment_list(validatedPsms[p],uniquePairwise[p])
    #             B = self.get_sum_intens_fragment_list(validatedPsms[p+1],uniquePairwise[p])

    #             validatedPsms[p].ratio = A/(A+B)
    #             validatedPsms[p+1].ratio = B/(A+B)  



    def __get_unique_fragments(self, proteo_a, proteo_b):
        "from two proteoform object get the list of unique fragments (frag with mass shift)"
        uniqueFragments = []
        for fragType in validatedPsms[p].proteoform.theoFrag.keys():
            for fragment in validatedPsms[p].proteoform.theoFrag[fragType]:
                if validatedPsms[p].proteoform.theoFrag[fragType][fragment] != validatedPsms[p+1].proteoform.theoFrag[fragType][fragment]:
                    uniqueFragments.append(fragment)

        return uniqueFragments

    def get_sum_intens_fragment_list(self, psm, fragments):
        #could be moved to psm object
        intensities = []
        for fragType in psm.annotation.values():
            for i, fragCode in enumerate(fragType["fragCode"]):
                if fragCode in fragments:
                    intensities.append(int(fragType["intens"][i]))

        if len(intensities) == 0:
            return 0
        else:
            return(sum(intensities))



    