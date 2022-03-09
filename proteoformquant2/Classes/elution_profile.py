
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
from scipy import integrate,stats
import numpy as np
import math
import warnings
from sklearn.linear_model import LinearRegression
from kneed import KneeLocator
from scipy.integrate import trapz, simps
from statistics import mean, median, stdev
from scipy.optimize import least_squares
from sympy import subsets
import sympy as sy


class ElutionProfile():

    def __init__(self):
        
        #Parameters of the function 
        self.param_estimated = [None]
        self.param_fitted = [None]

        #Goodness of fit 
        self.score_estimated = None
        self.score_fitted = None

        #List of PSMs
        self.psms_included = [] #all psms used to compute the elution profile
        self.psms_outliers = [] #psms considered as outliers in the modelling of the elution profile

        #Score threshold
        self.score_threshold = None

        pass
    
    
    # ---------------------------------- Getters --------------------------------- #

    def get_y_serie(self, x_values, method = "best"):
        """
        for a list of x values return the estimated y values given the fitted function
        use estimated parameters if fitted parameters are not determined or is method = best
        """
        if self.param_fitted[0] != None and method in ["fitted","best"]:
            return self.__skewnormal(x_values, *self.param_fitted), self.param_fitted
        elif (method in ["estimated","best"]):
            return self.__skewnormal(x_values, *self.param_estimated), self.param_estimated
        else:
            return [None], [None]

    def get_y(self, x, method = "best"):
        """
        for x return the estimated y value given the fitted function
        uses estimated parameters if fitted parameters are not determined
        """
        if self.param_fitted[0] != None and method in ["fitted","best"]:
            return self.__skewnormal(x, *self.param_fitted)
        elif (method in ["estimated","best"]):
            return self.__skewnormal(x, *self.param_estimated)
        else:
            return None

    def is_parameters_fitted(self):
        #return true if curve_fit results is available in the instance
        if self.param_fitted[0] == None:
            return False
        else:
            return True

    def get_parameters_fitted(self):
        #returns parameters of curve_fit

        return self.param_fitted

    def get_elution_profile_param_m(self, method):
        if method == "fitted":
            return self.param_estimated[0]
        elif method == "estimated":
            return self.parametersFitted[0]
        else: 
            print("Specify method (fitted or estimated)")

    def get_elution_profile_std(self, method ): #might be incorrect
        if method == "fitted":
            return self.param_estimated[1]
        elif method == "estimated":
            return self.parametersFitted[1]
        else:
            print("Specify method (fitted or estimated) ")

    def get_elution_profile_param_k(self, method): 
        if method == "fitted":
            return self.param_estimated[1]
        elif method == "estimated":
            return self.parametersFitted[1]
        else:
            print("Specify method (fitted or estimated) ")

    def get_auc(self, x_a, x_b):

        try:
            return integrate.quad(lambda x: self.__skewnormal(x, *self.param_fitted), x_a, x_b)[0]
        except(TypeError):
            return 0

    def get_x_at_max_y(self):
        """ Returns the x value for the peak of the elution profile  COULD EXIST A MATHEMATICAL SOLUTION"""
        m = self.get_elution_profile_param_m("fitted")
        k = self.get_elution_profile_param_k("fitted")

        y = self.get_y(m)
         
        search = True
        x = m

        if k > 0:
            while search:
                if y < self.get_y(x+1) :
                    x = x+1
                    y = self.get_y(x)
                else:
                    search = False

        elif k < 0:
            while search:
                if y < self.get_y(x-1) :
                    x = x-1
                    y = self.get_y(x)
                else:
                    search = False

        return x



        

    # ---------------------------------- Fitting --------------------------------- #


    def model_elution_profile(self, psms, score_threshold):

        self.psms_included  = psms
        self.score_threshold = score_threshold 

        self.data_x = np.array([psm.spectrum.get_rt() for psm in psms])
        self.data_y = np.array([psm.get_prec_intens_ratio() for psm in psms])
        
        #Fit model without outliers removal
        self.param_estimated, self.param_fitted, self.score_estimated, self.score_fitted = self.fit_skew_normal(self.data_x, self.data_y)
        
        # Fit model with outliers removal if below score threshold
        if self.score_fitted < self.score_threshold:
            self.param_estimated, self.param_fitted, self.score_estimated, self.score_fitted, self.psms_outliers, self.psms_included = self.exclude_outliers_mean_method()

        
        #Add to outliers:
        self.psms_outliers, self.psms_included = self.exclude_outlier_non_significant(0.5)



    def fit_skew_normal(self, data_x, data_y, param_init = [None], param_bounds = [None]):
        """ startPrevFit """

        #Transform to numpy array (might not be needed)
        data_x = np.array(data_x)
        data_y = np.array(data_y)


        #Determine parameters bounds
        if param_bounds[0] == None: #If no parameters bounds are suplied, estimate them
            param_bounds = self.__get_parameter_bounds_skewnormal(data_x, data_y)

        #Determine initial parameters
        if param_init[0] == None: #If no initial parameters are suplied, estimate them
            param_init = self.__estimate_initial_parameters(self.__skewnormal,data_x, data_y)

        #Save initial parameters as "estimated"
        param_estimated = param_init
        #goodness of the fit with estimated parameters
        score_estimated = self.__pearson_test(self.__skewnormal, param_estimated, data_x, data_y)


        try:

            param_bounds_curve_fit = tuple([ tuple([param_bounds[x][b] for x in range(len(param_bounds))]) for b in range(len(param_bounds[0]))]) #convert bounds for "curve_fit"
            param_fitted = least_squares(self.__skewnormal_residuals, param_estimated, bounds = param_bounds_curve_fit, loss='linear', f_scale=1, args=(data_x, data_y))
            param_fitted =param_fitted.x 

            # #Optimize curve fit
            # param_bounds_curve_fit = tuple([ tuple([param_bounds[x][b] for x in range(len(param_bounds))]) for b in range(len(param_bounds[0]))]) #convert bounds for "curve_fit"
            # param_fitted = curve_fit(self.__skewnormal, data_x, data_y, param_estimated, bounds=param_bounds_curve_fit)[0]
            # #goodness of the fit with fitted parameters
            score_fitted = self.__pearson_test(self.__skewnormal, param_fitted, data_x, data_y)
        except(RuntimeError):
            param_fitted = [None]
            score_fitted = 0
            print("Curve_fit failed")


        return param_estimated, param_fitted, score_estimated, score_fitted

    # def fit_skew_normal(self, data_x, data_y, param_init = [None], param_bounds = [None]):
    #     """ startPrevFit """

    #     #Transform to numpy array (might not be needed)
    #     data_x = np.array(data_x)
    #     data_y = np.array(data_y)


    #     #Determine parameters bounds
    #     if param_bounds[0] == None: #If no parameters bounds are suplied, estimate them
    #         param_bounds = self.__get_parameter_bounds_skewnormal(data_x, data_y)

    #     #Determine initial parameters
    #     if param_init[0] == None: #If no initial parameters are suplied, estimate them
    #         param_init = self.__estimate_initial_parameters(self.__skewnormal,data_x, data_y)

    #     #Save initial parameters as "estimated"
    #     param_estimated = param_init
    #     #goodness of the fit with estimated parameters
    #     score_estimated = self.__pearson_test(self.__skewnormal, param_estimated, data_x, data_y)

        
    #     try:
    #         #Optimize curve fit
    #         param_bounds_curve_fit = tuple([ tuple([param_bounds[x][b] for x in range(len(param_bounds))]) for b in range(len(param_bounds[0]))]) #convert bounds for "curve_fit"
    #         param_fitted = curve_fit(self.__skewnormal, data_x, data_y, param_estimated, bounds=param_bounds_curve_fit)[0]
    #         #goodness of the fit with fitted parameters
    #         score_fitted = self.__pearson_test(self.__skewnormal, param_fitted, data_x, data_y)
    #     except(RuntimeError):
    #         param_fitted = [None]
    #         score_fitted = 0
    #         print("Curve_fit failed")


    #     return param_estimated, param_fitted, score_estimated, score_fitted



    # ---------------------------------- Outlier detection --------------------------------- #

    def exclude_outlier_non_significant(self, auc_percent_threshold):
        "Exclude datapoint on the tails of the elution profile model that would represent less than X%  of the area under the curve"

        psms_included = self.psms_included
        psms_outliers = self.psms_outliers

     

        
        for psm in self.psms_included:

            auc_l_r = (self.get_auc(-np.inf, psm.spectrum.get_rt()), self.get_auc(psm.spectrum.get_rt(), np.inf))
            auc_tot = sum(auc_l_r)
            if auc_tot > 0:
                auc_min = min(auc_l_r)
                

                auc_percent = (auc_min/auc_tot)*100

                if auc_percent < auc_percent_threshold/2: #diveded by 2 because two tailed 
                    psms_outliers.append(psm)
            
        psms_included = [psm for psm in self.psms_included if psm not in psms_outliers]

        return psms_outliers, psms_included




    def exclude_outliers_mean_method(self ):
        "Try imrove the fit of the curve by iteratively removing datapoints the furthest from the RT mean"

        subsets_psms = [self.psms_included]
        subsets_rt = [[psm.spectrum.get_rt() for psm in self.psms_included]]
        subsets_index = [0]

        #print(subsets_psms[0][0].get_modification_brno())
        

        #Create a list of psms subsets
        for i in range(0,len(subsets_psms[0])-5): 
            m = mean(subsets_rt[i])
            outIndex = subsets_rt[i].index(max(subsets_rt[i], key=lambda x:abs(x-m)))
            subsets_psms.append( [v for y,v in enumerate(subsets_psms[i]) if y != outIndex] )
            subsets_rt.append( [v for y,v in enumerate(subsets_rt[i]) if y != outIndex] )
            subsets_index.append(i+1)

        #Compute fit score for each subsets
        subsets_scores = []
        for psmsSubset in subsets_psms:
            data_yT = np.array([psm.spectrum.get_rt() for psm in psmsSubset])
            yDataT = np.array([psm.get_prec_intens_ratio() for psm in psmsSubset])
            #refit the curve to subset
            param_estimated, param_fitted, score_estimated, score_fitted = self.fit_skew_normal(data_yT, yDataT)
            #store score and subset


            subsets_scores.append(score_fitted)

        #Find best score that retain the maximum number of psms
        kn = KneeLocator(subsets_index, subsets_scores, S=2, curve='concave', direction='increasing',interp_method= "polynomial", polynomial_degree=2)
        index = kn.knee
        print(index)
        #Return best fit results
        if index != None and index != 0:
            try:
                data_x = np.array([psm.spectrum.get_rt() for psm in subsets_psms[index]])
                data_y = np.array([psm.get_prec_intens_ratio() for psm in subsets_psms[index]])
                param_estimated, param_fitted, score_estimated, score_fitted = self.fit_skew_normal(data_x, data_y)
                psms_outliers = [psm for psm in self.psms_included if psm not in subsets_psms[index]] #TODO check whether outlier needs to be removed for proteoform.linkedPSMSs
                psms_included = [psm for psm in self.psms_included if psm in subsets_psms[index]]
            except(TypeError):
                print("TYPE ERROR: could not optimize fit by removing data points")
                psms_outliers = []
        else: 
            psms_outliers = []
            psms_included = self.psms_included
            print("NO OPTI: could not optimize fit by removing data points")
            pass

        print(subsets_index)
        print(subsets_scores)
        print(score_fitted)
        return param_estimated, param_fitted, score_estimated, score_fitted, psms_outliers, psms_included



    


    # ----------------------------- Model's function ----------------------------- #

    def __skewnormal(self, x, m, s, a, k):
        return a*np.exp( k*(x - m)/s - np.sqrt((x - m)/s*(x - m)/s+1))

    def __skewnormal_residuals(self, par, x, y):
        m, s, a, k = par
        return (a*np.exp( k*(x - m)/s - np.sqrt((x - m)/s*(x - m)/s+1))) - y

        

    def __skewnormal_cdf(self, x, m, s, a, k, range_start, range_end):
        values = []
        for value in x:
            integral = integrate.quad(lambda x: self.__skewnormal(x, m, s, a, k),range_start,value)[0]
            normalized = integral/integrate.quad(lambda x: self.__skewnormal(x, m, s, a, k),range_start,range_end)[0]
            values.append(normalized)
        return np.array(values)

    # ------------------------------ Error Functions ----------------------------- #
        
    def __sum_of_squared_error(self, theta, *data):
        warnings.filterwarnings("ignore") #do not print warnings by genetic algorithm
        data_x, data_y = data
        yDataPred = model(data_x, *theta)
        return np.sum((data_y - yDataPred) ** 2.0)

    def __coefficient_of_determination(self, theta, *data): #WIP
        warnings.filterwarnings("ignore")
        model, data_x, data_y = data
        yDataPred = model(data_x, *theta)
        return None

    def __MSPD(self, theta, *data):
        warnings.filterwarnings("ignore")
        model, data_x, data_y = data
        yDataPred = model(data_x, *theta)
        return 100 * math.sqrt( (1/(len(data_y)-4))*np.sum( ( (data_y-yDataPred)/data_y )**2 ) ) 

    
    # ------------------------------- Curve Fitting ------------------------------ #


    def __get_parameter_bounds_skewnormal(self,data_x,data_y):
        #TODO define clear param bound

        #parameters bounds
        parameterBounds = []
        parameterBounds.append([min(data_x)-((max(data_x)-min(data_x))) , max(data_x)+((max(data_x)-min(data_x)))+1]) # search bounds for m
        parameterBounds.append([0.1, stdev(data_x)*10+1]) # search bounds for s
        parameterBounds.append([0, max(data_y)*5+1]) # search bounds for a
        parameterBounds.append([-0.2, 0.8]) # search bounds for k

        return parameterBounds


    def __estimate_initial_parameters(self, model, data_x, data_y):

        parameter_estim = []
        parameter_estim.append(median(data_x)) # estimation for m
        parameter_estim.append(stdev(data_x)/2) # estimation for s
        parameter_estim.append(max(data_y)*2) # estimation for a
        parameter_estim.append(0) # estimation for k
        return np.array(parameter_estim)


    # def __estimate_initial_parameters(self, model, data_x, data_y):

    #     result = differential_evolution(self.__MSPD, self.__get_parameter_bounds_skewnormal(data_x, data_y), args =(model, data_x, data_y), seed=1)
    #     return result.x


    # -------------------------------- Scoring Fit ------------------------------- #

    #def __KSTest(self, parameters, data_x, data_y):
       # return stats.kstest(data_x, lambda x: self.__skewnormal_cdf(x, *parameters, min(data_x), max(data_x)))

    def __pearson_test(self, model, parameters, data_x, data_y): # !! NOT KS !! to be renamed 
        x = np.array(data_y).reshape((-1, 1))
        y = np.array([model(x, *parameters) for x in data_x])
        model = LinearRegression().fit(x, y)
        return model.score(x, y)
