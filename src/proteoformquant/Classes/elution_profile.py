from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
from scipy import integrate, stats
import numpy as np
import math
import warnings
from sklearn.linear_model import LinearRegression
from scipy.integrate import trapz, simps
from statistics import mean, median, stdev
from scipy.optimize import least_squares
from sympy import subsets, true
import sympy as sy
from scipy.optimize import fsolve as fsolve
from scipy.integrate import quad
import mpmath as mp
from scipy.stats import norm
from scipy.special import owens_t
from scipy.stats import spearmanr
import warnings
from numba import jit


class ElutionProfile:
    def __init__(self):

        # Parameters of the function
        self.param_estimated = [None]
        self.param_fitted = [None]

        # Goodness of fit
        self.score_estimated = None
        self.score_fitted = None

        # List of PSMs
        self.psms_included = []  # all psms used to compute the elution profile
        self.psms_outliers = []  # psms considered as outliers in the modelling of the elution profile

        # Score threshold
        self.score_threshold = None

        pass

    # ---------------------------------- Getters --------------------------------- #

    def is_parameters_fitted(self):
        # return true if curve_fit results is available in the instance
        if self.param_fitted[0] == None:
            return False
        else:
            return True

    def is_modeled(self):
        # return true if curve_fit results is available in the instance
        if self.param_fitted[0] == None:
            return False
        else:
            return True

    def get_y_serie(self, x_values, method="best"):
        """
        for a list of x values return the estimated y values given the fitted function
        use estimated parameters if fitted parameters are not determined or is method = best
        """
        if self.param_fitted[0] != None and method in ["fitted", "best"]:
            return self.skewnormal(x_values, *self.param_fitted), self.param_fitted
        elif method in ["estimated", "best"]:
            return self.skewnormal(x_values, *self.param_estimated), self.param_estimated
        else:
            return [None], [None]

    def get_y(self, x):
        """
        for x return the estimated y value given the fitted function
        uses estimated parameters if fitted parameters are not determined
        """
        if self.param_fitted[0] != None:  # and method in ["fitted", "best"]:
            return self.skewnormal(x, *self.param_fitted)
        # elif method in ["estimated", "best"]:
        #     return self.skewnormal(x, *self.param_estimated)
        else:
            return None

    def get_score(self):
        if self.is_parameters_fitted():
            return self.score_fitted
        else:
            return 0

    def get_parameters_fitted(self):
        # returns parameters of curve_fit
        return self.param_fitted

    def get_elution_profile_param_m(self, method):
        if method == "estimated":
            return self.param_estimated[0]
        elif method == "fitted":
            return self.param_fitted[0]
        else:
            print("Specify method (fitted or estimated)")

    def get_elution_profile_std(self, method):
        if method == "estimated":
            return self.param_estimated[1]
        elif method == "fitted":
            return self.param_fitted[1]
        else:
            print("Specify method (fitted or estimated) ")

    def get_elution_profile_param_k(self, method):
        if method == "estimated":
            return self.param_estimated[1]
        elif method == "fitted":
            return self.param_fitted[1]
        else:
            print("Specify method (fitted or estimated) ")

    def get_auc(self):

        try:
            return self.get_parameters_fitted()[2]
        except (TypeError):

            return 0

    def bounds_area_equation(self, x1x2, m, s, a, k, area_percent):
        # print("params of bounds area equation")
        # print(x1x2, m, s, a, k)

        x1, x2 = x1x2

        # percentage area
        total_auc = a
        bound_auc = self.skewnorm_cdf(x2, m, s, a, k) - self.skewnorm_cdf(x1, m, s, a, k)
        area_95 = total_auc * area_percent - bound_auc

        # x1 x2 at same height
        fx1_fx2_equal = self.skewnormal(x1, *[m, s, a, k]) - self.skewnormal(x2, *[m, s, a, k])

        return (area_95, fx1_fx2_equal)

    def get_bounds_area(self, area_percent=0.95):

        m, s, a, k = self.param_fitted

        # print(m, s, a, k)

        # negate warning and count number of warnings
        with warnings.catch_warnings(record=True) as w:
            res = fsolve(
                self.bounds_area_equation,
                [m, m],
                args=(m, s, a, k, area_percent),
                factor=0.2,  # TODO factor probably increase time maybe implement consecutive search with lower factor only if minima isn't reached
            )
        

        # print("m, s, a, k ", m, s, a, k)
        # print("res ", res)
        return res

    def get_total_auc(self):  # TODO hard-coded
        return self.get_auc(-10000, 20000)

    def get_x_at_max_y(self):
        """Returns the x value for the peak of the elution profile COULD EXIST A MATHEMATICAL SOLUTION"""
        m = self.get_elution_profile_param_m("fitted")
        k = self.get_elution_profile_param_k("fitted")

        y = self.get_y(m)

        search = True
        x = m

        if k > 0:
            while search:
                if y < self.get_y(x + 1):
                    x = x + 1
                    y = self.get_y(x)
                else:
                    search = False

        elif k < 0:
            while search:
                if y < self.get_y(x - 1):
                    x = x - 1
                    y = self.get_y(x)
                else:
                    search = False

        return x

    def get_y_range_central_auc(self, auc_percent=99):
        """get the range on the y axis were the auc representing auc_percent of the total auc"""
        print("WIP")
        pass

    def get_x_range_for_y(self, Y):
        """Return the estimation of x values for y = constante"""

        m = self.param_fitted[0]
        s = self.param_fitted[1]
        a = self.param_fitted[2]
        k = self.param_fitted[3]

        results_min, infodict, ier, mesg = fsolve(
            lambda x: self.skewnormal(x, *self.param_fitted) - Y, m - (s * (1 + np.abs(k))), full_output=True
        )

        if ier != 1:  # if no solution if found
            return (0, 0)

        results_max, infodict, ier, mesg = fsolve(
            lambda x: self.skewnormal(x, *self.param_fitted) - Y, m + (s * (1 + np.abs(k))), full_output=True
        )

        if ier != 1:  # if no solution if found
            return (0, 0)

        return (results_min[0], results_max[0])

    # ---------------------------------- Fitting --------------------------------- #

    def model_elution_profile(self, psms, score_threshold):

        self.psms_included = psms
        self.score_threshold = score_threshold

        self.data_x = np.array([psm.spectrum.get_rt() for psm in psms])
        self.data_y = np.array([psm.get_prec_intens_ratio() for psm in psms])

        # Fit model without outliers removal
        (
            self.param_estimated,
            self.param_fitted,
            self.score_estimated,
            self.score_fitted,
        ) = self.fit_skew_normal(self.data_x, self.data_y)

    def fit_skew_normal(self, data_x, data_y, param_init=[None], param_bounds=[None]):
        """startPrevFit"""

        # Transform to numpy array (might not be needed)
        data_x = np.array(data_x)
        data_y = np.array(data_y)

        # Determine parameters bounds
        if param_bounds[0] == None:  # If no parameters bounds are suplied, estimate them
            param_bounds = self.__get_parameter_bounds_skewnormal(data_x, data_y)

        # print(param_bounds)

        # Determine initial parameters
        if param_init[0] == None:  # If no initial parameters are suplied, estimate them
            param_init = self.__estimate_initial_parameters(self.skewnormal, data_x, data_y)

        # Save initial parameters as "estimated"
        param_estimated = param_init
        # goodness of the fit with estimated parameters
        score_estimated = self.__pearson_test(self.skewnormal, param_estimated, data_x, data_y)

        try:
            param_bounds_curve_fit = tuple(
                [
                    tuple([param_bounds[x][b] for x in range(len(param_bounds))])
                    for b in range(len(param_bounds[0]))
                ]
            )  # convert bounds for "curve_fit"

            # print("data for fitting")
            # print((data_x, data_y))

            param_fitted = least_squares(
                self.skewnormal_residuals,
                param_estimated,
                bounds=param_bounds_curve_fit,
                loss="linear",
                f_scale=1,
                args=(data_x, data_y),
            )

            param_fitted = param_fitted.x

            # #Optimize curve fit
            # param_bounds_curve_fit = tuple([ tuple([param_bounds[x][b] for x in range(len(param_bounds))]) for b in range(len(param_bounds[0]))]) #convert bounds for "curve_fit"
            # param_fitted = curve_fit(self.__skewnormal, data_x, data_y, param_estimated, bounds=param_bounds_curve_fit)[0]
            # #goodness of the fit with fitted parameters

            score_fitted = self.__pearson_test(self.skewnormal, param_fitted, data_x, data_y)
            # print(score_fitted)

        except (RuntimeError, ValueError):
            param_fitted = [None]
            score_fitted = 0
            # print("Curve_fit failed")

        return param_estimated, param_fitted, score_estimated, score_fitted

    # ----------------------------- Model's function ----------------------------- #

    def skewnormal(self, x, m, s, a, k):
        u = (x - m) / s
        return a * ((2 / s) * norm.pdf(u) * norm.cdf(k * u))

    def skewnorm_cdf(self, x, m, s, a, k):
        u = (x - m) / s
        return a * (norm.cdf(u) - 2 * owens_t(u, k))

    def skewnormal_residuals(self, par, x, y):
        m, s, a, k = par
        return self.skewnormal(x, m, s, a, k) - y

    # ------------------------------ Error Functions ----------------------------- #

    def __sum_of_squared_error(self, theta, *data):
        warnings.filterwarnings("ignore")  # do not print warnings by genetic algorithm
        data_x, data_y = data
        yDataPred = model(data_x, *theta)
        return np.sum((data_y - yDataPred) ** 2.0)

    def __coefficient_of_determination(self, theta, *data):  # WIP
        warnings.filterwarnings("ignore")
        model, data_x, data_y = data
        yDataPred = model(data_x, *theta)
        return None

    def __MSPD(self, theta, *data):
        warnings.filterwarnings("ignore")
        model, data_x, data_y = data
        yDataPred = model(data_x, *theta)
        return 100 * math.sqrt((1 / (len(data_y) - 4)) * np.sum(((data_y - yDataPred) / data_y) ** 2))

    # ------------------------------- Curve Fitting ------------------------------ #

    def __get_parameter_bounds_skewnormal(self, data_x, data_y):
        # TODO define clear param bound

        # parameters bounds
        parameterBounds = []
        parameterBounds.append(
            [min(data_x) - ((max(data_x) - min(data_x))), max(data_x) + ((max(data_x) - min(data_x))) + 1]
        )  # search bounds for m
        parameterBounds.append([0.1, stdev(data_x) * 1500 + 1])  # search bounds for s
        parameterBounds.append([0, 15000])  # search bounds for a
        parameterBounds.append([-7.5, 7.5])  # search bounds for k

        return parameterBounds

    def __estimate_initial_parameters(self, model, data_x, data_y):

        parameter_estim = []
        parameter_estim.append(median(data_x))  # estimation for m
        parameter_estim.append(stdev(data_x))  # estimation for s
        parameter_estim.append(max(data_y) * 10)  # estimation for a
        parameter_estim.append(0)  # estimation for k
        return np.array(parameter_estim)

    # -------------------------------- Scoring Fit ------------------------------- #

    def __pearson_test(self, model, parameters, data_x, data_y):  
        x = np.array(data_y).reshape((-1, 1))
        y = np.array([model(x, *parameters) for x in data_x])

        
        
        try:
            cor_results = spearmanr(x, y)
        except (Warning):
            print("Issue in spearmanR scoring, returning score of 0")
            return 0

        if cor_results.correlation < -1 or cor_results.correlation > 1:
            print("Issue in spearmanR scoring, returning score of 0")
            return 0

        if math.isnan(cor_results.correlation) or cor_results.correlation in ["nan", np.nan, math.nan]:
            return 0

        return cor_results.correlation

    def scoring_return(self, model, parameters, data_x, data_y):  # for testing purpose return pred vs expect
        x = np.array(data_y)
        # x = np.log2(x, out=np.zeros_like(x), where=(x != 0))
        y = np.array([model(x, *parameters) for x in data_x])
        # y = np.log2(y, out=np.zeros_like(y), where=(y != 0))
        # print(x)
        # print(y)
        return x, y
