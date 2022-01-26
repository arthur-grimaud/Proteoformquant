from tkinter import Y
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
from scipy import integrate,stats
import numpy as np
import math
import warnings

class Envelope():

    def __init__(self, psms):

        self.xData = np.array([psm.spectrum.getRt() for psm in psms])
        self.yData = np.array([psm.getPrecIntensRatio() for psm in psms])
        
        self.estimatedParam = [None]
        self.fittedParam = [None]

        self.KsEstimated = None
        self.KsFitted = None


        self.estimatedParam, self.fittedParam = self.fitSkewNormal()


        #print(self.estimatedParam)
        pass
    
    
    #Getters

    def getEnvelopeSerie(self, xData, method = "best"):
        """for a list of x values return the estimated y values given the fitted function
        use estimated parameters if fitted parameters are not determined"""
        #print(self.__dict__)


        if self.fittedParam[0] != None and method in ["fitted","best"]:
            return self.__skewNormal(xData, *self.fittedParam), self.fittedParam
        elif (method in ["estimated","best"]):
            return self.__skewNormal(xData, *self.estimatedParam), self.estimatedParam
        else:
            return [None], [None]

    #Fitting


    def fitSkewNormal(self):
        
        geneticParameters = self.__generate_Initial_Parameters()
        #print(geneticParameters)

        self.KsEstimated = self.__KSTest(geneticParameters)
        #print(self.KsEstimated)
        try:
            fittedParameters = curve_fit(self.__skewNormal, self.xData, self.yData, geneticParameters)
            fittedParameters = fittedParameters[0]
            self.KsFitted = self.__KSTest(fittedParameters)
            #print("curve_fit succeeded")
        except(RuntimeError,TypeError):
            fittedParameters = [None]
            #print("curve_fit failed")

        return geneticParameters, fittedParameters



    def __skewNormal(self, x, m, s, a, k):
        return a*np.exp(k*((x - m))/s - np.sqrt(((x - m))/s*((x - m))/s+1)) 

    def __skewNormalCdf(self, x, m, s, a, k, range_start, range_end):
        values = []
        for value in x:
            integral = integrate.quad(lambda x: self.__skewNormal(x, m, s, a, k),range_start,value)[0]
            normalized = integral/integrate.quad(lambda x: self.__skewNormal(x, m, s, a, k),range_start,range_end)[0]
            values.append(normalized)
        return np.array(values)
        
    def __sumOfSquaredError(self, theta):
        warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
        val = self.__skewNormal(self.xData, *theta)
        return np.sum((self.yData - val) ** 2.0)


    def __generate_Initial_Parameters(self):

        # parameters bounds
        parameterBounds = []
        parameterBounds.append( [self.xData[np.argmin(self.yData)]-200, self.xData[np.argmin(self.yData)]+200] ) # search bounds for m
        parameterBounds.append( [(max(self.xData)-min(self.xData))/20, (max(self.xData)-min(self.xData))/5] ) # search bounds for s
        parameterBounds.append( [(max(self.yData)-min(self.yData)), (max(self.yData)-min(self.yData))*10] ) # search bounds for a
        #parameterBounds.append( [0, 0] ) # search bounds for b
        parameterBounds.append( [0, 0.999999]) # search bounds for k

        #print(parameterBounds)
        # "seed" the numpy random number generator for repeatable results
        result = differential_evolution(self.__sumOfSquaredError, parameterBounds, seed=1)
        return result.x
    
    def __KSTest(self, parameters):
        return stats.kstest(self.xData, lambda x: self.__skewNormalCdf(x, *parameters, min(self.xData), max(self.xData)))