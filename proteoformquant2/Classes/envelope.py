
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
from scipy import integrate,stats
import numpy as np
import math
import warnings
from sklearn.linear_model import LinearRegression
from kneed import KneeLocator
from scipy.integrate import trapz, simps
from statistics import mean

class Envelope():

    def __init__(self, psms, scoreThreshold):
        

        self.psms = psms

        self.scoreThreshold = scoreThreshold #threshold below which trimming is triggered

        self.xData = np.array([psm.spectrum.getRt() for psm in psms])
        self.yData = np.array([psm.getPrecIntensRatio() for psm in psms])
        
        self.estimatedParam = [None]
        self.fittedParam = [None]

        self.scoreEstimated = None
        self.scoreFitted = None

        self.splitScores = []
        self.split = None

        self.psmsOutliers = []


        self.estimatedParam, self.fittedParam, self.scoreEstimated, self.scoreFitted = self.fitSkewNormal(self.xData, self.yData)


        if self.scoreFitted < self.scoreThreshold:
            self.excludeOutliersMeanMethod()
            
        #else: #both steps could be performed ? 


        pass
    
    
    # ---------------------------------- Getters --------------------------------- #

    def getEnvelopeSerie(self, xData, method = "best"):
        """for a list of x values return the estimated y values given the fitted function
        use estimated parameters if fitted parameters are not determined or is method = best"""


        if self.fittedParam[0] != None and method in ["fitted","best"]:
            return self.__skewNormal(xData, *self.fittedParam), self.fittedParam
        elif (method in ["estimated","best"]):
            return self.__skewNormal(xData, *self.estimatedParam), self.estimatedParam
        else:
            return [None], [None]

    def getY(self, x, method = "best"):
        """for x return the estimated y value given the fitted function
        uses estimated parameters if fitted parameters are not determined"""
        if self.fittedParam[0] != None and method in ["fitted","best"]:
            return self.__skewNormal(x, *self.fittedParam)
        elif (method in ["estimated","best"]):
            return self.__skewNormal(x, *self.estimatedParam)
        else:
            return None

    def getEnvelopeMean(self, method):
        if method == "fitted":
            m, s, a, k =  estimatedParam[0], estimatedParam[1], estimatedParam[2], estimatedParam[3]
        elif method == "estimated":
            m, s, a, k =  parametersFitted[0], parametersFitted[1], parametersFitted[2], parametersFitted[3]
        else: 
            print("Specify method (fitted or estimated)")

        return m + (math.sqrt(2/math.pi)) * ((s*a)/math.sqrt(1+(a**2)) ) 

    def getEnvelopeStd(self, method ): #might be incorrect
        if method == "fitted":
            return estimatedParam[1]
        elif method == "estimated":
            return parametersFitted[1]
        else:
            print("Specify method (fitted or estimated) ")

    def getAUC(self):


        try:
            return integrate.quad(lambda x: self.__skewNormal(x, *self.fittedParam),0,10000)[0]
        except(TypeError):
            return 0
        

    # ---------------------------------- Fitting --------------------------------- #

    def fitSkewNormal(self, xData, yData, startPrevFit = False):
        """ startPrevFit  """


        if startPrevFit == True and self.fittedParam != [None]:

            estimatedParam = self.fittedParam

        else:

            scoreEstimated = 0
            scoreFitted = 0
            #Get Estimated Parameters
            geneticParameters = self.__generate_Initial_Parameters(xData, yData)
            scoreEstimated = self.__PearsonTest(geneticParameters, xData, yData)
            #estimated param from genetic algo
            estimatedParam = geneticParameters
        

        #Optimize model
        bounds = self.__getParameterBounds(xData, yData) #Get param bounds
        bounds = tuple([ tuple([bounds[x][b] for x in range(len(bounds))])  for b in range(len(bounds[0]))]) #convert to "curve_fit" bounds format

        try:
            fittedParameters = curve_fit(self.__skewNormal, xData, yData, estimatedParam, bounds=bounds)
            fittedParameters = fittedParameters[0]
            scoreFitted = self.__PearsonTest(fittedParameters, xData, yData, )
        except(RuntimeError,TypeError):
            fittedParameters = [None]

        return geneticParameters, fittedParameters, scoreEstimated, scoreFitted

    # ---------------------------------- Outlier detection --------------------------------- #

    def excludeOutliersRightMethod(self):
        """Try to improve the fit of a curve by removing data from the highest RT to lowest"""

        print("Start subset")
        print("psms")

        scores=[self.scoreFitted] #list of scores for each subset
        indexes=[0]
        psmsSubsets = [self.psms]

        if len(self.xData) >= 7: #limit prob 7
            for n in range(1,len(self.xData)-5):
                    #subset of spectra
                    psmsSubset = self.psms[:-n]
                    xDataT = np.array([psm.spectrum.getRt() for psm in psmsSubset])
                    yDataT = np.array([psm.getPrecIntensRatio() for psm in psmsSubset])
                    #refit the curve using subset
                    estimatedParam, fittedParam, scoreEstimated, scoreFitted = self.fitSkewNormal(xDataT, yDataT)
                    #store score of subset
                    scores.append(scoreFitted)
                    indexes.append(n)
                    psmsSubsets.append(psmsSubset)

        
            print(scores)

            kn = KneeLocator(indexes, scores, S=2, curve='concave', direction='increasing',interp_method= "polynomial", polynomial_degree=2)
            index = kn.knee

            print(index)

            if index != None and index != 0:
                
                print(psmsSubsets[index])
                try:
                    self.xData = np.array([psm.spectrum.getRt() for psm in psmsSubsets[index]])
                    self.yData = np.array([psm.getPrecIntensRatio() for psm in psmsSubsets[index]])
                    self.estimatedParam, self.fittedParam, self.scoreEstimated, self.scoreFitted = self.fitSkewNormal(self.xData, self.yData)
                    self.psmsOutliers = [psm for psm in self.psms if psm not in psmsSubsets[index]] #TODO check whether outlier needs to be removed for proteoform.linkedPSMSs
                except(TypeError):
                    print("could not optimize fit by removing data points")
                    self.psmsOutliers = []
            else:   
                print("could not optimize fit by removing data points")
                pass 

    def excludeOutliersMeanMethod(self):
        "Try imrove the fit of the curve by iteratively removing datapoints the furthest from the RT mean"

        scores=[self.scoreFitted] #list of scores for each subset
        indexes=[0]
        psmsSubsets = [self.psms]
        print(".........................")
        print(self.psms[0].proteoform.getModificationBrno())
        print(self.psms)
        if len(self.xData) >= 8: #limit prob 7
            for n in range(1,len(self.xData)-5):
                psmsSubset = psmsSubsets[n-1]
                #get index of furthest point form mean
                xDataT = [psm.spectrum.getRt() for psm in psmsSubset]
                xMean = mean(xDataT)
                outPsmIndex = xDataT.index(max(xDataT, key=lambda x:abs(x-xMean)))
                #exclude psm
                psmsSubset.pop(outPsmIndex)
                #subset 
                xDataT = np.array([psm.spectrum.getRt() for psm in psmsSubset])
                yDataT = np.array([psm.getPrecIntensRatio() for psm in psmsSubset])
                #refit the curve to subset
                estimatedParam, fittedParam, scoreEstimated, scoreFitted = self.fitSkewNormal(xDataT, yDataT)
                #store score and subset
                scores.append(scoreFitted)
                indexes.append(n)
                psmsSubsets.append(psmsSubset)

            print(scores)
            print(indexes)
            print(psmsSubsets)

            kn = KneeLocator(indexes, scores, S=2, curve='concave', direction='increasing',interp_method= "polynomial", polynomial_degree=2)
            index = kn.knee
            print(index)

            if index != None and index != 0:
                
                print(psmsSubsets[index])
                try:
                    self.xData = np.array([psm.spectrum.getRt() for psm in psmsSubsets[index]])
                    self.yData = np.array([psm.getPrecIntensRatio() for psm in psmsSubsets[index]])
                    self.estimatedParam, self.fittedParam, self.scoreEstimated, self.scoreFitted = self.fitSkewNormal(self.xData, self.yData)
                    self.psmsOutliers = [psm for psm in self.psms if psm not in psmsSubsets[index]] #TODO check whether outlier needs to be removed for proteoform.linkedPSMSs
                except(TypeError):
                    print("could not optimize fit by removing data points")
                    self.psmsOutliers = []
            else:   
                print("could not optimize fit by removing data points")
                pass 

        

    def excludeOutlierNonSignificant():
        "Exclude datapoint on the tails of the envellope that would represent less than X%  of the area under the curve"

    # ----------------------------- Model's function ----------------------------- #

    def __skewNormal(self, x, m, s, a, k):
        return a*np.exp( k*(x - m)/s - np.sqrt((x - m)/s*(x - m)/s+1)) 

    def __skewNormalCdf(self, x, m, s, a, k, range_start, range_end):
        values = []
        for value in x:
            integral = integrate.quad(lambda x: self.__skewNormal(x, m, s, a, k),range_start,value)[0]
            normalized = integral/integrate.quad(lambda x: self.__skewNormal(x, m, s, a, k),range_start,range_end)[0]
            values.append(normalized)
        return np.array(values)

    # ------------------------------ Error Functions ----------------------------- #
        
    def __sumOfSquaredError(self, theta, *data):
        warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
        xData, yData = data
        yDataPred = self.__skewNormal(xData, *theta)
        return np.sum((yData - yDataPred) ** 2.0)

    def __coefficientOfDetermination(self, theta, *data): #WIP
        warnings.filterwarnings("ignore")
        xData, yData = data
        yDataPred = self.__skewNormal(xData, *theta)
        return None

    def __MSPD(self, theta, *data):
        warnings.filterwarnings("ignore")
        xData, yData = data
        yDataPred = self.__skewNormal(xData, *theta)
        return 100 * math.sqrt( (1/(len(yData)-4))*np.sum( ( (yData-yDataPred)/yData )**2 ) ) 

    
    # ------------------------------- Curve Fitting ------------------------------ #


    def __getParameterBounds(self,xData,yData):
        #TODO define clear param bound

        # parameters bounds
        parameterBounds = []
        parameterBounds.append([ min(xData)-((max(xData)-min(xData))/2) , max(xData)+((max(xData)-min(xData))/2) ]) # search bounds for m
        parameterBounds.append([0.1, 50 ]) # search bounds for s
        parameterBounds.append( [0, max(yData)*3 ] ) # search bounds for a
        parameterBounds.append( [-0.2, 0.8] ) # search bounds for k

        return parameterBounds

    def __generate_Initial_Parameters(self, xData, yData):

        # "seed" the numpy random number generator for repeatable results
        result = differential_evolution(self.__MSPD, self.__getParameterBounds( xData, yData), args =(xData, yData), seed=1)
        return result.x


    # -------------------------------- Scoring Fit ------------------------------- #

    #def __KSTest(self, parameters, xData, yData):
       # return stats.kstest(xData, lambda x: self.__skewNormalCdf(x, *parameters, min(xData), max(xData)))

    def __PearsonTest(self, parameters, xData, yData): # !! NOT KS !! to be renamed 
        x = np.array(yData).reshape((-1, 1))
        y = np.array([self.__skewNormal(x, *parameters) for x in xData])
        model = LinearRegression().fit(x, y)
        return model.score(x, y)
