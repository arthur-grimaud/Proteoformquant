
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
from scipy import integrate,stats
import numpy as np
import math
import warnings
from sklearn.linear_model import LinearRegression
from kneed import KneeLocator
from scipy.integrate import trapz, simps

class Envelope():

    def __init__(self, psms):


        self.corThreshold = 0.6 #threshold below which trimming is triggered

        self.xData = np.array([psm.spectrum.getRt() for psm in psms])
        self.yData = np.array([psm.getPrecIntensRatio() for psm in psms])
        
        self.estimatedParam = [None]
        self.fittedParam = [None]

        self.KsEstimated = None
        self.KsFitted = None

        self.splitScores = []
        self.split = None

        self.psmsOutliers = []


        self.estimatedParam, self.fittedParam, self.corEstimated, self.corFitted = self.fitSkewNormal(self.xData, self.yData)


        if self.corFitted < 0.6:
            nRemoved = []
            splitScores = []

            for n in range(0,len(self.xData)-5):
                if len(self.xData[:-n]) > 5:

                    #print("subset spectrum of proteoform")
                    
                    xDataT = self.xData[:-n]
                    yDataT = self.yData[:-n]

                    #print(xDataT)
                    #print(yDataT)

                    estimatedParam, fittedParam, corEstimated, corFitted = self.fitSkewNormal(xDataT, yDataT)
                    
                    nRemoved.append(n)        
                    splitScores.append(corFitted)

            #print(nRemoved)
            print(splitScores)
            if len(splitScores) >3:
                kn = KneeLocator(nRemoved, splitScores,S=3, curve='concave', direction='increasing',interp_method= "polynomial", polynomial_degree=2)
                cutoffVal = kn.knee
                if cutoffVal != None:


                    try:
                        self.xData = self.xData[:-cutoffVal]
                        self.yData = self.yData[:-cutoffVal]

                        self.estimatedParam, self.fittedParam, self.corEstimated, self.corFitted = self.fitSkewNormal(self.xData, self.yData)
                    
                    
                        self.psmsOutliers = psms[-cutoffVal:]
                    except(TypeError):
                        print("could not optimize fit by removing data points")
                        self.psmsOutliers = []

                #print([psm.spectrum.getRt() for psm in self.psmsOutliers])
                
            else:
                 self.envelopes = [] #delete envelope
        #else: #both steps could be performed ? 


        pass
    
    
    # ---------------------------------- Getters --------------------------------- #

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
 

    def fitSkewNormal(self, xData, yData):

        KsEstimated = None
        KsFitted = None
        scoreEstimated = 0
        scoreFitted = 0

        #Get Estimated Parameters
        geneticParameters = self.__generate_Initial_Parameters(xData, yData)
        scoreEstimated = self.__KSTest(geneticParameters, xData, yData)

        #Optimize model
        bounds = self.__getParameterBounds(xData, yData)
        bounds = tuple([ tuple([bounds[x][b] for x in range(len(bounds))])  for b in range(len(bounds[0]))]) #convert to curve fit bounds format
        try:
            fittedParameters = curve_fit(self.__skewNormal, xData, yData, geneticParameters,bounds=bounds)
            fittedParameters = fittedParameters[0]
            scoreFitted = self.__KSTest(fittedParameters, xData, yData)
        except(RuntimeError,TypeError):
            fittedParameters = [None]

        return geneticParameters, fittedParameters, scoreEstimated, scoreFitted

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
        parameterBounds.append( [min(xData)-500, max(xData)] ) # search bounds for m
        parameterBounds.append( [0.1, 10000 ] ) # search bounds for s
        parameterBounds.append( [0, max(yData)*2 ] ) # search bounds for a
        #parameterBounds.append( [0, 0] ) # search bounds for b
        parameterBounds.append( [-0.2, 0.9] ) # search bounds for k
        #print(parameterBounds)
        return parameterBounds

    def __generate_Initial_Parameters(self, xData, yData):

        # "seed" the numpy random number generator for repeatable results
        result = differential_evolution(self.__MSPD, self.__getParameterBounds( xData, yData), args =(xData, yData), seed=1)
        return result.x


    # -------------------------------- Scoring Fit ------------------------------- #

    #def __KSTest(self, parameters, xData, yData):
       # return stats.kstest(xData, lambda x: self.__skewNormalCdf(x, *parameters, min(xData), max(xData)))

    def __KSTest(self, parameters, xData, yData): # !! NOT KS !! to be renamed 
        x = np.array(yData).reshape((-1, 1))
        y = np.array([self.__skewNormal(x, *parameters) for x in xData])
        model = LinearRegression().fit(x, y)
        return model.score(x, y)
