# from scipy.optimize import curve_fit
# from scipy.optimize import differential_evolution
# import numpy as np
# import math


# class Envelope():

#     def __init__(self):


#         pass
    
    


#     def skewNormal(self, x, theta):

#         m = theta[1]
#         s = theta[2]
#         a = theta[3]
#         b = theta[4]
#         k = theta[5]

#         return a*math.exp(k*((x - m))/s - math.sqrt(((x - m))/s*((x - m))/s+1)) + b


#     def fitSkewNormal(self, psms):


#         xData = [ for  in ]
#         ydata = [ for  in ]

#         geneticParameters = self.generate_Initial_Parameters(xData, ydata)
#         param = curve_fit(self.skewNormal, xData, ydata)


    
#     def sumOfSquaredError(self, xData, yData, theta):
#         warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
#         val = self.skewNormal(xData, theta)
#         return np.sum((yData - val) ** 2.0)


#     def generate_Initial_Parameters(self, xData, yData):
#         # min and max used for bounds
#         maxX = max(xData)
#         minX = min(xData)
#         maxY = max(yData)
#         minY = min(yData)

#         parameterBounds = []
#         parameterBounds.append([x[which.max(y)], maxX]) # search bounds for m
#         parameterBounds.append([(max(x)-min(x))/10, maxX]) # search bounds for s
#         parameterBounds.append([0.0, maxY]) # search bounds for a
#         parameterBounds.append([0.0, maxY]) # search bounds for a
#         parameterBounds.append([0.0, maxY]) # search bounds for a

#         # "seed" the numpy random number generator for repeatable results
#         result = differential_evolution(self.sumOfSquaredError, parameterBounds, seed=3)
#         return result.x