'''
Created on Aug 30, 2018

@author: Vincent
'''
from numpy import average
from statistics import stdev


#import alpha_vantage

#from alpha_vantage.timeseries import TimeSeries
#ts = TimeSeries(key='UQV2D45VFXJGX00J', output_format='pandas')
#data, meta_data = ts.get_intraday(symbol = 'GE', interval = '1min')

no_vol = [17.93, 14.47, 14.43, 15.38, 15.34, 14.02, 14.87, 15.49, 15.51, 15.87, 13.69, 15.92, 15.15,
          14.89, 14.305, 17.075, 14.5, 15.56, 15.905, 16.25, 15.29, 14.73, 16.19, 13.47, 14.56, 15.36, 
          13.5, 15.48, 14.86, 15.46]
vol = [16.69, 16.57, 14.8, 15.25, 17.03, 14.79, 16.82, 15.1, 16.72, 15.64, 15.72, 16.21, 15.13, 16.0, 15.15,
       15.94, 15.62, 17.83, 14.4, 15.52, 16.7, 13.56, 17.32, 17.34, 15.96, 16.05, 15.58, 15.63, 16.07, 16.47]

ten_gens = [30.80, 29.62, 32.46, 33.20, 32.98, 33.47, 34.59, 35.55, 33.77, 35.08, 32.15, 34.25, 33.99, 32.36, 31.90, 
            33.96, 30.64, 34.37, 35.27, 31.33, 33.44, 32.90, 30.52, 22.29, 33.20, 27.38, 32.87, 34.44, 34.75,33.99]
one_gen = [29.48, 28.17, 28.60, 34.96, 28.57, 32.03, 27.28, 30.28, 32.66, 27.76, 27.89, 30.68, 29.28, 30.83, 32.31, 
           33.45, 33.51, 37.42, 30.78, 27.75, 30.70, 31.61, 29.97, 30.81, 29.97, 32.95, 30.24, 33.51, 29.49, 31.53]
three_gen = [49.58, 55.95, 48.90, 50.65, 47.25, 53.43, 50.48, 56.84, 49.62, 48.46, 52.57, ]
print("no vol avg: " + str(round(average(no_vol),2)) + 
      "\n std dev: " + str(round(stdev(no_vol),2))) 
print("vol avg: " + str(round(average(vol),2)) + 
      "\n std dev: " + str(round(stdev(vol),2)))

print("ten generations avg: " + str(round(average(ten_gens),2)) + 
      "\n std dev: " + str(round(stdev(ten_gens),2)))
print("one generation avg: " + str(round(average(one_gen),2)) + 
      "\n std dev: " + str(round(stdev(one_gen),2)))
