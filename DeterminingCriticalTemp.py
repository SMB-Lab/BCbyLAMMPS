# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 08:56:40 2024

@author: Dani

This script calculates the critical temperature 
using Equation 6 from Dignon et al, 2018.
Takes the same input file as the python phase
fitting script from Brady et al, 2017:
*.txt file that is tab delimited with no header
Temperature in Kelvin
Low Density (mg/mL)
Low Density Std
High Density (mg/mL)
High Density Std
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

InputFile = 'file.txt'

os.chdir(r'FILE LOCATION')
#%%
data = np.loadtxt(InputFile, delimiter='\t')
data_fit = data
#print(data)
#print(data.shape[0])

k=0
for i in range(0,data.shape[0]):
    if data[i,1] == 0.01:
        data_fit = np.delete(data_fit, k, axis=0)
        k = k - 1
    k = k + 1   
print(data)
print(data_fit)
#%%
# Define the function you want to fit (example: a linear function with two parameters)
def fun(Temp, A, T_c):
    return A * (T_c - Temp)**(0.325)

den_H = data_fit[:,3]
den_L = data_fit[:,1]
Temp = data_fit[:,0]

den_diff = den_H - den_L

# Perform the curve fit
params, covariance = curve_fit(fun, Temp, den_diff, p0=[5, 400])

# Extract the fitted parameters
a, b = params
print(a, b)

Temp_range = np.arange(300, b, 1)
# Create a fitted line using the fitted parameters
fitted_line = fun(Temp_range, a, b)

# Plot the original data and the fitted line
plt.scatter(den_diff, Temp, label='Original Data')
plt.plot(fitted_line, Temp_range, label=f'Fitted Line (A={a:.2f}, T_c={b:.2f})')
plt.ylabel('Temperature (K)')
plt.xlabel('Density Difference (mg/mL)')
plt.title(InputFile)
plt.legend()
plt.show()
