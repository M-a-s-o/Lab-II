#!/usr/bin/env python3

# 2022 - 03 - 16
import numpy as np
import sympy as smp
import pandas as pd

####### Dati utente #######
dati = 1 
variables = 2

# calculate the covariance matrix, mean and std dev of each set of measurements
cov = np.array([[653.9833567828732**2, 0], [0, 802.4612219083949**2]])
mean = np.array([11375.12901881004, 14355.853446031802])

# create variables and function
x = [smp.symbols('x%d' % i) for i in range(variables)]
####### Funzione utente #######
f = x[0]/x[1]

# calculate standard deviation
subs = [(x[i], mean[i]) for i in range(variables)]
coeff = np.array([smp.diff(f,i).subs(subs) for i in x])
mean_f = f.subs(subs)
var = float(coeff.dot(cov).dot(coeff.T))
dstd = np.sqrt(var)

# print
print("mean mis: {}\nsubs: {}\ncoeff: {}\nmean func: {}\nvar: {}\ndstd: {}".format(mean, subs, coeff, mean_f, var, dstd))