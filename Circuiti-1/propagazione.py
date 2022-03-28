#!/usr/bin/env python3

# 2022 - 03 - 16
import numpy as np
import sympy as smp
import pandas as pd

####### Dati utente #######
dati = 1 
variables = 2

# read measurements from file using run number
col_list = ["x{} {}".format(i+1, dati) for i in range(variables)]
df = pd.read_csv("dstd.csv", usecols=col_list, sep=";", decimal=".")

# manipulate data into numpy array
df = df.dropna()
misure = np.array([df['x{} {}'.format(i+1, dati)] for i in range(variables)])

# calculate the covariance matrix, mean and std dev of each set of measurements
cov = np.cov(misure)
mean = np.array([np.mean(i) for i in misure])
dstd_mis = np.array([np.std(i, ddof = 1) for i in misure])

# create variables and function
x = [smp.symbols('x%d' % i) for i in range(variables)]
####### Funzione utente #######
f = (1/x[0]-1/x[1])**(-1)

# calculate standard deviation
subs = [(x[i], mean[i]) for i in range(variables)]
coeff = np.array([smp.diff(f,i).subs(subs) for i in x])
mean_f = f.subs(subs)
var = float(coeff.dot(cov).dot(coeff.T))
dstd = np.sqrt(var)

# print
print("misure: {}\nmean mis: {}\ndstd mis: {}\nsubs: {}\ncoeff: {}\nmean func: {}\nvar: {}\ndstd: {}".format(misure, mean, dstd_mis, subs, coeff, mean_f, var, dstd))