"""
Propagate
=========
Libreria per la propagazione degli errori utilizzando le espressioni simboliche di SymPy.
"""
import numpy as np
import pandas as pd
import sympy as smp

def LoadData(file, dati, variables) -> np.ndarray:
    if file is None:
        print("Defaulting to \"./dstd.csv\" file.")
        file = "./dstd.csv"
    # read measurements from file using run number
    col_list = [f"x{i+1} {dati}" for i in range(variables)]
    df = pd.read_csv(file, usecols=col_list, sep=",", decimal=".")

    # manipulate data into numpy array
    df = df.dropna()
    return np.array([df[f'x{i+1} {dati}'] for i in range(variables)])

def GenerateXs(variables) -> list:
    return [smp.symbols('x%d' % i) for i in range(variables)]

def PropagateFile(file, dati, variables = 2, func=None, x=None) -> None:
    misure = LoadData(file, dati, variables)
    # calculate the covariance matrix, mean and std dev of each set of measurements
    cov = np.cov(misure)
    mean = np.array([np.mean(i) for i in misure])
    dstd_mis = np.array([np.std(i, ddof = 1) for i in misure])

    if func is None:
        # create variables and function
        x = GenerateXs(variables)
        ####### Funzione utente #######
        func = x[0]-x[1]

    # calculate standard deviation
    subs = [(x[i], mean[i]) for i in range(variables)]
    coeff = np.array([smp.diff(func,i).subs(subs) for i in x])
    mean_f = func.subs(subs)
    var = float(coeff.dot(cov).dot(coeff.T))
    dstd = np.sqrt(var)

    print(f"mean measures: {mean}\nstd dev: {dstd_mis}\nderivatives: {coeff}\nmean function: {mean_f}\nstd dev: {dstd}")

def Propagate(mean, cov, variables = 2, func=None, x=None) -> None:
    if func is None:
        # create variables and function
        x = GenerateXs(variables)
        ####### Funzione utente #######
        func = x[0]-x[1]

    # calculate standard deviation
    subs = [(x[i], mean[i]) for i in range(variables)]
    coeff = np.array([smp.diff(func,i).subs(subs) for i in x])
    mean_f = func.subs(subs)
    var = float(coeff.dot(cov).dot(coeff.T))
    dstd = np.sqrt(var)
    #print(f"mean measures: {mean}\nstd dev: {np.sqrt(np.diag(cov))}\nderivatives: {coeff}\nmean function: {mean_f}\nstd dev: {dstd}")
    return [mean_f, dstd]

def CovarX (var, cov) -> np.ndarray:
    if var == cov.shape[0]:
        return np.diag(cov)
    X = np.zeros((var,var))
    X[np.triu_indices(X.shape[0], k=0)] = cov
    X=X+X.T-np.diag(np.diag(X))
    return X