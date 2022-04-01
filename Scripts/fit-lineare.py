#!/usr/bin/env python3

# 2022 - 04 - 01
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from lmfit import Model
import pandas as pd
import sys
from tkinter import Tk
import statsmodels.api as sm
import csvtotxt

def fit_lineare(file_name, dati, onlyB, error, method):
    # metodi
    metodi = {
        "1": "leastsq",
        "2": "least_squares",
        "3": "nelder",
        "4": "lbfgsb",
        "5": "basinhopping",
        "6": "ampgo",
        "7": "powell",
        "8": "cg",
        "9": "slsqp",
    }

    if method in metodi:
        method = metodi[method]
    else:
        print("Scegliere un metodo:")
        for keys, values in metodi.items():
            print("{}: {}".format(keys,values))
        sys.exit(1)

    # read x, y and y error from file using run number
    col_list = ["x {}".format(dati), "y {}".format(dati), "yerr {}".format(dati)]
    df = pd.read_csv(file_name, usecols=col_list, sep=",", decimal=".")

    # manipulate data into numpy array
    df = df.dropna()
    xdata = df['x {}'.format(dati)]
    ydata = df['y {}'.format(dati)]
    yerr = df['yerr {}'.format(dati)]

    #### fit ####
    if onlyB == "0":
        def func(x, A, B):
            return A + B*x
    elif onlyB == "1":
        def func(x, B):
            return B*x
    else:
        print("L'argomento dev'essere 0 oppure 1.")
        sys.exit(1)

    # set model's parameters
    fmodel = Model(func)

    if onlyB == "0" or onlyB == "1":
        fmodel.set_param_hint('B', value=1)
        if onlyB == "0":
            fmodel.set_param_hint('A', value=0)
    else:
        print("L'argomento dev'essere 0 oppure 1.")
        sys.exit(1)

    params = fmodel.make_params()

    # do fit, with either leastsq, least_squares or nelder
    if error == "0":
        result = fmodel.fit(ydata, params, x=xdata, method=method)
    elif error == "1":
        result = fmodel.fit(ydata, params, x=xdata, weights=1/yerr, method=method)
    else:
        print("L'argomento dev'essere 0 oppure 1.")
        sys.exit(1)

    print(result.fit_report())

    best_params = np.array(list(dict.values(result.best_values)))
    # dict.items() for array of tuples
    try:
        clipboard = Tk()
        clipboard.clipboard_clear()
        print("\nParametri:")
        for i, j, k in zip([item[0] for item in list(dict.items(params))], best_params, np.sqrt(np.diag(result.covar))):
            print("{}: {} +/- {}".format(i, j, k))
            clipboard.clipboard_append(repr(j) + " ")
    except:
        print("\nNo uncertainties.")

    clipboard.update()

    print("Matrice di covarianza:\n{}".format(result.covar))

    # Plot absolute residuals
    #print(result.residual)
    #plt.figure(figsize=(30,15))
    #result.plot_residuals()
    #plt.show()
    #plt.close()

    # Plot normalised residuals
    #plt.figure(figsize=(30,15))
    #plt.axhline(y=0)
    #plt.plot(np.array(xdata), result.residual, 'bo')
    #plt.show()
    #plt.close()

    # Normality tests for residuals
    #sm.qqplot(result.residual, fit=True, line='s')
    #plt.show()
    #print(stats.shapiro(result.residual))

    #### plot ####
    xdum=np.linspace(np.amin(xdata),np.amax(xdata), 1000)
    plt.figure(figsize=(30,15))
    #plt.plot(xdata, ydata, 'b*', label='data')
    plt.errorbar(xdata, ydata, yerr=yerr, fmt='bo', label='data')
    #plt.plot(xdata, result.best_fit, 'r-', label='best fit')
    plt.plot(xdum, func(xdum, *best_params), 'g--', label='best fit')
    plt.xlabel('xdata')
    plt.ylabel('ydata')
    plt.legend()

    try:
        plt.savefig('fit_func.png', bbox_inches='tight')
        plt.show()
    except:
        pass
    plt.close()

    print("Degrees of freedom = {}".format(len(xdata)-len(params)))
    print("p-value = {}".format(1-stats.chi2.cdf(result.chisqr,len(xdata)-len(params))))

    csvtotxt.csvtotxt(dati, file_name)

if __name__ == '__main__':
    if (len(sys.argv) != 6):
        print("Aspettati cinque argomenti: file name, run number, only B (0/1), incertezza y (0/1), metodo (1,2,3,...).")
        sys.exit(1)

    # run number, boolean for only B, boolean for error, fitting method
    fit_lineare(sys.argv[1:])