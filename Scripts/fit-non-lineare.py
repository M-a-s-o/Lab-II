#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from lmfit import Model
import pandas as pd
import sys
from tkinter import Tk
import statsmodels.api as sm
import csvtotxt

def fit_non_lineare(file_name, dati, interpol, error, method):
    # Legge x, y e yerr da file usando numero dati
    col_list = [f"x {dati}", f"y {dati}", f"yerr {dati}"]
    df = pd.read_csv(file_name, usecols=col_list, sep=",", decimal=".")

    # Separa il dataframe nei rispettivi array
    df = df.dropna()
    xdata = df[f'x {dati}']
    ydata = df[f'y {dati}']
    yerr = df[f'yerr {dati}']

    #### Interpolazioni ####
    if interpol=="1":
        ## prima interpolazione ## 
        def func(xval, I0, fattore, A):
            return I0*(np.exp(fattore*xval)-1)+A
        pars = [('I0', 1), ('fattore', 38.6), ('A', 0)]
    elif interpol == "3":
        ## RC
        def func(xval, V_g, tau, A):
            periodo = 1./.8
            return V_g*(1-2/(1+np.exp(-periodo/(2*tau)))*np.exp(-xval/tau))+A
        pars = [('V_g', 9.875198121607973), ('tau', 0.11148566083571734), ('A', 9.826661090020343)] # Nelder-Mead
    elif interpol == "4":
        ## RL
        def func(xval, V_g, tau, A):
            periodo = 1/100
            return V_g*(2*np.exp(-xval/tau)/(1+np.exp(-periodo/(2*tau))))+A
        pars = [('V_g', 1.8897933832058271), ('tau', 0.07711593596900626), ('A', 0.11063333194641872)] # Nelder-Mead
    elif interpol == "5":
        ## RL reale
        def func(xval, V_g, R_L, tau, A):
            periodo = 1/100
            R = 992
            return V_g*(2*np.exp(-xval/tau)/(1+np.exp(-periodo/(2*tau)))+R_L/(R_L+R)*(1-2*np.exp(-xval/tau)/(1+np.exp(-periodo/(2*tau)))))+A
        pars = [('V_g', 1.0323245746163687), ('R_L', 57.769650453274096), ('tau', 7.711618152518606e-05), ('A', 0.05382326137423557)] # Nelder-Mead
    elif interpol == "6":
        ## RLC sotto
        def func(xval, I0, gamma, omega, A, B):
            R = 10.
            return R*I0*np.exp(-gamma*xval)*np.sin(omega*xval+B)+A
        pars = [('I0', 1.842206038819525), ('gamma', 6.25077381318831), ('omega', 310.3087890994233), ('A', 19.344835155124766), ('B', 4.660252861131004)] # Nelder-Mead
        #pars = [('I0', 1.842206038819525), ('gamma', 6250.77381318831), ('omega', 310308.7890994233), ('A', 19.344835155124766), ('B', 4.660252861131004)] # Nelder-Mead
    elif interpol == "7":
        ## RLC crit
        def func(xval, Q0, gamma):
            R = 54.77
            return R*Q0*gamma**2*xval*np.exp(-gamma*xval)
        pars = [('I0', 1), ('gamma', 1)]
    elif interpol == "8":
        ## RLC sovr
        def func(xval, Q0, gamma, omega):
            R = 1
            return R*Q0*(gamma**2+omega**2)/(2*omega)*(np.exp(-(gamma-omega)*xval)-np.exp(-(gamma+omega)*xval))
        pars = [('I0', 1), ('gamma', 1), ('omega', 1)]
    else:
        print("Scegliere una interpolazione.")
        sys.exit(1)

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
            print(f"{keys}: {values}")
        sys.exit(1)

    # Costruisce il modello
    fmodel = Model(func)
    #fmodel.set_param_hint('I0', value=1, min=0., max=10.)
    #fmodel.set_param_hint('I0', value=1, min=0)
    for i, j in pars:
        fmodel.set_param_hint(i, value=j)
    params=fmodel.make_params()

    # Interpolazione
    if error=="0":
        result = fmodel.fit(ydata, params, xval=xdata, method=method)
    elif error=="1":
        result = fmodel.fit(ydata, params, xval=xdata, weights=1/yerr, method=method)
    else:
        print("L'argomento dell'incertezza dev'essere 0 oppure 1.")
        sys.exit(1)

    print(result.fit_report())

    # Stampa nuovamente i parametri con più cifre significative
    best_params=np.array(list(dict.values(result.best_values)))
    try:
        clipboard = Tk()
        clipboard.withdraw()
        clipboard.clipboard_clear()
        print("\nParametri:")
        for i, j, k in zip([item[0] for item in list(dict.items(params))], best_params, np.sqrt(np.diag(result.covar))):
            print(f"{i}: {j} +/- {k}")
            clipboard.clipboard_append(repr(j) + " ")
    except:
        print("\nNessuna incertezza.")

    clipboard.update()

    # Stampa altri dati riguardo l'interpolazione
    print(f"Degrees of freedom = {len(xdata)-len(params)}")
    print(f"p-value = {1-stats.chi2.cdf(result.chisqr,len(xdata)-len(params))}")

    print(f"Matrice di covarianza:\n{result.covar}")

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

    csvtotxt.csvtotxt(dati, file_name)

if __name__ == '__main__':
    if (len(sys.argv) != 6):
        print("Aspettati cinque argomenti: nome file, numero dati, numero interpolazione, incertezze y (0/1), metodo interpolazione (1,2,3,...).")
        sys.exit(1)

    # Nome file, numero dati, numero interpolazione, incertezze y (0/1), metodo interpolazione (1,2,3,...)
    fit_non_lineare(*sys.argv[1:])