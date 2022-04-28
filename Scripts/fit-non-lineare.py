#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import sys

import csvtotxt

from scipy import stats
from lmfit import Model
from tkinter import Tk

def fit_non_lineare(file_name, dati, interpol, error, method):
    # Legge x, y e yerr da file usando numero dati
    col_list = [f"x {dati}", f"y {dati}", f"yerr {dati}"]
    df = pd.read_csv(file_name, usecols=col_list, sep=",", decimal=".")

    # Separa il dataframe nei rispettivi array
    df = df.dropna()
    xdata = df[f'x {dati}']
    ydata = df[f'y {dati}']
    yerr = df[f'yerr {dati}']
    xdata = np.array(xdata, dtype=np.float128)
    ydata = np.array(ydata, dtype=np.float128)
    yerr = np.array(yerr, dtype=np.float128)

    #### Interpolazioni ####
    if interpol=="1":
        ## prima interpolazione ## 
        def func(xval, I0, fattore, A):
            return I0*(np.exp(fattore*xval)-1)+A
        pars = [('I0', 1), ('fattore', 38.6), ('A', 0)]
    elif interpol == "3":
        ## RC
        def func(xval, V_g, tau, A):
            periodo = 1./.8 # secondi
            return V_g*(1-2/(1+np.exp(-periodo/(2*tau)))*np.exp(-xval/tau))+A
        pars = [('V_g', 9.875198121607973), ('tau', 0.11148566083571734), ('A', 9.826661090020343)] # Nelder-Mead
    elif interpol == "4":
        ## RL
        def func(xval, V_g, tau, A):
            #periodo = 1/100
            periodo = 10 # milli secondi
            return V_g*(2*np.exp(-xval/tau)/(1+np.exp(-periodo/(2*tau))))+A
        pars = [('V_g', 0.9763507751615479), ('tau', 0.07711593596900626), ('A', 0.11063333194641872)] # Nelder-Mead
    elif interpol == "5":
        ## RL reale
        def func(xval, V_g, R_L, tau, A):
            periodo = 1/100 # secondi
            R = 992+50 # ohm
            return V_g*(2*np.exp(-xval/tau)/(1+np.exp(-periodo/(2*tau)))+R_L/(R_L+R)*(1-2*np.exp(-xval/tau)/(1+np.exp(-periodo/(2*tau)))))+A
        pars = [('V_g', 1.0323245746163687), ('R_L', 57.769650453274096), ('tau', 7.711618152518606e-05), ('A', 0.05382326137423557)] # Nelder-Mead
    elif interpol == "6":
        ## RLC sotto
        def func(xval, I0, gamma, omega, A, B):
            R = 1#R = 10+50+60 # ohm
            return R*I0*np.exp(-gamma*xval)*np.sin(omega*xval+B)+A
        pars = [('I0', 100), ('gamma', 1.2), ('omega', 5.22), ('A', -2), ('B', 0)] # Nelder-Mead
    elif interpol == "7":
        ## RLC crit
        def func(xval, gamma, A, B):
            R = 1 # ohm
            return R*A*xval*np.exp(-gamma*xval)+B
        pars = [('gamma', 5), ('A', 17), ('B', 0)]
    elif interpol == "8":
        ## RLC sovr
        def func(xval, A, B, C, D, E):
            R = 1 # ohm
            return R*A*np.exp(B*xval)-R*C*np.exp(D*xval)+E
        pars = [('A', 2), ('B', -1), ('C', 3), ('D', -1), ('E', 0)]
        # A = I0*B, B = -gamma+beta,    C = I0*C,   D=-gamma-beta
        # Le lettere dopo ogni "=" si riferiscono al return seguente
        #def func(xval, I0, gamma, beta, A, B, C):
            # R = 1
            #return R*I0*np.exp(-gamma*xval)*(B*np.exp(beta*xval)-C*np.exp(-beta*xval))+A
        #pars = [('I0', 6), ('gamma', 39), ('beta', 38), ('A', 0), ('B', 1), ('C', 1)]
    elif interpol == "9":
        #xdata *= 2*np.pi # se xdata è la frequenza dell'onda, bisogna propagare le incertezze
        #xdata = 2*np.pi/xdata # se xdata è il periodo dell'onda, bisogna propagare le incertezze
        ## Trasferimento da V_g a V_C. Modulo. Circuito RC
        def func(xval, C): # xval = omega
            V_g = 1 # volt
            R = 2002 # ohm
            return V_g/np.sqrt(1+(xval*R*C)*(xval*R*C))
        pars = [('C', 1e-6)]
    elif interpol == "10":
        ## Trasferimento da V_g a V_L. Modulo. Circuito RL
        def func(xval, L): # xval = omega
            V_g = 1 # volt
            R = 992 # ohm
            return V_g*xval*L/np.sqrt(R*R+(xval*L)*(xval*L))
        pars = [('L', 1)]
    elif interpol == "11":
        ## Trasferimento da V_g a V_L. Modulo. Circuito RL reale
        def func(xval, L, R_L): # xval = omega
            V_g = 1 # volt
            R = 992 # ohm
            return V_g*np.sqrt((R_L*R_L+(xval*L)*(xval*L))/((R+R_L)*(R+R_L)+(xval*L)*(xval*L)))
        pars = [('L', 1), ('R_L', 1)]
    elif interpol == "12":
        ## Trasferimento da V_g a V_R. Modulo. Circuito RLC
        def func(xval, C, L): # xval = omega
            V_g = 1 # volt
            R = 992 # ohm
            return V_g*R/np.sqrt(R*R+np.square(xval*L-1/(xval*C)))
        pars = [('C', 1), ('L', 1)]
    elif interpol == "13":
        ## Trasferimento da V_g a V_C. Modulo. Circuito RLC
        def func(xval, C, L): # xval = omega
            V_g = 1 # volt
            R = 992 # ohm
            #return V_g/np.sqrt(np.square(xval*C*R)+np.square((xval*xval*C*L-1)))
            return V_g*1/(xval*C)*1/np.sqrt(R*R+np.square(xval*L-1/(xval*C)))
        pars = [('C', 1), ('L', 1)]
    elif interpol == "14":
        ## Trasferimento da V_g a V_L. Modulo. Circuito RLC
        def func(xval, C, L): # xval = omega
            V_g = 1 # volt
            R = 992 # ohm
            return V_g*xval*L/np.sqrt(R*R+np.square(xval*L-1/(xval*C)))
        pars = [('C', 1), ('L', 1)]
    elif interpol == "15":
        ## Trasferimento da V_g a V_R. Come 12, R_L reale
        def func(xval, C, L, R_L): # xval = omega
            V_g = 1 # volt
            R = 992 # ohm
            return V_g*R/np.sqrt(np.square(R+R_L)+np.square(xval*L-1/(xval*C)))
        pars = [('C', 1), ('L', 1), ('R_L', 1)]
    elif interpol == "16":
        ## Trasferimento da V_g a V_C. Come 13, R_L reale
        def func(xval, C, L, R_L): # xval = omega
            V_g = 1 # volt
            R = 992 # ohm
            #return V_g/np.sqrt(np.square(xval*C*R)+np.square((xval*xval*C*L-1)))
            return V_g*1/(xval*C)*1/np.sqrt((np.square(R+R_L)+np.square(xval*L-1/(xval*C))))
        pars = [('C', 1), ('L', 1), ('R_L', 1)]
    elif interpol == "17":
        ## Trasferimento da V_g a V_L. Come 14, R_L reale
        def func(xval, C, L, R_L): # xval = omega
            V_g = 1 # volt
            R = 992 # ohm
            return V_g*np.sqrt((R_L*R_L+np.square(xval*L))/(np.square(R+R_L)+np.square(xval*L-1/(xval*C))))
        pars = [('C', 1), ('L', 1), ('R_L', 1)]
    #elif interpol == "18":
    #    def func(xval):
    #        return xval
    #    pars = []
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
    plt.figure(figsize=(30,15))
    plt.axhline(y=0)
    plt.plot(np.array(xdata), result.residual, 'bo')
    plt.show()
    plt.close()

    ## Normality tests for residuals
    sm.qqplot(result.residual, fit=True, line='s')
    plt.show()
    print(stats.shapiro(result.residual))

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
        plt.savefig('fit_python.png', bbox_inches='tight')
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