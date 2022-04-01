import numpy as np
import matplotlib.pyplot as plt
import sympy as smp
from scipy import stats
from lmfit import Model
import pandas as pd
import sys

# metodi
metodi = {
    1: "leastsq",
    2: "least_squares",
    3: "nelder",
    4: "emcee",
}

def LoadDataframe(file, dati):
    if file is None:
        print("Defaulting to \"./Dati.csv\" file.")
        file = "./Dati.csv"
    # read x, y and y error from file using run number
    col_list = [f"x {dati}", f"y {dati}", f"yerr {dati}"]
    df = pd.read_csv(file, usecols=col_list, sep=",", decimal=".")
    
    # manipulate data into numpy array
    df = df.dropna()
    xdata = df[f'x {dati}']
    ydata = df[f'y {dati}']
    yerr = df[f'yerr {dati}']
    return xdata, ydata, yerr

def DefineFunction(interpol, func):
    #### fit ####
    if interpol == 0:
        return func
    elif interpol == 1:
        def func(x, B):
            return B*x
    elif interpol == 2:
        def func(x, A, B):
            return A + B*x
    elif interpol == 3:
        def func(x, I0, fattore):
            return I0*(np.exp(fattore*x)-1)
    else:
        print("Funzione non disponibile.")
        return
    return func

def SetParameters(fmodel, interpol, par):
    if interpol == 0:
        for i,j in par:
            fmodel.set_param_hint(i, value = j)
    elif interpol == 1:
        fmodel.set_param_hint('B', value=1)
    elif interpol == 2:
        fmodel.set_param_hint('A', value=0)
        fmodel.set_param_hint('B', value=1)
    elif interpol == 3:
        fmodel.set_param_hint('I0', value = 1)
        fmodel.set_param_hint('fattore', value = 38.6)
    else:
        print("Interpolazione non disponibile.")
        return

    return fmodel.make_params()

def Plot(xdata, ydata, yerr, func, best_params, show, name):
    #### plot ####
    xdum = np.linspace(np.amin(xdata),np.amax(xdata), 1000)
    plt.figure(figsize=(30,15))
    #plt.plot(xdata, ydata, 'b*', label='data')
    plt.errorbar(xdata, ydata, yerr=yerr, fmt='bo', label='data')
    #plt.plot(xdata, result.best_fit, 'r-', label='best fit')
    plt.plot(xdum, func(xdum, *best_params), 'g--', label='best fit')
    plt.xlabel('xdata')
    plt.ylabel('ydata')
    plt.legend()

    try:
        plt.savefig(name, bbox_inches='tight')
        if show:
            plt.show()
    except:
        pass
    plt.close()

# Aspettati quattro argomenti: run number, numero interpolazione (1, 2, 3, ...), incertezza y (0/1), metodo interpolazione (1,2,3,...)
# run number: numero corrispondente al set di dati
# numero interpolazione: 0 è l'interpolazione a piacere, 1 è proporzionalità diretta, 2 è lineare, 3 legge di Shockley
# incertezza sulle y: sì o no
# metodo dell'interpolazione tramite lmfit: 1 è Levenberg-Marquardt, 2 è least_squares, 3 è Nelder-Mead
def Fitter(dati, interpol, error, file = "Dati.csv", method = 1, func=None, par=None, show=False, name="fit_func.png"):
    if method in metodi:
        method = metodi[method]
    else:
        print("Scegliere un metodo:")
        for keys, values in metodi.items():
            print(f"{keys}: {values}")
        return
    
    xdata, ydata, yerr = LoadDataframe(file, dati)
    func = DefineFunction(interpol, func)

    # set model's parameters
    fmodel = Model(func)
    par = [('I0', 1), ('fattore', 38.6)]
    params = SetParameters(fmodel, interpol, par)

    # do fit, with either leastsq, least_squares or nelder
    if error == 0:
        result = fmodel.fit(ydata, params, x=xdata, method=method)
    elif error == 1:
        result = fmodel.fit(ydata, params, x=xdata, weights=1/yerr, method=method)
    else:
        print("L'argomento dev'essere 0 oppure 1.")
        return

    print(result.fit_report())

    best_params = np.array(list(dict.values(result.best_values)))
    # dict.items() for array of tuples

    try:
        print("\nParametri:")
        for i, j, k in zip([item[0] for item in list(dict.items(params))], best_params, np.sqrt(np.diag(result.covar))):
            print(f"{i}: {j} +/- {k}")
    except:
        print("\nNo uncertainties.")

    print(f"Matrice di covarianza:\n{result.covar}")
    print(f"Degrees of freedom = {len(xdata)-len(params)}")
    print(f"p-value = {1-stats.chi2.cdf(result.chisqr,len(xdata)-len(params))}")

    Plot(xdata, ydata, yerr, func, best_params, show, name)

def LoadData(file, dati, variables):
    if file is None:
        print("Defaulting to \"./dstd.csv\" file.")
        file = "./dstd.csv"
    # read measurements from file using run number
    col_list = ["x{} {}".format(i+1, dati) for i in range(variables)]
    df = pd.read_csv(file, usecols=col_list, sep=";", decimal=".")

    # manipulate data into numpy array
    df = df.dropna()
    return np.array([df['x{} {}'.format(i+1, dati)] for i in range(variables)])

def GenerateXs(variables):
    return [smp.symbols('x%d' % i) for i in range(variables)]

def PropagateFile(file, dati, variables = 2, func=None, x=None):
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

def Propagate(mean, cov, variables = 2, func=None, x=None):
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
    print(f"mean measures: {mean}\nstd dev: {np.sqrt(np.diag(cov))}\nderivatives: {coeff}\nmean function: {mean_f}\nstd dev: {dstd}")

def CovarX (var, cov):
    if var == cov.shape[0]:
        return np.diag(cov)
    X = np.zeros((var,var))
    X[np.triu_indices(X.shape[0], k=0)] = cov
    X=X+X.T-np.diag(np.diag(X))
    return X

if __name__ == "__main__":
    fit = int(sys.argv[1])
    if fit == 1:
        Fitter(int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[2], int(sys.argv[6]))
    elif fit == 0:
        pass
    else:
        pass