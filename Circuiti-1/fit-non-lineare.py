#!/usr/bin/env python3

# 2022 - 03 - 27
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from lmfit import Model
import pandas as pd
import sys
from tkinter import Tk
import statsmodels.api as sm

if (len(sys.argv) != 5):
    print("Aspettati quattro argomenti: numero dati, numero interpolazione, incertezze y (0/1), metodo interpolazione (1,2,3,...).")
    sys.exit(1)

# Numero dati, numero interpolazione, incertezze y (0/1), metodo interpolazione (1,2,3,...)
dati = sys.argv[1]
interpol = sys.argv[2]
error = sys.argv[3]
method = sys.argv[4]

# Legge x, y e yerr da file usando numero dati
col_list = ["x {}".format(dati), "y {}".format(dati), "yerr {}".format(dati)] 
df = pd.read_csv("Dati.csv", usecols=col_list, sep=",", decimal=".")

# Separa il dataframe nei rispettivi array
df = df.dropna()
xdata = df["x {}".format(dati)]
ydata = df["y {}".format(dati)]
yerr = df["yerr {}".format(dati)]

#### Interpolazioni ####
if interpol=="1":
    ## prima interpolazione ## 
    def func(tensione, I0, fattore, A):
        return I0*(np.exp(fattore*tensione)-1)+A
else:
    print("Scegliere una interpolazione.")
    sys.exit(1)

# metodi
metodi = {
    "1": "leastsq",
    "2": "least_squares",
    "3": "nelder",
	"4": "emcee",
}

if method in metodi:
	method = metodi[method]
else:
	print("Scegliere un metodo:")
	for keys, values in metodi.items():
		print("{}: {}".format(keys,values))
	sys.exit(1)

# Costruisce il modello
fmodel = Model(func)
fmodel.set_param_hint('I0', value=1, min=0., max=10.)
params=fmodel.make_params(fattore=38.6, A=0)

# Interpolazione
if error=="0":
    result = fmodel.fit(ydata, params, tensione=xdata, method=method)
elif error=="1":
    result = fmodel.fit(ydata, params, tensione=xdata, weights=1/yerr, method=method)
else:
    print("L'argomento dell'incertezza dev'essere 0 oppure 1.")
    sys.exit(1)

print(result.fit_report())

# Stampa nuovamente i parametri con più cifre significative
best_params=np.array(list(dict.values(result.best_values)))
try:
    clipboard = Tk()
    clipboard.clipboard_clear()
    print("\nParametri:")
    for i, j, k in zip([item[0] for item in list(dict.items(params))], best_params, np.sqrt(np.diag(result.covar))):
        print("{}: {} +/- {}".format(i, j, k))
        clipboard.clipboard_append(repr(j) + " ")
except:
    print("\nNessuna incertezza.")

clipboard.update()

# Stampa altri dati riguardo l'interpolazione
print("\nGradi di libertà = {}".format(len(xdata)-len(params)))
print("p-value chi2 = {}".format(1-stats.chi2.cdf(result.chisqr,len(xdata)-len(params))))

print("Matrice di covarianza:\n{}".format(result.covar))

# Plot absolute residuals
print(result.residual)
plt.figure(figsize=(30,15))
result.plot_residuals()
plt.show()
plt.close()

# Plot normalised residuals
plt.figure(figsize=(30,15))
plt.axhline(y=0)
plt.plot(np.array(xdata), result.residual, 'bo')
plt.show()
plt.close()

# Normality tests for residuals
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
    plt.savefig('fit_func.png', bbox_inches='tight')
    plt.show()
except:
    pass
plt.close()
