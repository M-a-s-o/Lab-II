import pandas as pd
import sys

if (len(sys.argv) != 2):
    print("Specificare colonna dati.")
    sys.exit(1)

dati = sys.argv[1]

# read x, y and y error from file using run number
col_list = [f"x {dati}", f"xerr {dati}", f"y {dati}", f"yerr {dati}"]
df = pd.read_csv("Dati.csv", usecols=col_list, sep=",", decimal=".")

# manipulate data into numpy array
df = df.dropna()
#df = df.sort_values(by=df.columns[0]).reset_index()
df = df.sort_values([df.columns[0], df.columns[1]]).reset_index()
xdata = df[f"x {dati}"]
ydata = df[f"y {dati}"]
xerr = df[f"xerr {dati}"]
yerr = df[f"yerr {dati}"]

file = open("data.txt", "w")
for i in range(len(xdata)):
    file.write(repr(xdata[i]) + " " + repr(ydata[i]) + " " + repr(xerr[i]) + " " + repr(yerr[i]) + "\n")
file.close