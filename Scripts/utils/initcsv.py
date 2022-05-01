from tkinter import Tk
import sys

def initcsv(start_datum, end_datum, num_vars, name=None):
    vars_string = ""
    if name is not None:
        file = open(name, "w")
        vars_string += "\t"
    for var in range(1, num_vars+1):
        #vars_string += ","
        vars_string += "x{1} {0}\t".format("{0}", repr(var))
    cols_names_str = ""
    for datum in range(start_datum, end_datum+1):
        cols_names_str += vars_string.format(repr(datum)) + "\t"
    if name is not None:
        file.write(cols_names_str)
        file.close
    cp = Tk()
    cp.withdraw()
    cp.clipboard_clear()
    cp.clipboard_append(cols_names_str)
    cp.update()

if __name__ == '__main__':
    if len(sys.argv) > 5 or len(sys.argv) < 2:
        print("Aspettati tre/quattro argomenti: inizio dati, fine dati, numero di variabili, nome file (opzionale).")
        sys.exit(1)
    elif len(sys.argv) == 4:
        arr = [int(x) for x in sys.argv[1:4]]
        initcsv(*arr)
    else:
        arr = [int(x) for x in sys.argv[1:4]]
        initcsv(*arr, sys.argv[4])