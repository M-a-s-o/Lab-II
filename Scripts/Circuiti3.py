import numpy as np
import sys

from utils import propagate as pr
from utils.repeatandscatter import ScatterPointsRepeat

class RepeatDataCircuiti3(ScatterPointsRepeat):
    r"""
        Format
        ------
        Per ogni frequenza, i dati ripetuti sono raccolti in cinque colonne in un file chiamato \"dstd.csv\": V_A, V_B, V_(A-B), Delta phi\', Delta phi\'\'. Ogni colonna ha titolo x1 #, x2 #, x3 #, x4 #, x5 #; dove # è un numero ordinale che identifica la frequenza a cui si sono raccolti i dati.

        Output
        ------
        La classe calcola direttamente quanto desiderato e fornisce l'output in appositi file (temporanei) come due colonne: valore medio e sua deviazione standard.
    """
    def __init__(self, file, dati, variables) -> None:
        super().__init__(file, dati, variables)
        self.dati = dati

        thrshld = 1.58e-16
        std_volt = self.stdevs.T[0:3]
        # check for mean > 1
        std_volt = np.where((self.means.T[0:3] >= 1) & (std_volt <= thrshld), 0.01/np.sqrt(3), std_volt)
        # check for zero std dev
        std_volt = np.where(std_volt <= thrshld, 0.001/np.sqrt(3), std_volt)

        std_ang = self.stdevs.T[3:5]
        # check for mean > 10
        std_ang = np.where((np.abs(self.means.T[3:5]) < 10) & (std_ang <= thrshld), 0.01/np.sqrt(3), std_ang)
        # check for zero std dev
        std_ang = np.where(std_ang <= thrshld, 0.1/np.sqrt(3), std_ang)
        self.stdevs = np.concatenate((std_volt.T, std_ang.T), axis=1)

    def propagate(self, dati, *var):
        x = pr.GenerateXs(2)
        func = x[1]/x[0]
        means = [self.means[dati-self.dati[0], i-1] for i in var]
        covar = self.columns[dati-self.dati[0]].get_covar_mat_part([i-1 for i in var])
        variance = np.square(np.array([self.stdevs[dati-self.dati[0], i-1] for i in var]))
        np.fill_diagonal(covar, variance)
        return pr.Propagate(means, covar, 2, func, x)

    def get_VAB_on_VA(self):
        return np.array([self.propagate(i, 1, 3) for i in range(*tuple([self.dati[0], self.dati[1]+1]))])
    
    def get_VB_on_VA(self):
        return np.array([self.propagate(i, 1, 2) for i in range(*tuple([self.dati[0], self.dati[1]+1]))])
    
    def get_phi_VAB_wrt_VA(self): # Delta phi prime
        return np.radians(np.array([self.means.T[3], self.stdevs.T[3]]).T)

    def get_phi_VB_wrt_VA(self): # Delta phi double prime
        return np.radians(np.array([self.means.T[4], self.stdevs.T[4]]).T)

def main(start, end = None):
    if end is None:
        if start in range(1, 6):
            RawDataColumns = {
                1: (1, 12),
                2: (13, 25),
                3: (26, 34),
                4: (35, 44),
                5: (45, 54),
            }
            start, end = RawDataColumns[start]
        else:
            print("Valori preimpostati da 1 a 5.")
            sys.exit(1)
    
    dati = RepeatDataCircuiti3("../Circuiti-3/dstd.csv", [start, end], 5)
    np.savetxt("../Circuiti-3/VAB-VA.txt", dati.get_VAB_on_VA(), delimiter="\t")
    np.savetxt("../Circuiti-3/VB-VA.txt", dati.get_VB_on_VA(), delimiter="\t")
    np.savetxt("../Circuiti-3/phi-VAB-VA.txt", dati.get_phi_VAB_wrt_VA(), delimiter="\t")
    np.savetxt("../Circuiti-3/phi-VB-VA.txt", dati.get_phi_VB_wrt_VA(), delimiter="\t")

if __name__ == '__main__':
    #if (len(sys.argv) < 3):
    #    print("Due argomenti: inizio e fine dei set di dati di cui calcolare VAB/VA e VB/VA.")
    #    sys.exit(1)
    if len(sys.argv) not in range(2,4):
        print("Due argomenti: inizio e fine dei set di dati di cui calcolare VAB/VA e VB/VA.")
        print("Oppure un argomento con inizio e fine già preimpostate.")
        sys.exit(1)
    main(*[int(x) for x in sys.argv[1:]])