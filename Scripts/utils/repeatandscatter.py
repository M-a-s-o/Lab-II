import numpy as np
import propagate as pr

class RepeatData():
    r"""
        file : string
            String with file path specified.

        dati : list of int, tuple of int
            Start and end of scatter points.

        variables : int
            Number of variables in each set of repeated data x y z etc.
    """
    def __init__(self, file, dati, variables) -> None:
        self.variables = variables
        self.data = pr.LoadData(file, dati, variables)
    
    def get_mean(self, var):
        return np.mean(self.data[var])

    def get_means(self):
        return np.array([self.get_mean(i) for i in range(self.variables)])
    
    def get_stdev(self, var):
        return np.std(self.data[var], ddof = 1)/np.sqrt(len(self.data[var]))
    
    def get_stdevs(self):
        return np.array([self.get_stdev(i) for i in range(self.variables)])

    def get_covar_mat(self):
        return np.cov(self.data)/len(self.data[0])
    
    def get_covar_mat_part(self, *var):
        mat = self.get_covar_mat()
        if isinstance(var[0], list):
            var = var[0]
        mat = np.take(mat, list(var), axis=0)
        mat = np.take(mat, list(var), axis=1)
        return mat
    
    def get_covar(self, *var):
        return self.get_covar_mat()[var[0], var[1]]

class ScatterPointsRepeat():
    r"""
        file : string
            String with file path specified.

        dati : list of int, tuple of int
            Start and end of scatter points.

        variables : int
            Number of variables in each set of repeated data x y z etc.
    """
    def __init__(self, file, dati, variables) -> None:
        dati = (dati[0], dati[1]+1)
        self.columns = np.array([RepeatData(file, i, variables) for i in range(*tuple(dati))])
        self.means = np.array([col.get_means() for col in self.columns])
        self.stdevs = np.array([col.get_stdevs() for col in self.columns])
        self.covars = np.array([col.get_covar_mat() for col in self.columns])