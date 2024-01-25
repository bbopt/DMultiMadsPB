import copy
from joblib import Parallel, delayed
import matlab.engine
from multiprocessing import Pool
import math
import numpy as np
from tqdm import tqdm

# Require pymoo 0.4.2.2
from pymoo.model.problem import Problem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.optimize import minimize

def run_nsgaii_solver_pb(name):
    # Initialize engine
    eng = matlab.engine.start_matlab()
    eng.addpath("../../problems")

    # initialize dictionnary
    dict_problems = dict()

    # BK1
    dict_problems["BK1"] = (eng.BK1, 2, 2, np.array([-5, -5]), np.array([10, 10]))

    # CL1
    F_p = 10
    sigma_p = 10
    lb = np.array(
            [
                F_p / sigma_p,
                math.sqrt(2) * F_p / sigma_p,
                math.sqrt(2) * F_p / sigma_p,
                F_p / sigma_p,
                ]
            )
    dict_problems["CL1"] = (eng.CL1, 4, 2, np.copy(lb), (3 * F_p / sigma_p) * np.ones(4))

    # Deb41 function
    lb = np.array([0.1, 0.0])
    dict_problems["Deb41"] = (eng.Deb41, 2, 2, lb, np.ones(2))

    # Deb512a function
    dict_problems["Deb512a"] = (eng.Deb512a, 2, 2, np.zeros(2), np.ones(2))

    # Deb512b function
    dict_problems["Deb512b"] = (eng.Deb512b, 2, 2, np.zeros(2), np.ones(2))

    # Deb512c function
    dict_problems["Deb512c"] = (eng.Deb512c, 2, 2, np.zeros(2), np.ones(2))

    # Deb513 function:
    dict_problems["Deb513"] = (eng.Deb513, 2, 2, np.zeros(2), np.ones(2))

    # Deb521a function
    dict_problems["Deb521a"] = (eng.Deb521a, 2, 2, np.zeros(2), np.ones(2))

    # Deb521b function
    dict_problems["Deb521b"] = (eng.Deb521b, 2, 2, np.zeros(2), np.ones(2))

    # Deb53 function
    dict_problems["Deb53"] = (eng.Deb53, 2, 2, np.zeros(2), np.ones(2))

    # DG01 function
    dict_problems["DG01"] = (eng.DG01, 1, 2, -10 * np.ones(1), 13 * np.ones(1))

    # DPAM1 function
    dict_problems["DPAM1"] = (eng.DPAM1, 10, 2, -0.3 * np.ones(10), 0.3 * np.ones(10))

    # DTLZ1 function
    dict_problems["DTLZ1"] = (eng.DTLZ1, 7, 3, np.zeros(7), np.ones(7))

    # DTLZ1n2 function
    dict_problems["DTLZ1n2"] = (eng.DTLZ1n2, 2, 2, np.zeros(2), np.ones(2))

    # DTLZ2 function
    dict_problems["DTLZ2"] = (eng.DTLZ2, 12, 3, np.zeros(12), np.ones(12))

    # DTLZn2 function
    dict_problems["DTLZ2n2"] = (eng.DTLZ2n2, 2, 2, np.zeros(2), np.ones(2))

    # DTLZ3 function
    dict_problems["DTLZ3"] = (eng.DTLZ3, 12, 3, np.zeros(12), np.ones(12))

    # DTLZ3n2 function
    dict_problems["DTLZ3n2"] = (eng.DTLZ3n2, 2, 2, np.zeros(2), np.ones(2))

    # DTLZ4 function
    dict_problems["DTLZ4"] = (eng.DTLZ4, 12, 3, np.zeros(12), np.ones(12))

    # DTLZ4n2 function
    dict_problems["DTLZ4n2"] = (eng.DTLZ4n2, 2, 2, np.zeros(2), np.ones(2))

    # DTLZ5 function
    dict_problems["DTLZ5"] = (eng.DTLZ5, 12, 3, np.zeros(12), np.ones(12))

    # DTLZ5n2 function
    dict_problems["DTLZ5n2"] = (eng.DTLZ5n2, 2, 2, np.zeros(2), np.ones(2))

    # DTLZ6 function
    dict_problems["DTLZ6"] = (eng.DTLZ6, 22, 3, np.zeros(22), np.ones(22))

    # DTLZ6n2 function
    dict_problems["DTLZ6n2"] = (eng.DTLZ6n2, 2, 2, np.zeros(2), np.ones(2))

    # ex005 function
    dict_problems["ex005"] = (eng.ex005, 2, 2, np.array([-1, 1]), np.array([2, 2]))

    # Far1 function
    dict_problems["Far1"] = (eng.Far1, 2, 2, -1 * np.ones(2), np.ones(2))

    # FES1 function
    dict_problems["FES1"] = (eng.FES1, 10, 2, np.zeros(10), np.ones(10))

    # FES2 function
    dict_problems["FES2"] = (eng.FES2, 10, 3, np.zeros(10), np.ones(10))

    # FES3 function
    dict_problems["FES3"] = (eng.FES3, 10, 4, np.zeros(10), np.ones(10))

    # Fonseca function
    dict_problems["Fonseca"] = (eng.Fonseca, 2, 2, -4 * np.ones(2), 4 * np.ones(2))

    # I1 function
    dict_problems["I1"] = (eng.I1, 8, 3, np.zeros(8), np.ones(8))

    # I2 function
    dict_problems["I2"] = (eng.I2, 8, 3, np.zeros(8), np.ones(8))

    # I3 function
    dict_problems["I3"] = (eng.I3, 8, 3, np.zeros(8), np.ones(8))

    # I4 function
    dict_problems["I4"] = (eng.I4, 8, 3, np.zeros(8), np.ones(8))

    # I5 function
    dict_problems["I5"] = (eng.I5, 8, 3, np.zeros(8), np.ones(8))

    # IKK1 function
    dict_problems["IKK1"] = (eng.IKK1, 2, 3, -50 * np.ones(2), 50 * np.ones(2))

    # IM1 function
    dict_problems["IM1"] = (eng.IM1, 2, 2, np.ones(2), np.array([4, 2]))

    # Jin1 function
    dict_problems["Jin1"] = (eng.Jin1, 2, 2, np.zeros(2), np.ones(2))

    # Jin2 function
    dict_problems["Jin2"] = (eng.Jin2, 2, 2, np.zeros(2), np.ones(2))

    # Jin3 function
    dict_problems["Jin3"] = (eng.Jin3, 2, 2, np.zeros(2), np.ones(2))

    # Jin4 function
    dict_problems["Jin4"] = (eng.Jin4, 2, 2, np.zeros(2), np.ones(2))

    # Kursawe function
    dict_problems["Kursawe"] = (eng.Kursawe, 3, 2, -5 * np.ones(3), 5 * np.ones(3))

    # L1ZDT4 function
    dict_problems["L1ZDT4"] = (
            eng.L1ZDT4,
            10,
            1,
            np.concatenate([np.array([0]), -5 * np.ones(9)]),
            np.concatenate([np.array([1]), 5 * np.ones(9)]),
            )

    # L2ZDT1 function
    dict_problems["L2ZDT1"] = (eng.L2ZDT1, 30, 2, np.zeros(30), np.ones(30))

    # L2ZDT2 function
    dict_problems["L2ZDT2"] = (eng.L2ZDT2, 30, 2, np.zeros(30), np.ones(30))

    # L2ZDT3 function
    dict_problems["L2ZDT3"] = (eng.L2ZDT3, 30, 2, np.zeros(30), np.ones(30))

    # L2ZDT4 function
    dict_problems["L2ZDT4"] = (eng.L2ZDT4, 30, 2, np.zeros(30), np.ones(30))

    # L2ZDT6 function
    dict_problems["L2ZDT6"] = (eng.L2ZDT6, 10, 2, np.zeros(10), np.ones(10))

    # L3ZDT1 function
    dict_problems["L3ZDT1"] = (eng.L3ZDT1, 30, 2, np.zeros(30), np.ones(30))

    # L3ZDT2 function
    dict_problems["L3ZDT2"] = (eng.L3ZDT2, 30, 2, np.zeros(30), np.ones(30))

    # L3ZDT3 function
    dict_problems["L3ZDT3"] = (eng.L3ZDT3, 30, 2, np.zeros(30), np.ones(30))

    # L3ZDT4 function
    dict_problems["L3ZDT4"] = (eng.L3ZDT4, 30, 2, np.zeros(30), np.ones(30))

    # L3ZDT6 function
    dict_problems["L3ZDT6"] = (eng.L3ZDT6, 10, 2, np.zeros(10), np.ones(10))

    # LE1 function
    dict_problems["LE1"] = (eng.LE1, 2, 2, np.zeros(2), np.ones(2))

    # lovison1 function
    dict_problems["lovison1"] = (eng.lovison1, 2, 2, np.zeros(2), 3 * np.ones(2))

    # lovison2 function
    dict_problems["lovison2"] = (eng.lovison2, 2, 2, -0.5 * np.ones(2), np.array([0, 0.5]))

    # lovison3 function
    dict_problems["lovison3"] = (eng.lovison3, 2, 2, np.array([0, -4]), np.array([6, 4]))

    # lovison4 function
    dict_problems["lovison4"] = (eng.lovison4, 2, 2, np.array([0, -1]), np.array([6, 1]))

    # lovison5 function
    dict_problems["lovison5"] = (eng.lovison5, 3, 3, -1 * np.ones(3), 4 * np.ones(3))

    # lovison6 function
    dict_problems["lovison6"] = (eng.lovison6, 3, 3, -1 * np.ones(3), 4 * np.ones(3))

    # LRS1 function
    dict_problems["LRS1"] = (eng.LRS1, 2, 2, -50 * np.ones(2), 50 * np.ones(2))

    # MHHM1 function
    dict_problems["MHHM1"] = (eng.MHHM1, 1, 3, np.zeros(1), np.ones(1))

    # MHHM2 function
    dict_problems["MHHM2"] = (eng.MHHM2, 2, 3, np.zeros(2), np.ones(2))

    # MLF1 function
    dict_problems["MLF1"] = (eng.MLF1, 1, 2, np.zeros(1), 20 * np.ones(1))

    # MLF2 function
    dict_problems["MLF2"] = (eng.MLF2, 2, 2, -2 * np.ones(2), 2 * np.ones(2))

    # MOP1 function
    dict_problems["MOP1"] = (eng.MOP1, 1, 2, -10 ** (-5) * np.ones(1), 10 ** (5) * np.ones(1))

    # MOP2 function
    dict_problems["MOP2"] = (eng.MOP2, 4, 2, -4 * np.ones(4), 4 * np.ones(4))

    # MOP3 function
    dict_problems["MOP3"] = (eng.MOP3, 2, 2, -math.pi * np.ones(2), math.pi * np.ones(2))

    # MOP4 function
    dict_problems["MOP4"] = (eng.MOP4, 3, 2, -5 * np.ones(3), 5 * np.ones(3))

    # MOP5 function
    dict_problems["MOP5"] = (eng.MOP5, 2, 3, -30 * np.ones(2), 30 * np.ones(2))

    # MOP6 function
    dict_problems["MOP6"] = (eng.MOP6, 2, 2, np.zeros(2), np.ones(2))

    # MOP7 function
    dict_problems["MOP7"] = (eng.MOP7, 2, 3, -400 * np.ones(2), 400 * np.ones(2))

    # OKA1 function
    lb = np.array([6 * math.sin(math.pi / 12), -2 * math.pi * math.sin(math.pi / 12)])
    ub = np.array(
            [
                6 * math.sin(math.pi / 12) + 2 * math.pi * math.cos(math.pi / 12),
                6 * math.cos(math.pi / 12),
                ]
            )
    dict_problems["OKA1"] = (eng.OKA1, 2, 2, np.copy(lb), np.copy(ub))

    # OKA2 function
    dict_problems["OKA2"] = (
            eng.OKA2,
            3,
            2,
            np.array([-math.pi, -5, -5]),
            np.array([math.pi, 5, 5]),
            )

    # QV1 function
    dict_problems["QV1"] = (eng.QV1, 10, 2, -5.12 * np.ones(10), 5.12 * np.ones(10))

    # Sch1 function
    dict_problems["Sch1"] = (eng.Sch1, 1, 2, np.zeros(1), 5 * np.ones(1))

    # SK1 function
    dict_problems["SK1"] = (eng.SK1, 1, 2, np.array([-10]), np.array([10]))

    # SK2 function
    dict_problems["SK2"] = (eng.SK2, 4, 2, -10 * np.ones(4), 10 * np.ones(4))

    # SP1 function
    dict_problems["SP1"] = (eng.SP1, 2, 2, -1 * np.ones(2), 5 * np.ones(2))

    # SSFYY1 function
    dict_problems["SSFYY1"] = (eng.SSFYY1, 2, 2, -100 * np.ones(2), 100 * np.ones(2))

    # SSFYY2 function
    dict_problems["SSFYY2"] = (eng.SSFYY2, 1, 2, -100 * np.ones(1), 100 * np.ones(1))

    # TKLY1 function
    dict_problems["TKLY1"] = (eng.TKLY1, 4, 2, np.array([0.1, 0, 0, 0]), np.ones(4))

    # VFM1 function
    dict_problems["VFM1"] = (eng.VFM1, 2, 3, -2 * np.ones(2), 2 * np.ones(2))

    # VU1 function
    dict_problems["VU1"] = (eng.VU1, 2, 2, -3 * np.ones(2), 3 * np.ones(2))

    # VU2 function
    dict_problems["VU2"] = (eng.VU2, 2, 2, -3 * np.ones(2), 3 * np.ones(2))

    # WFG1 function
    dict_problems["WFG1"] = (eng.WFG1, 8, 3, np.zeros(8), 2 * np.linspace(1, 8, 8))

    # WFG2 function
    dict_problems["WFG2"] = (eng.WFG2, 8, 3, np.zeros(8), 2 * np.linspace(1, 8, 8))

    # WFG3 function
    dict_problems["WFG3"] = (eng.WFG3, 8, 3, np.zeros(8), 2 * np.linspace(1, 8, 8))

    # WFG4 function
    dict_problems["WFG4"] = (eng.WFG4, 8, 3, np.zeros(8), 2 * np.linspace(1, 8, 8))

    # WFG5 function
    dict_problems["WFG5"] = (eng.WFG5, 8, 3, np.zeros(8), 2 * np.linspace(1, 8, 8))

    # WFG6 function
    dict_problems["WFG6"] = (eng.WFG6, 8, 3, np.zeros(8), 2 * np.linspace(1, 8, 8))

    # WFG7 function
    dict_problems["WFG7"] = (eng.WFG7, 8, 3, np.zeros(8), 2 * np.linspace(1, 8, 8))

    # WFG8 function
    dict_problems["WFG8"] = (eng.WFG8, 8, 3, np.zeros(8), 2 * np.linspace(1, 8, 8))

    # WFG9 function
    dict_problems["WFG9"] = (eng.WFG9, 8, 3, np.zeros(8), 2 * np.linspace(1, 8, 8))

    # ZDT1 function
    dict_problems["ZDT1"] = (eng.ZDT1, 30, 2, np.zeros(30), np.ones(30))

    # ZDT2 function
    dict_problems["ZDT2"] = (eng.ZDT2, 30, 2, np.zeros(30), np.ones(30))

    # ZDT3 function
    dict_problems["ZDT3"] = (eng.ZDT3, 30, 2, np.zeros(30), np.ones(30))

    # ZDT4 function
    dict_problems["ZDT4"] = (
            eng.ZDT4,
            10,
            2,
            np.concatenate([np.array([0]), -5 * np.ones(9)]),
            np.concatenate([np.array([1]), 5 * np.ones(9)]),
            )

    # ZDT6 function
    dict_problems["ZDT6"] = (eng.ZDT6, 10, 2, np.zeros(10), np.ones(10))

    # ZLT1 function
    dict_problems["ZLT1"] = (eng.ZLT1, 10, 3, -1000 * np.ones(10), 1000 * np.ones(10))


    class MXProblem(Problem):
        def __init__(self, id_problem, type_constraints):
            super().__init__()

            # Initialize history data structures
            self.id_problem = id_problem
            self.x_cache = []
            self.f_cache = []
            self.c_cache = []

            # Add problems properties
            probproperties = dict_problems[id_problem]
            self.n_var = probproperties[1]
            self.n_obj = probproperties[2]
            if type_constraints in [1,2,5]:
                self.n_constr = self.n_var - 2
            elif type_constraints in [3,4]:
                self.n_constr = self.n_var - 1
            else:
                self.n_constr = 1

            self.xl = probproperties[3]
            self.xu = probproperties[4]

            # Type of constraints
            self.type_constraints = type_constraints


        # x represents a population array
        def _evaluate(self, x, out, *args, **kwargs):
            # Convert list of x to list of matlab element
            l_x_matlab = [
                    matlab.double(copy.copy(x[i, :].tolist())) for i in range(x.shape[0])
                    ]

            self.x_cache += x.tolist()

            # Get result of the matlab function for the list of x
            l_f_matlab = []
            for elt in l_x_matlab:
                l_f_matlab += [dict_problems[self.id_problem][0](eng.transpose(elt), nargout=1)]

            for elt in l_f_matlab:
                self.f_cache += [[elt[i][0] for i in range(elt.size[0])]]

            # Convert result to objective function f
            l_f = []
            for elt in l_f_matlab:
                l_f += [[elt[i][0] for i in range(elt.size[0])]]

            l_x_matlab = []
            l_f_matlab = []

            out["F"] = np.asarray(l_f)

            # Constraints
            # Put into an 2d array: each column represents a constraint <=
            n = x.shape[1]
            if self.type_constraints == 1:
                out["G"] = (3 - 2 * x[:, 1:n-1]) * x[:, 1:n-1] - x[:, 0:n-2] - 2 * x[:,2:n] + 1
            elif self.type_constraints == 2:
                out["G"] = (3 - 2 * x[:, 1:n-1]) * x[:, 1:n-1] - x[:, 0:n-2] - 2 * x[:,2:n] + 2.5
            elif self.type_constraints == 3:
                out["G"] = x[:, 0:n-1]**2 + x[:, 1:n]**2 + x[:, 0:n-1]*x[:,1:n] - 2 * x[:,0:n-1] - 2 * x[:,1:n] + 1
            elif self.type_constraints == 4:
                out["G"] = x[:, 0:n-1]**2 + x[:, 1:n]**2 + x[:, 0:n-1]*x[:,1:n] - 1
            elif self.type_constraints == 5:
                out["G"] = (3 - 0.5 * x[:, 1:n-1]) * x[:, 1:n-1] - x[:, 0:n-2] - 2 * x[:,2:n] + 1
            else: # self.type_constraints == 6
                out["G"] = np.sum((3 - 0.5 * x[:, 1:n-1]) * x[:, 1:n-1] - x[:, 0:n-2] - 2 * x[:,2:n] + 1, axis=1).reshape(x.shape[0],1)

            self.c_cache += out["G"].tolist()


    def write_cache(problem, filename):
        with open(filename, "w") as f:
            f.write(str(problem.n_var) + " " + str(problem.n_obj) + "\n")
            for elt_x, elt_f, elt_g in zip(problem.x_cache, problem.f_cache, problem.c_cache):
                for x in elt_x:
                    f.write(str(x) + " ")
                for y in elt_f:
                    f.write(str(y) + " ")
                for z in elt_g:
                    f.write(str(z) + " ")
                f.write("\n")

    # The folder where you want to store your caches
    foldername="./"
    if dict_problems[name][1] <= 2:
        return 0

    feasible_pbs = [
                    "CL1_1",
                    "CL1_2",
                    "CL1_5",
                    "CL1_6",
                    "DPAM1_4",
                    "DTLZ1_1",
                    "DTLZ1_3",
                    "DTLZ1_4",
                    "DTLZ2_1",
                    "DTLZ2_3",
                    "DTLZ2_4",
                    "DTLZ3_1",
                    "DTLZ3_3",
                    "DTLZ3_4",
                    "DTLZ4_1",
                    "DTLZ4_3",
                    "DTLZ4_4",
                    "DTLZ5_1",
                    "DTLZ5_3",
                    "DTLZ5_4",
                    "DTLZ6_1",
                    "DTLZ6_3",
                    "DTLZ6_4",
                    "FES1_1",
                    "FES1_3",
                    "FES1_4",
                    "FES2_1",
                    "FES2_3",
                    "FES2_4",
                    "FES3_1",
                    "FES3_3",
                    "FES3_4",
                    "I1_1",
                    "I1_3",
                    "I1_4",
                    "I2_1",
                    "I2_3",
                    "I2_4",
                    "I3_1",
                    "I3_3",
                    "I3_4",
                    "I4_1",
                    "I4_3",
                    "I4_4",
                    "I5_1",
                    "I5_3",
                    "I5_4",
                    "Kursawe_1",
                    "Kursawe_2",
                    "Kursawe_3",
                    "Kursawe_4",
                    "Kursawe_5",
                    "Kursawe_6",
                    "L1ZDT4_1",
                    "L1ZDT4_2",
                    "L1ZDT4_3",
                    "L1ZDT4_4",
                    "L1ZDT4_5",
                    "L1ZDT4_6",
                    "L2ZDT1_1",
                    "L2ZDT1_3",
                    "L2ZDT1_4",
                    "L2ZDT2_1",
                    "L2ZDT2_3",
                    "L2ZDT2_4",
                    "L2ZDT3_1",
                    "L2ZDT3_3",
                    "L2ZDT3_4",
                    "L2ZDT4_1",
                    "L2ZDT4_3",
                    "L2ZDT4_4",
                    "L2ZDT6_1",
                    "L2ZDT6_3",
                    "L2ZDT6_4",
                    "L3ZDT1_1",
                    "L3ZDT1_3",
                    "L3ZDT1_4",
                    "L3ZDT2_1",
                    "L3ZDT2_3",
                    "L3ZDT2_4",
                    "L3ZDT3_1",
                    "L3ZDT3_3",
                    "L3ZDT3_4",
                    "L3ZDT4_1",
                    "L3ZDT4_3",
                    "L3ZDT4_4",
                    "L3ZDT6_1",
                    "L3ZDT6_3",
                    "L3ZDT6_4",
                    "MOP2_1",
                    "MOP2_2",
                    "MOP2_3",
                    "MOP2_4",
                    "MOP2_5",
                    "MOP2_6",
                    "MOP4_1",
                    "MOP4_2",
                    "MOP4_3",
                    "MOP4_4",
                    "MOP4_5",
                    "MOP4_6",
                    "OKA2_1",
                    "OKA2_2",
                    "OKA2_3",
                    "OKA2_4",
                    "OKA2_5",
                    "OKA2_6",
                    "QV1_1",
                    "QV1_2",
                    "QV1_3",
                    "QV1_4",
                    "QV1_5",
                    "QV1_6",
                    "SK2_1",
                    "SK2_2",
                    "SK2_3",
                    "SK2_4",
                    "SK2_5",
                    "SK2_6",
                    "TKLY1_1",
                    "TKLY1_3",
                    "TKLY1_4",
                    "TKLY1_5",
                    "TKLY1_6",
                    "WFG1_1",
                    "WFG1_2",
                    "WFG1_3",
                    "WFG1_4",
                    "WFG1_5",
                    "WFG1_6",
                    "WFG2_1",
                    "WFG2_2",
                    "WFG2_3",
                    "WFG2_4",
                    "WFG2_5",
                    "WFG2_6",
                    "WFG3_1",
                    "WFG3_2",
                    "WFG3_3",
                    "WFG3_4",
                    "WFG3_5",
                    "WFG3_6",
                    "WFG4_1",
                    "WFG4_2",
                    "WFG4_3",
                    "WFG4_4",
                    "WFG4_5",
                    "WFG4_6",
                    "WFG5_1",
                    "WFG5_2",
                    "WFG5_3",
                    "WFG5_4",
                    "WFG5_5",
                    "WFG5_6",
                    "WFG6_1",
                    "WFG6_2",
                    "WFG6_3",
                    "WFG6_4",
                    "WFG6_5",
                    "WFG6_6",
                    "WFG7_1",
                    "WFG7_2",
                    "WFG7_3",
                    "WFG7_4",
                    "WFG7_5",
                    "WFG7_6",
                    "WFG8_1",
                    "WFG8_2",
                    "WFG8_3",
                    "WFG8_4",
                    "WFG8_5",
                    "WFG8_6",
                    "WFG9_1",
                    "WFG9_2",
                    "WFG9_3",
                    "WFG9_4",
                    "WFG9_5",
                    "WFG9_6",
                    "ZDT1_1",
                    "ZDT1_3",
                    "ZDT1_4",
                    "ZDT2_1",
                    "ZDT2_3",
                    "ZDT2_4",
                    "ZDT3_1",
                    "ZDT3_3",
                    "ZDT3_4",
                    "ZDT4_1",
                    "ZDT4_2",
                    "ZDT4_3",
                    "ZDT4_4",
                    "ZDT4_5",
                    "ZDT4_6",
                    "ZDT6_1",
                    "ZDT6_3",
                    "ZDT6_4",
                    "ZLT1_1",
                    "ZLT1_2",
                    "ZLT1_3",
                    "ZLT1_4",
                    "ZLT1_5",
                    "ZLT1_6",
                    "lovison5_1",
                    "lovison5_2",
                    "lovison5_3",
                    "lovison5_4",
                    "lovison5_5",
                    "lovison5_6",
                    "lovison6_1",
                    "lovison6_2",
                    "lovison6_3",
                    "lovison6_4",
                    "lovison6_5",
                    "lovison6_6"]

    print(name)
    for id_const in [1,2,3,4,5,6]:
        if name + "_" + str(id_const) in feasible_pbs:
            for seed in [1, 4734, 6652, 3507, 1121, 3500, 5816, 2006, 9622, 6117,
                         1571, 4117, 3758, 8045, 6554, 8521, 4889, 5893, 9123, 7238,
                         3089, 4197, 6489, 8, 7551, 7621, 4809, 8034, 3487, 3680]: # 30 runnings
                        #  5363, 7786, 104, 8171, 7505, 7075, 9656, 4569, 3121, 4187,
                        #  201, 951, 1093, 3147, 1978, 1475, 2896, 9814, 8803, 1084]:

                method = NSGA2(pop_size=100)

                problem = MXProblem(name, id_const)
                minimize(
                    problem, method, termination=("n_gen", 300), seed=seed, save_history=True, disp=False
                )
                write_cache(problem, foldername + name + "_" + str(id_const) + "_nsgaii_" + str(seed) + ".txt")
    print("Done!")
    return 0


if __name__ == '__main__':
    benchmark_names = ["BK1",
                       "CL1",
                       "Deb41",
                       "Deb512a",
                       "Deb512b",
                       "Deb512c",
                       "Deb513",
                       "Deb521a",
                       "Deb521b",
                       "Deb53",
                       "DG01",
                       "DPAM1",
                       "DTLZ1",
                       "DTLZ1n2",
                       "DTLZ2",
                       "DTLZ2n2",
                       "DTLZ3",
                       "DTLZ3n2",
                       "DTLZ4",
                       "DTLZ4n2",
                       "DTLZ5",
                       "DTLZ5n2",
                       "DTLZ6",
                       "DTLZ6n2",
                       "ex005",
                       "Far1",
                       "FES1",
                       "FES2",
                       "FES3",
                       "Fonseca",
                       "I1",
                       "I2",
                       "I3",
                       "I4",
                       "I5",
                       "IKK1",
                       "IM1",
                       "Jin1",
                       "Jin2",
                       "Jin3",
                       "Jin4",
                       "Kursawe",
                       "L1ZDT4",
                       "L2ZDT1",
                       "L2ZDT2",
                       "L2ZDT3",
                       "L2ZDT4",
                       "L2ZDT6",
                       "L3ZDT1",
                       "L3ZDT2",
                       "L3ZDT3",
                       "L3ZDT4",
                       "L3ZDT6",
                       "LE1",
                       "lovison1",
                       "lovison2",
                       "lovison3",
                       "lovison4",
                       "lovison5",
                       "lovison6",
                       "LRS1",
                       "MHHM1",
                       "MHHM2",
                       "MLF1",
                       "MLF2",
                       "MOP1",
                       "MOP2",
                       "MOP3",
                       "MOP4",
                       "MOP5",
                       "MOP6",
                       "MOP7",
                       "OKA1",
                       "OKA2",
                       "QV1",
                       "Sch1",
                       "SK1",
                       "SK2",
                       "SP1",
                       "SSFYY1",
                       "SSFYY2",
                       "TKLY1",
                       "VFM1",
                       "VU1",
                       "VU2",
                       "WFG1",
                       "WFG2",
                       "WFG3",
                       "WFG4",
                       "WFG5",
                       "WFG6",
                       "WFG7",
                       "WFG8",
                       "WFG9",
                       "ZDT1",
                       "ZDT2",
                       "ZDT3",
                       "ZDT4",
                       "ZDT6",
                       "ZLT1"]

    for name in tqdm(benchmark_names):
        run_nsgaii_solver_pb(name)
    #  Parallel(n_jobs=3)(delayed(run_nsgaii_solver_pb)(name) for name in tqdm(benchmark_names))
