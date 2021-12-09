import numpy as np
import subprocess
from datetime import datetime
from tqdm import tqdm

from pymoo.model.problem import Problem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.optimize import minimize

class Solar9Problem(Problem):
    def __init__(self):
        super().__init__(n_var = 22,
                         n_obj = 2,
                         n_constr = 17,
                         xl=np.array([1.0, 1.0, 20.0, 1.0, 1.0, 1.0, 0.0, 1.0, 793.0, 1.0, 1.0, 0.01, 0.01, 495.0, 0.01, 0.0050, 0.006, 0.007, 0.5, 0.0050, 0.006, 0.15]),
                         xu=np.array([40.0, 40.0, 250.0, 30.0, 30.0, 89.0, 20.0, 20.0, 995.0, 50.0, 30.0, 5.00, 5.00, 650.0, 5.00, 0.1000, 0.100, 0.200, 10.0, 0.1000, 0.100, 0.40]),
                         elementwise_evaluation=True)
        self.x_cache = []
        self.f_cache = []
        self.c_cache = []

    # x represents a population array
    def _evaluate(self, x, out, *args, **kwargs):

        self.x_cache += [x.tolist()]

        # Call the blackbox on each point
        tag = datetime.now()
        tag = tag.strftime("%d%m%Y%H%M%S")
        input_coordinates_filename = "solar9_tmp_x_" + tag + ".txt"

        # Real entries
        inputs_solar9 = np.zeros(29)
        inputs_solar9[0:5] = x[0:5]
        inputs_solar9[6:15] = x[5:14]
        inputs_solar9[16:24] = x[14:22]

        # Integer entries: as this version does not support integer entries, they are
        # fixed to the x0 proposed
        inputs_solar9[5] = 1000 # Maximum number of heliostats
        inputs_solar9[15] = 500 # Receiver number of tubes
        inputs_solar9[24] = 3 # Exchanger number of tubes
        inputs_solar9[25] = 12000 # Exchanger number of tubes
        inputs_solar9[26] = 1 # Exchanger number of shells
        inputs_solar9[27] = 2 # Exchanger number of passes per shell
        inputs_solar9[28] = 2 # Type of turbine

        with open(input_coordinates_filename, "w") as f:
            for elt in inputs_solar9.tolist():
                f.write(str(elt) + "\n")

        # Launch simulation
        result = subprocess.run(["./solar_bb.exe", "9", input_coordinates_filename],
                                capture_output=True,
                                encoding="UTF-8")
        # Collect outputs
        outputs = [float(elt) for elt in result.stdout.rsplit()]
        #  print(outputs)

        # Delete temporary file
        subprocess.run(["rm", input_coordinates_filename])

        # Convert to objective and constraints functions
        self.f_cache += [[outputs[i] for i in range(2)]]
        self.c_cache += [[outputs[i] for i in range(2, len(outputs))]]

        out["F"] = [outputs[i] for i in range(2)]
        out["G"] = [outputs[i] for i in range(2, len(outputs))]


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

def run_solar9_nsgaii(seed):
    method = NSGA2(pop_size=100)

    problem = Solar9Problem()
    minimize(
        problem, method, termination=("n_gen", 50), seed=seed, save_history=True, disp=False
    )
    write_cache(problem, "solar9" + "_nsgaii_" + str(seed) + ".txt")
    return

if __name__ == '__main__':
    for seed in [1, 4734, 6652, 3507, 1121, 3500, 5816, 2006, 9622, 6117]:
        print("problem Solar9 ", "seed = " + str(seed))
        run_solar9_nsgaii(seed)
        print("done")
