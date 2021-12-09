import numpy as np
import subprocess
from datetime import datetime
from joblib import Parallel, delayed
from tqdm import tqdm

from pymoo.model.problem import Problem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.optimize import minimize

class Solar8Problem(Problem):
    def __init__(self):
        super().__init__(n_var = 11,
                         n_obj = 2,
                         n_constr = 9,
                         xl=np.array([1.0, 1.0, 20.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.01, 0.005, 0.0060]),
                         xu=np.array([40.0, 40.0, 250.0, 30.0, 30.0, 89.0, 20.0, 20.0, 5.00, 0.100, 0.1000]),
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
        input_coordinates_filename = "solar8_tmp_x_" + tag + ".txt"

        # Real entries
        inputs_solar8 = np.zeros(13)
        inputs_solar8[0:5] = x[0:5]
        inputs_solar8[6:9] = x[5:8]
        inputs_solar8[10:13] = x[8:11]

        # Integer entries: as this version does not support integer entries, they are
        # fixed to the x0 proposed
        inputs_solar8[5] = 2650 # Maximum number of heliostats
        inputs_solar8[9] = 36 # Receiver number of tubes

        with open(input_coordinates_filename, "w") as f:
            for elt in inputs_solar8.tolist():
                f.write(str(elt) + "\n")

        # Launch simulation
        result = subprocess.run(["./solar_bb.exe", "8", input_coordinates_filename],
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

def run_solar8_nsgaii(seed):
    method = NSGA2(pop_size=100)

    problem = Solar8Problem()
    minimize(
        problem, method, termination=("n_gen", 50), seed=seed, save_history=True, disp=False
    )
    write_cache(problem, "solar8" + "_nsgaii_" + str(seed) + ".txt")
    return

if __name__ == '__main__':
    for seed in [1, 4734, 6652, 3507, 1121, 3500, 5816, 2006, 9622, 6117]:
        print("problem Solar8 ", "seed = " + str(seed))
        run_solar8_nsgaii(seed)
        print("done")
