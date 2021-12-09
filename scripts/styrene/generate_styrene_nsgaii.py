import numpy as np
import subprocess
from datetime import datetime
from tqdm import tqdm

from pymoo.model.problem import Problem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.optimize import minimize

class StyreneProblem(Problem):
    def __init__(self):
        super().__init__(n_var = 8,
                         n_obj = 3,
                         n_constr = 9,
                         xl=np.zeros(8),
                         xu=100*np.ones(8),
                         elementwise_evaluation=True)
        self.x_cache = []
        self.f_cache = []
        self.c_cache = []
        self.seed_tag = None

    # x represents a population array
    def _evaluate(self, x, out, *args, **kwargs):

        self.x_cache += [x.tolist()]

        # Call the blackbox on each point
        tag = datetime.now()
        tag = tag.strftime("%d%m%Y%H%M%S")
        input_coordinates_filename = "styrene_tmp_x_" + str(self.seed_tag) + "_" + tag + ".txt"

        # Real entries
        inputs_styrene = np.copy(x)

        with open(input_coordinates_filename, "w") as f:
            for elt in inputs_styrene.tolist():
                f.write(str(elt) + "\n")

        # Launch simulation
        result = subprocess.run(["./truth.exe", input_coordinates_filename],
                                capture_output=True,
                                encoding="UTF-8")
        outputs = None
        if "ERROR" in result.stdout:
            outputs = [np.inf for i in range(12)]
        else:
            # Collect outputs
            outputs = [float(elt) for elt in result.stdout.rsplit()]

        # Delete temporary file
        subprocess.run(["rm", input_coordinates_filename])

        # Convert to objective and constraints functions
        # outputs[0:4] boolean constraints: set to 1 if simulation fails, 0 otherwise.
        # outputs[4] : minimal purity of produced styrene (f2)
        # outputs[5] : minimal purity of produced benzene
        # outputs[6] : overall ethylbenzene conversion into styrene (f3)
        # outputs[7:11] : Four constraints relating to payout time, cashflow, investment and annual costs
        # outputs[11] : net present value of the project (f1)
        self.f_cache += [[outputs[i] for i in [11,4,6]]]
        self.c_cache += [[outputs[i] for i in [0,1,2,3,5,7,8,9,10]]]

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

def run_styrene_nsgaii(seed):
    method = NSGA2(pop_size=100)

    problem = StyreneProblem()
    problem.seed_tag=seed
    minimize(
        problem, method, termination=("n_gen", 200), seed=seed, save_history=True, disp=False
    )
    write_cache(problem, "styrene" + "_nsgaii_" + str(seed) + ".txt")
    return

if __name__ == '__main__':
    for seed in [1, 4734, 6652, 3507, 1121, 3500, 5816, 2006, 9622, 6117]:
        print("problem STYRENE ", "seed = " + str(seed))
        run_styrene_nsgaii(seed)
        print("done")
