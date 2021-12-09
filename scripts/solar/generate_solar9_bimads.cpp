#include <iostream>
#include "nomad.hpp"
#include "math.h"
using namespace NOMAD;
using namespace std;

/*----------------------------------------*/
/*               Solar                  */
/*----------------------------------------*/
class SolarEvaluator : public NOMAD::Multi_Obj_Evaluator
{

public:
  SolarEvaluator(const NOMAD::Parameters &p) : NOMAD::Multi_Obj_Evaluator(p)
    {}

  ~SolarEvaluator(void) {}
};

/*------------------------------------------*/
/*             Main program                 */
/*------------------------------------------*/
int main(int argc, char **argv)
{

  // display:
  NOMAD::Display out(std::cout);
  out.precision(NOMAD::DISPLAY_PRECISION_STD);

  try
  {

    // NOMAD initializations:
    NOMAD::begin(argc, argv);

    // parameters creation:
    NOMAD::Parameters p(out);

    // dimensions of the blackbox
    int n = 29;
    int m = 17;

    p.set_DIMENSION(n); // number of variables

    vector<NOMAD::bb_output_type> bbot(m); // definition of output types
    // It has two objectives
    for (int i = 0; i < 2; ++i) {
        bbot[i] = OBJ;
    }
    for (int i = 2; i < 17; ++i) {
        bbot[i] = PB;
    }

    // And 9 inequality constraints
    p.set_BB_OUTPUT_TYPE(bbot);

    // p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.

    // parameters
    NOMAD::Point lb(n);
    lb[0] = 1.0;
    lb[1] = 1.0;
    lb[2] = 20.0;
    lb[3] = 1.0;
    lb[4] = 1.0;
    lb[5] = 1.0;
    lb[6] = 1.0;
    lb[7] = 0.0;
    lb[8] = 1.0;
    lb[9] = 793.0;
    lb[10] = 1.0;
    lb[11] = 1.0;
    lb[12] = 0.01;
    lb[13] = 0.01;
    lb[14] = 495.0;
    lb[15] = 1.0;
    lb[16] = 0.01;
    lb[17] = 0.0050;
    lb[18] = 0.006;
    lb[19] = 0.007;
    lb[20] = 0.5;
    lb[21] = 0.0050;
    lb[22] = 0.006;
    lb[23] = 0.15;
    lb[24] = 2.0;
    lb[25] = 1.0;
    lb[26] = 1.0;
    lb[27] = 1.0;
    lb[28] = 1.0;

    NOMAD::Point ub(n);
    ub[0] = 40.0;
    ub[1] = 40.0;
    ub[2] = 250.0;
    ub[3] = 30.0;
    ub[4] = 30.0;
    ub[6] = 89.0;
    ub[7] = 20.0;
    ub[8] = 20.0;
    ub[9] = 995.00;
    ub[10] = 50.0;
    ub[11] = 30.0;
    ub[12] = 5.0;
    ub[13] = 5.0;
    ub[14] = 650.0;
    ub[16] = 5.0;
    ub[17] = 0.1000;
    ub[18] = 0.100;
    ub[19] = 0.200;
    ub[20] = 10.0;
    ub[21] = 0.1000;
    ub[22] = 0.100;
    ub[23] = 0.40;
    ub[26] = 10.0;
    ub[27] = 9.0;
    ub[28] = 8.0;

    p.set_LOWER_BOUND(lb); // all var. >= lb
    p.set_UPPER_BOUND(ub); //all var <= ub

    // bb input types
    p.set_BB_INPUT_TYPE(5, INTEGER);
    p.set_BB_INPUT_TYPE(15, INTEGER);

    p.set_BB_INPUT_TYPE(24, INTEGER);
    p.set_BB_INPUT_TYPE(25, INTEGER);
    p.set_BB_INPUT_TYPE(26, INTEGER);
    p.set_BB_INPUT_TYPE(27, INTEGER);
    p.set_BB_INPUT_TYPE(28, INTEGER);

    // starting point
    p.set_X0("x0_solar9.txt");

    // p.set_FIXED_VARIABLE(5);
    // p.set_FIXED_VARIABLE(15);
    //
    // p.set_FIXED_VARIABLE(24);
    // p.set_FIXED_VARIABLE(25);
    // p.set_FIXED_VARIABLE(26);
    // p.set_FIXED_VARIABLE(27);
    // p.set_FIXED_VARIABLE(28);

    //p.set_X0(x0); // starting poin

    //p.set_LH_SEARCH(20,-1);

    //p.set_DISPLAY_DEGREE(NOMAD::MINIMAL_DISPLAY);
    p.set_DISPLAY_STATS("obj");

    p.set_MULTI_OVERALL_BB_EVAL(20000); //(5000); // the algorithm terminates after
    // 100 black-box evaluations
    
    // p.set_MULTI_NB_MADS_RUNS(30);

    //p.set_MULTI_OVERALL_BB_EVAL (100); // the algorithm terminates after
    // 100 black-box evaluations

    // p.set_TMP_DIR ("/tmp");      // directory for
    // temporary files

    int seed = 0;
    // int seed = 1;
    p.set_SEED(seed);
    
    p.set_HISTORY_FILE("SOLAR9_bimads_"+to_string(seed)+".txt");
    p.set_STATS_FILE("test_SOLAR9.txt", "BBE OBJ");

    // // disable models
    // p.set_DISABLE_MODELS();
    //
    // // disable Nelder Mead search
    // p.set_NM_SEARCH(false);

    p.set_BB_EXE("$./solar_bb.exe $9");

    // parameters validation:
    p.check();

    SolarEvaluator ev(p);

    // algorithm creation and execution:
    Mads mads(p, &ev);
    mads.multi_run();
  }
  catch (exception &e)
  {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  Slave::stop_slaves(out);
  end();

  return EXIT_SUCCESS;
}
