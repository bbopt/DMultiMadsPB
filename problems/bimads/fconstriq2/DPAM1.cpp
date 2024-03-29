#include <iostream>
#include "nomad.hpp"
#include "math.h"
using namespace NOMAD;

const double PI = 3.141592653589793238463;

/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
using namespace std;

// using namespace NOMAD; avoids putting NOMAD:: everywhere

/*----------------------------------------*/
/*                  DPAM1                */
/*----------------------------------------*/
class DPAM1 : public Multi_Obj_Evaluator
{

public:
  DPAM1(const NOMAD::Parameters &p) : NOMAD::Multi_Obj_Evaluator(p) {}

  ~DPAM1(void) {}

  bool eval_x(NOMAD::Eval_Point &x,
              const NOMAD::Double &h_max,
              bool &count_eval) const
  {
    int n = 10;

    // uniform matrix (-1, 1) 10 X 10
    double A[10][10] = {
        {0.218418, -0.620254, 0.843784, 0.914311, -0.788548, 0.428212,
         0.103064, -0.47373, -0.300792, -0.185507},
        {0.330423, 0.151614, 0.884043, -0.272951, -0.993822, 0.511197,
         -0.0997948, -0.659756, 0.575496, 0.675617},
        {0.180332, -0.593814, -0.492722, 0.0646786, -0.666503, -0.945716,
         -0.334582, 0.611894, 0.281032, 0.508749},
        {-0.0265389, -0.920133, 0.308861, -0.0437502, -0.374203, 0.207359,
         -0.219433, 0.914104, 0.184408, 0.520599},
        {-0.88565, -0.375906, -0.708948, -0.37902, 0.576578, 0.0194674,
         -0.470262, 0.572576, 0.351245, -0.480477},
        {0.238261, -0.1596, -0.827302, 0.669248, 0.494475, 0.691715,
         -0.198585, 0.0492812, 0.959669, 0.884086},
        {-0.218632, -0.865161, -0.715997, 0.220772, 0.692356, 0.646453,
         -0.401724, 0.615443, -0.0601957, -0.748176},
        {-0.207987, -0.865931, 0.613732, -0.525712, -0.995728, 0.389633,
         -0.064173, 0.662131, -0.707048, -0.340423},
        {0.60624, 0.0951648, -0.160446, -0.394585, -0.167581, 0.0679849,
         0.449799, 0.733505, -0.00918638, 0.00446808},
        {0.404396, 0.449996, 0.162711, 0.294454, -0.563345, -0.114993,
         0.549589, -0.775141, 0.677726, 0.610715}};

    // transformation
    NOMAD::Point y(n, 0);
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
      {
        y[i] += A[i][j] * x[j];
      }
    }

    NOMAD::Double g = 1 + 10 * (n - 1);
    for (int i = 1; i < n; ++i)
    {
      g += y[i] * y[i] - 10 * cos(4 * PI * y[i].value());
    }
    NOMAD::Double f1 = y[0];
    NOMAD::Double f2 = g.value() * exp(-y[0].value() / g.value());

    x.set_bb_output(0, f1); // objective 1
    x.set_bb_output(1, f2); // objective 2

    int l = n - 2;
    for (int j =0; j < l; ++j) {
        NOMAD::Double c = (3 - 2 * x[j+1]) * x[j+1] - x[j] - 2 * x[j+2] + 2.5;
        x.set_bb_output(j+2, c); // constraints
    }

    count_eval = true; // count a black-box evaluation

    return true; // the evaluation succeeded
  }
};

/*------------------------------------------*/
/*            NOMAD main function           */
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
    int n = 10;
    int m = 2;
    int l = n - 2;

    p.set_DIMENSION(n); // number of variables

    vector<NOMAD::bb_output_type> bbot(m+l); // definition of output types
    for (int i = 0; i < m; ++i)
    {
      bbot[i] = OBJ;
    }
    for (int i = m; i < m+l; ++i)
    {
      bbot[i] = PB;
    }
    p.set_BB_OUTPUT_TYPE(bbot);

    //    p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.

    NOMAD::Point lb(n, -0.3);
    NOMAD::Point ub(n, 0.3);

    p.set_LOWER_BOUND(lb); // all var. >= 0.0
    p.set_UPPER_BOUND(ub); //all var <= 1

    for (int j = 0; j < n; ++j)
    {
      NOMAD::Point x0(n, 0);
      for (int i = 0; i < n; ++i)
      {
        x0[i] = lb[i] + j * (ub[i] - lb[i]) / (n - 1);
      }
      p.set_X0(x0);
    }
    //p.set_X0(x0); // starting point

    //    p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.

    // p.set_LH_SEARCH(20, -1);

    //p.set_DISPLAY_DEGREE(NOMAD::MINIMAL_DISPLAY);
    p.set_DISPLAY_STATS("obj");

    p.set_MULTI_OVERALL_BB_EVAL(30000); // the algorithm terminates after
                                        // 100 black-box evaluations

    //p.set_MULTI_NB_MADS_RUNS(30);

    // p.set_TMP_DIR ("/tmp");      // directory for
    // temporary files

    int seed = 0;
    p.set_SEED(seed);

    p.set_HISTORY_FILE("DPAM1_2_bimads_"+to_string(seed)+".txt");
    p.set_STATS_FILE("test_DPAM1.txt", "BBE OBJ");

    // // disable models
    // p.set_DISABLE_MODELS();
    //
    // // disable NM search
    // p.set_NM_SEARCH(false);

    // parameters validation:
    p.check();

    // custom evaluator creation:
    DPAM1 ev(p);

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
