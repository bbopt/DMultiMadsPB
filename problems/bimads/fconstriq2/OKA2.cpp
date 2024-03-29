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
/*                  OKA2                */
/*----------------------------------------*/
class OKA2 : public Multi_Obj_Evaluator
{

public:
    OKA2(const NOMAD::Parameters &p) : NOMAD::Multi_Obj_Evaluator(p) {}

    ~OKA2(void) {}

    bool eval_x(NOMAD::Eval_Point &x,
                const NOMAD::Double &h_max,
                bool &count_eval) const
    {
        int n = 3;

        NOMAD::Double f1 = x[0];
        NOMAD::Double f2 = 1 - (x[0] + PI) * (x[0] + PI) / (4 * PI * PI) 
            + pow( abs(x[1].value() - 5 * cos(x[0].value())), 1.0/3 )
            + pow( abs(x[2].value() - 5 * sin(x[0].value())), 1.0/3 );

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
        int n = 3;
        int m = 2;
        int l = n-2;

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

        NOMAD::Point lb(n, -5);
        lb[0] = -PI;
        NOMAD::Point ub(n, 5);
        ub[0] = PI;

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

        // p.set_MULTI_NB_MADS_RUNS(30);

        // p.set_TMP_DIR ("/tmp");      // directory for
        // temporary files

    int seed = 0;
    p.set_SEED(seed);

        p.set_HISTORY_FILE("OKA2_2_bimads_"+to_string(seed)+".txt");
        p.set_STATS_FILE("test_OKA2.txt", "BBE OBJ");

        // // disable models
        // p.set_DISABLE_MODELS();
        //
        // // disable NM search
        // p.set_NM_SEARCH(false);

        // parameters validation:
        p.check();

        // custom evaluator creation:
        OKA2 ev(p);

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
