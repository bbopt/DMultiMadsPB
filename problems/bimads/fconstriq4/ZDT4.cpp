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
/*                  ZDT4                */
/*----------------------------------------*/
class ZDT4 : public Multi_Obj_Evaluator {

    public:
        ZDT4  ( const NOMAD::Parameters & p ) :
            NOMAD::Multi_Obj_Evaluator ( p ) {}

        ~ZDT4 ( void ) {}

        bool eval_x ( NOMAD::Eval_Point   & x,
                const NOMAD::Double & h_max,
                bool                & count_eval   ) const
        {
            int n = 10;
            NOMAD::Double g = 0;
            for (int i = 1; i < n; ++i) {
                g += x[i] * x[i] - 10 * cos(4 * PI * x[i].value());
            }
            g += 1 + 10 * (n - 1);

            NOMAD::Double f1 = x[0];
            NOMAD::Double h = 1 - sqrt(f1.value() / g.value());

            x.set_bb_output  ( 0 , f1); // objective 1
            x.set_bb_output  ( 1 , g * h ); // objective 2

            int l = n - 1;
            for (int j =0; j < l; ++j) {
                NOMAD::Double c = x[j]*x[j] + x[j+1]*x[j+1] + x[j]*x[j+1] - 1;
                x.set_bb_output(j+2, c); // constraints
            }

            count_eval = true; // count a black-box evaluation

            return true;       // the evaluation succeeded
        }

};

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

    // display:
    NOMAD::Display out ( std::cout );
    out.precision ( NOMAD::DISPLAY_PRECISION_STD );

    try {
        // NOMAD initializations:
        NOMAD::begin(argc, argv);

        // parameters creation:
        NOMAD::Parameters p(out);

        // dimensions of the black box
        int n = 10;
        int m = 2;
        int l = n-1;

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

        NOMAD::Point lb(n, -5.0);
        lb[0] = 0.0;
        NOMAD::Point ub(n, 5.0);
        ub[0] = 1.0;

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

        p.set_HISTORY_FILE("ZDT4_4_bimads_"+to_string(seed)+".txt");
        p.set_STATS_FILE("test_ZDT4.txt", "BBE OBJ");


        // // disable models
        // p.set_DISABLE_MODELS();
        //
        // // disable NM search
        // p.set_NM_SEARCH(false);

        // parameters validation:
        p.check();

        // custom evaluator creation:
        ZDT4 ev(p);

        // algorithm creation and execution:
        Mads mads(p, &ev);
        mads.multi_run();
    }
    catch ( exception & e ) {
        cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }

    Slave::stop_slaves ( out );
    end();

    return EXIT_SUCCESS;
}

