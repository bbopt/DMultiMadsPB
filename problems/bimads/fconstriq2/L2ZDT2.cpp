#include <iostream>
#include "nomad.hpp"
#include "math.h"
using namespace NOMAD;

/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
using namespace std;

// using namespace NOMAD; avoids putting NOMAD:: everywhere

/*----------------------------------------*/
/*                  L2ZDT2                */
/*----------------------------------------*/
class L2ZDT2 : public Multi_Obj_Evaluator
{

    public:
        L2ZDT2(const NOMAD::Parameters &p) : NOMAD::Multi_Obj_Evaluator(p) {}

        ~L2ZDT2(void) {}

        bool eval_x(NOMAD::Eval_Point &x,
                const NOMAD::Double &h_max,
                bool &count_eval) const
        {
            int n = 30;

            // random matrix M 30 * 30
            double M[30][30] = {
                {0.218418, -0.620254, 0.843784, 0.914311, -0.788548, 0.428212, 0.103064, -0.47373, -0.300792, -0.185507, 0.330423, 0.151614, 0.884043, -0.272951, -0.993822, 0.511197, -0.0997948, -0.659756, 0.575496, 0.675617, 0.180332, -0.593814, -0.492722, 0.0646786, -0.666503, -0.945716, -0.334582, 0.611894, 0.281032, 0.508749},
                {-0.0265389, -0.920133, 0.308861, -0.0437502, -0.374203, 0.207359, -0.219433, 0.914104, 0.184408, 0.520599, -0.88565, -0.375906, -0.708948, -0.37902, 0.576578, 0.0194674, -0.470262, 0.572576, 0.351245, -0.480477, 0.238261, -0.1596, -0.827302, 0.669248, 0.494475, 0.691715, -0.198585, 0.0492812, 0.959669, 0.884086},
                {-0.218632, -0.865161, -0.715997, 0.220772, 0.692356, 0.646453, -0.401724, 0.615443, -0.0601957, -0.748176, -0.207987, -0.865931, 0.613732, -0.525712, -0.995728, 0.389633, -0.064173, 0.662131, -0.707048, -0.340423, 0.60624, 0.0951648, -0.160446, -0.394585, -0.167581, 0.0679849, 0.449799, 0.733505, -0.00918638, 0.00446808},
                {0.404396, 0.449996, 0.162711, 0.294454, -0.563345, -0.114993, 0.549589, -0.775141, 0.677726, 0.610715, 0.0850755, 0.0419388, -0.323614, -0.973719, -0.680238, -0.270873, -0.209617, 0.968436, 0.908798, 0.975851, -0.994918, -0.0621977, 0.628171, 0.761228, 0.34372, -0.792042, -0.144765, -0.965748, 0.0133606, -0.0260565},
                {-0.742377, 0.426391, 0.408202, 0.633885, -0.0351053, -0.723444, -0.577654, 0.0276004, 0.0712472, -0.622791, 0.155451, 0.442717, -0.792786, 0.925785, 0.670266, -0.865566, -0.638281, 0.333094, 0.477628, 0.47261, 0.23151, 0.82132, -0.589803, 0.796275, 0.57713, 0.101149, 0.970191, 0.532821, 0.814769, -0.0687269},
                {0.712758, -0.191812, -0.390938, 0.952828, 0.921519, 0.923094, 0.93011, -0.945394, -0.0934027, 0.964123, -0.795609, -0.289563, 0.614236, -0.670585, 0.466877, 0.144597, -0.206416, 0.6937, -0.967958, -0.0951247, -0.942473, -0.610767, -0.655472, -0.0960986, -0.302779, -0.734976, -0.342188, -0.315861, -0.912834, 0.24499},
                {0.0969326, 0.089775, -0.241157, 0.0835558, -0.420236, -0.686633, -0.711276, -0.00325377, 0.435196, -0.710002, 0.00283691, -0.168757, -0.134045, -0.655235, 0.172361, 0.998291, 0.376291, -0.962215, -0.363174, -0.88777, -0.519929, -0.560554, -0.984415, 0.601529, -0.984103, -0.228237, -0.578066, 0.307023, 0.606123, 0.959635},
                {0.00225943, 0.0101814, 0.441456, 0.0633629, 0.406631, -0.0100638, -0.177972, -0.491075, 0.537035, -0.924987, -0.699424, 0.742285, 0.0181443, 0.718971, -0.0308272, 0.086931, 0.524476, 0.956457, 0.143024, 0.616481, 0.217909, -0.128427, -0.262427, -0.938208, -0.52479, 0.12919, 0.721925, 0.766492, 0.470845, -0.0976466},
                {0.507807, 0.804148, 0.963269, 0.357128, -0.832565, -0.312441, 0.327779, 0.184745, 0.246139, -0.936814, -0.931734, -0.0327827, 0.319293, 0.044473, -0.641645, 0.596118, -0.293934, -0.63373, 0.409658, 0.759892, -0.257078, 0.939616, -0.227661, 0.115754, 0.10964, -0.240557, 0.66842, 0.855535, -0.451536, 0.264961},
                {-0.61366, -0.204783, -0.842476, -0.249524, -0.0985226, 0.0671501, -0.527707, -0.509489, -0.883254, 0.14851, -0.906465, 0.496238, -0.853211, -0.779234, -0.979515, 0.827175, 0.228969, -0.402829, -0.970118, 0.762559, 0.506495, 0.460303, 0.897304, 0.686003, 0.739986, 0.15731, 0.281697, -0.922955, -0.780824, 0.449716},
                {0.125225, 0.487649, 0.147046, 0.679639, 0.593707, -0.311828, -0.797099, -0.35815, 0.95808, 0.907244, 0.772426, 0.720574, -0.873217, 0.371431, -0.826029, 0.942716, 0.70609, -0.658158, -0.782185, -0.806743, -0.627986, -0.405551, -0.258495, -0.796524, 0.222498, 0.087545, -0.0917108, -0.62542, -0.110256, 0.0417576},
                {0.24476, 0.941339, -0.613783, 0.402772, 0.300775, -0.820314, -0.894233, -0.405896, 0.0735439, 0.486645, -0.394355, 0.125097, -0.316386, -0.701215, -0.845742, 0.2065, -0.413743, 0.406725, -0.423813, -0.941255, -0.558804, 0.312326, 0.345314, 0.319143, -0.644653, -0.0408415, 0.176461, 0.740113, 0.470737, -0.914927},
                {-0.591523, -0.606614, -0.181873, 0.692975, 0.50208, -0.536704, 0.359652, 0.839082, 0.56817, -0.0776788, -0.00332785, 0.459538, -0.518313, -0.270738, -0.629958, -0.755084, -0.721573, 0.431107, -0.221877, 0.32543, 0.163743, 0.0759916, 0.695064, -0.656856, 0.074163, 0.264319, -0.73174, 0.731548, -0.489341, 0.678946},
                {0.0271269, 0.804879, -0.402973, 0.800373, 0.760082, -0.878364, 0.176801, -0.548932, -0.225601, -0.164912, -0.208143, 0.7768, -0.542743, -0.156021, 0.671736, 0.878648, -0.419588, -0.0752896, 0.0299447, -0.494459, -0.72415, 0.35978, -0.32646, -0.96605, 0.0127605, 0.563174, -0.814853, -0.949609, -0.526794, -0.801902},
                {-0.753397, 0.617418, 0.689874, 0.983384, 0.668786, 0.0304653, -0.625221, -0.13318, 0.827343, -0.101358, -0.999522, -0.0525574, -0.458319, 0.587409, -0.334639, 0.0759643, 0.0255827, 0.128944, 0.17317, -0.284309, 0.287161, -0.550725, -0.433083, -0.242821, 0.878879, 0.691699, -0.660499, 0.389985, 0.599856, -0.711442},
                {-0.798697, -0.244945, -0.942649, 0.402856, -0.494672, 0.439941, -0.88216, 0.170196, 0.650734, -0.0982391, -0.468732, 0.342133, -0.838071, -0.832362, 0.658177, -0.565361, 0.149473, 0.69331, -0.491848, 0.74916, 0.526025, -0.155339, 0.0998096, 0.468761, 0.324649, 0.128488, 0.544144, -0.495222, 0.965229, -0.79314},
                {-0.545421, -0.500243, 0.154371, 0.170017, -0.259108, -0.868862, -0.50731, -0.848317, 0.835712, 0.616391, -0.442608, -0.158, 0.313451, 0.703748, -0.755984, -0.249443, 0.491564, 0.985068, 0.678644, 0.808324, 0.81975, -0.435823, -0.839855, 0.00282368, -0.569165, 0.0884339, -0.222144, 0.499412, -0.565198, 0.64824},
                {0.956914, -0.0620912, 0.634479, 0.928617, 0.464664, 0.377022, 0.63047, -0.198619, -0.576153, 0.565373, -0.524245, -0.187299, -0.614524, 0.429316, -0.491171, 0.399495, -0.333898, -0.646636, -0.0189709, -0.339605, -0.798791, 0.0494081, 0.367012, 0.852545, 0.43557, 0.150039, -0.0454542, 0.604861, -0.598288, -0.500696},
                {0.249008, 0.370711, -0.633174, -0.0121906, 0.42006, 0.169373, -0.975542, -0.0297852, 0.80481, 0.638317, -0.670967, 0.935792, -0.35605, 0.175773, 0.878601, -0.275168, -0.932517, -0.372497, -0.0732907, -0.185493, -0.357004, 0.314786, -0.229239, 0.530256, -0.51327, 0.44187, 0.940309, -0.240334, -0.0276121, 0.74383},
                {-0.630329, -0.763376, 0.62538, 0.818945, 0.891598, 0.680494, 0.471868, -0.769787, -0.878099, -0.973724, 0.354362, -0.1792, -0.225034, -0.44548, 0.598865, 0.544005, -0.478962, 0.327193, -0.525784, 0.903179, -0.899248, 0.156514, 0.154329, 0.499808, -0.836327, -0.802627, 0.378082, -0.112673, -0.47926, -0.3355},
                {-0.699445, 0.237731, -0.324597, -0.800406, -0.42585, -0.710739, -0.144068, -0.828545, -0.800912, 0.184654, -0.63675, -0.16696, 0.240427, -0.513443, 0.812664, 0.744943, 0.970612, 0.00172899, -0.726378, -0.0985012, 0.224232, 0.16495, 0.560077, -0.813112, 0.112894, -0.0955366, 0.0187107, 0.913887, 0.123076, 0.550338},
                {0.400334, -0.367816, 0.198455, -0.983183, 0.278976, 0.714817, 0.307911, 0.812861, -0.403497, -0.784382, -0.161823, -0.120835, 0.323172, 0.583739, 0.732924, -0.220603, -0.594121, 0.935093, -0.216736, 0.659318, -0.750417, -0.284773, -0.271496, 0.491731, -0.712174, -0.763681, 0.0781023, 0.951666, 0.734031, 0.826912},
                {0.57488, -0.361951, -0.0739728, 0.91438, -0.391653, 0.0193887, 0.412634, -0.169813, 0.471794, 0.660792, -0.350906, -0.612644, 0.347876, 0.112573, -0.501126, 0.456761, -0.109004, 0.289352, -0.566504, 0.585042, 0.584934, 0.923676, 0.895312, -0.161036, -0.995895, 0.0853141, -0.583368, -0.157612, 0.234119, 0.875043},
                {0.430805, 0.706102, 0.423887, 0.296828, -0.265607, 0.338806, -0.15829, 0.642516, 0.355126, 0.174447, -0.975015, 0.869905, -0.145285, -0.484002, -0.475966, -0.67704, 0.996452, -0.0685748, -0.851985, 0.416498, 0.791047, -0.211323, -0.302819, 0.640735, -0.317908, -0.116586, -0.896382, -0.817317, -0.948837, -0.597427},
                {0.975863, -0.971371, -0.124115, 0.4339, -0.254671, 0.298068, -0.349803, -0.73185, 0.488106, -0.0495073, 0.253969, 0.168116, 0.148772, 0.889593, -0.512213, -0.165437, 0.666773, -0.976304, -0.170024, 0.905794, 0.473908, -0.855725, -0.0413591, -0.508661, 0.443453, 0.842925, -0.144503, 0.936699, -0.443935, -0.182996},
                {0.803564, 0.960386, -0.0323329, 0.638181, -0.895684, -0.360502, 0.0646258, -0.202449, -0.717228, 0.970489, 0.404608, -0.0861868, -0.879417, -0.866462, -0.938336, -0.799685, 0.213464, -0.932344, -0.668236, 0.751366, -0.22712, -0.407783, 0.657463, 0.0970092, -0.579109, -0.868866, -0.504041, 0.926483, 0.169978, -0.00842563},
                {-0.530324, 0.282745, 0.0255867, 0.287686, 0.410417, -0.766576, -0.536566, -0.628288, 0.69665, 0.820713, -0.506162, -0.404114, 0.640099, -0.956123, -0.576586, 0.435502, -0.470676, -0.367062, -0.831765, -0.294942, 0.518991, 0.922338, 0.337886, -0.67474, -0.725667, 0.916684, 0.39175, 0.759081, 0.496979, -0.200691},
                {0.0417966, -0.687391, 0.438773, 0.287357, 0.316636, -0.262311, -0.0755541, -0.442313, 0.621378, 0.670105, 0.060982, 0.944162, 0.643442, -0.750684, -0.639973, 0.217424, 0.592823, 0.929094, -0.239135, -0.41628, 0.570893, -0.0798988, -0.917135, -0.749545, -0.982047, 0.0626998, -0.977963, 0.660401, 0.470569, -0.0528868},
                {-0.00138645, 0.931065, -0.748519, 0.304188, -0.266153, 0.672524, -0.105179, -0.874749, -0.154355, -0.774656, -0.69654, 0.433098, 0.615897, -0.387919, -0.429779, 0.650202, 0.122306, -0.237727, 0.626817, -0.227929, 0.405916, 0.483328, 0.282047, -0.262206, 0.784123, 0.83125, -0.662272, 0.702768, 0.875814, -0.701221},
                {0.553793, 0.471795, 0.769147, 0.059668, -0.841617, -0.191179, -0.972471, -0.825361, 0.779826, -0.917201, 0.43272, 0.10301, 0.358771, 0.793448, -0.0379954, -0.870112, 0.600442, -0.990603, 0.549151, 0.512146, -0.795843, 0.490091, 0.372046, -0.549437, 0.0964285, 0.753047, -0.86284, -0.589688, 0.178612, -0.720358}};

            NOMAD::Point y(n, 0);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    y[i] += M[i][j] * x[j];
                }
            }

            NOMAD::Double f1 = y[0] * y[0];

            NOMAD::Double g = 1;
            for (int i = 1; i < n; ++i)
            {
                g += (9.0 / (n - 1)) * y[i] * y[i];
            }
            NOMAD::Double h = 1 - (f1.value() / g.value()) * (f1.value() / g.value());
            NOMAD::Double f2 = g * h;

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
        int n = 30;
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

        NOMAD::Point lb(n, 0.0);
        NOMAD::Point ub(n, 1.0);

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

        p.set_HISTORY_FILE("L2ZDT2_2_bimads_"+to_string(seed)+".txt");
        p.set_STATS_FILE("test_L2ZDT2.txt", "BBE OBJ");

        // // disable models
        // p.set_DISABLE_MODELS();
        //
        // // disable NM search
        // p.set_NM_SEARCH(false);

        // parameters validation:
        p.check();

        // custom evaluator creation:
        L2ZDT2 ev(p);

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
