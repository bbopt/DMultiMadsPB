import Pkg
import Dates
Pkg.activate("../../../DMultiMadsPB")

import DMultiMadsPB

function styrene(x, seed::Int)

    # Check x dimension
    if length(x) != 8
        error("Error : n != 8")
    end

    # Create input file
    # If used in a distributed environment, the name should not be the same
    tag = string(Dates.DateTime(Dates.now()))
    tag = replace(tag,"-" => "")
    tag = replace(tag,":" => "")
    input_coordinates_filename = "styrene_tmp_x_PB_v28_" * string(seed) * "_" * tag * ".txt"

    open(input_coordinates_filename, "w") do f
        for elt in x
            write(f, "$elt\n")
        end
    end

    output_coordinates_filename = "styrene_tmp_outputs_v28_" * string(seed) * "_" * tag * ".txt"
    # Launch simulation
    cmd = run(pipeline(`./truth.exe $input_coordinates_filename`,
                      stdout=output_coordinates_filename))
    flag = cmd.exitcode

    # Execution fails
    if (flag == 1)
        return Inf * ones(12)
    end

    # Get the result
    outputs_bb = open(output_coordinates_filename) do f
        read(f, String)
    end

    # Check the outputs
    v = begin
        if occursin("ERROR", outputs_bb)
            Inf * ones(12)
        else
            map(elt-> parse(Float64, elt), split(outputs_bb))
        end
    end

    # Delete the temporary files
    rm(input_coordinates_filename)
    rm(output_coordinates_filename)

    # Reorganize a bit the constraints to facilitate the analysis
    # v[1:4] boolean constraints: set to 1 if simulation fails, 0 otherwise.
    # v[5] : minimal purity of produced styrene (f2)
    # v[6] : minimal purity of produced benzene
    # v[7] : overall ethylbenzene conversion into styrene (f3)
    # v[8:11] : Four constraints relating to payout time, cashflow, investment and annual costs
    # v[12] : net present value of the project (f1)

    perm = [12,5,7,1,2,3,4,6,8,9,10,11]
    permute!(v, perm)

    return v
end

for seed in [1234]
    for variant in ["PB", "EB", "Penalty"]
        # The blackbox makes a scaling of the inputs before running the simulation.
        # It is then normal than lower and upper bounds diverge from the original description of
        # the problem in Mads-VNS.
        bbStyrene = DMultiMadsPB.BBProblem(x -> styrene(x, seed),
                                        8, 12,
                                        [repeat([DMultiMadsPB.OBJ], 3); repeat([DMultiMadsPB.CSTR], 9)],
                                        lvar=zeros(8),
                                        uvar=100*ones(8),
                                        name="Styrene")

        model = DMultiMadsPB.MadsModel(bbStyrene)
        model.options.neval_bb_max = 20000
        model.params.seed = seed
        model.options.display = true

        if variant == "EB"
            model.params.h_init = 0
        elseif variant == "Penalty"
            model.options.use_penalty_approach = true
        end

        DMultiMadsPB.solve!(model, [[54.0; 66.0; 86.0; 8.0; 29.0; 51.0; 32.0; 15.0]])

        DMultiMadsPB.save_cache(model.cache, "STYRENE_dmultimads_" * variant * "_" * string(seed) * ".txt")

    end
end
