import Pkg
import Dates
Pkg.activate("../../../DMultiMadsPB")

import DMultiMadsPB

function solar8(x, seed::Int)

    # Check x dimension
    if length(x) != 11
        error("Error : n != 11")
    end

    # Create input file
    # If used in a distributed environment, the name should not be the same
    tag = string(Dates.DateTime(Dates.now()))
    tag = replace(tag,"-" => "")
    tag = replace(tag,":" => "")
    input_coordinates_filename = "solar8_tmp_x_PB_" * string(seed) * "_" * tag * ".txt"

    inputs_solar8 = Array{Float64}(undef, 13)
    # Real entries
    inputs_solar8[1:5] .= x[1:5]
    inputs_solar8[7:9] .= x[6:8]
    inputs_solar8[11:13] .= x[9:11]
    # Integer entries: as this version does not support integer entries, they are fixed
    # to the x0 entries proposed.
    inputs_solar8[6] = 2650 # Maximum number of heliostats
    inputs_solar8[10] = 36 # Receiver number of tubes

    open(input_coordinates_filename, "w") do f
        for elt in inputs_solar8
            write(f, "$elt\n")
        end
    end

    output_coordinates_filename = "solar8_tmp_outputs_PB_" * string(seed) * "_" * tag * ".txt"
    # Launch simulation
    cmd = run(pipeline(`./solar_bb.exe 8 $input_coordinates_filename`,
                      stdout=output_coordinates_filename))
    flag = cmd.exitcode

    # Execution fails
    if (flag == 1)
        return Inf * ones(11)
    end

    # Get the result
    outputs_bb = open(output_coordinates_filename) do f
        read(f, String)
    end

    # Check the outputs
    v = begin
        if occursin("ERROR", outputs_bb)
            Inf * ones(11)
        else
            map(elt-> parse(Float64, elt), split(outputs_bb))
        end
    end

    # Delete the temporary files
    rm(input_coordinates_filename)
    rm(output_coordinates_filename)

    # For NaN values
    if any(isnan.(v))
        return Inf * ones(11)
    end

    return v
end

for seed in [0]
    for strategy in ["PB", "EB", "Penalty"]
        bbSolar8 = DMultiMadsPB.BBProblem(x -> solar8(x, seed),
                                          11, 11,
                                          [repeat([DMultiMadsPB.OBJ], 2); repeat([DMultiMadsPB.CSTR], 9)],
                                          lvar=[1.0, 1.0, 20.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.01, 0.005, 0.0060],
                                          uvar=[40.0, 40.0, 250.0, 30.0, 30.0, 89.0, 20.0, 20.0, 5.00, 0.100, 0.1000],
                                          name="Solar8")

        model = DMultiMadsPB.MadsModel(bbSolar8)
        model.options.neval_bb_max = 5000
        model.params.seed = seed
        model.options.display = true

        if strategy == "EB"
            model.params.h_init = 0
        elseif strategy == "Penalty"
            model.options.use_penalty_approach = true
        end

        DMultiMadsPB.solve!(model, [[11.0, 11.0, 200.0, 10.0, 10.0, 89.0, 0.5, 8.0, 0.30, 0.020, 0.0216]])

        DMultiMadsPB.save_cache(model.cache, "SOLAR8_dmultimads_" * strategy * "_" * string(seed) * ".txt")
    end
end
