import Pkg
import Dates
Pkg.activate("./../../../DMultiMadsPB")

import DMultiMadsPB

function solar9(x, seed::Int)

    # Check x dimension
    if length(x) != 22
        error("Error : n != 22")
    end

    # Create input file
    # If used in a distributed environment, the name should not be the same
    tag = string(Dates.DateTime(Dates.now()))
    tag = replace(tag,"-" => "")
    tag = replace(tag,":" => "")
    input_coordinates_filename = "solar9_tmp_x_" * string(seed) * "_" * tag * ".txt"

    inputs_solar9 = Array{Float64}(undef, 29)
    # Real entries
    inputs_solar9[1:5] .= x[1:5]
    inputs_solar9[7:15] .= x[6:14]
    inputs_solar9[17:24] .= x[15:22]
    # Integer entries: as this version does not support integer entries, they are fixed
    # to the x0 entries proposed.
    inputs_solar9[6] = 1000 # Maximum number of heliostats
    inputs_solar9[16] = 500 # Receiver number of tubes
    inputs_solar9[25] = 3 # Exchanger number of baffles
    inputs_solar9[26] = 12000 # Exchanger number of tubes
    inputs_solar9[27] = 1 # Exchanger number of shells
    inputs_solar9[28] = 2 # Exchanger number of passes per shell
    inputs_solar9[29] = 2 # Type of turbine

    open(input_coordinates_filename, "w") do f
        for elt in inputs_solar9
            write(f, "$elt\n")
        end
    end

    output_coordinates_filename = "solar9_tmp_outputs_" * string(seed) * "_" * tag * ".txt"
    # Launch simulation
    cmd = run(pipeline(`./solar_bb.exe 9 $input_coordinates_filename`,
                      stdout=output_coordinates_filename))
    flag = cmd.exitcode

    # Execution fails
    if (flag == 1)
        return Inf * ones(19)
    end

    # Get the result
    outputs_bb = open(output_coordinates_filename) do f
        read(f, String)
    end

    # Check the outputs
    v = begin
        if occursin("ERROR", outputs_bb)
            Inf * ones(19)
        else
            map(elt-> parse(Float64, elt), split(outputs_bb))
        end
    end

    # Delete the temporary files
    rm(input_coordinates_filename)
    rm(output_coordinates_filename)

    # For NaN values
    if any(isnan.(v))
        open("curious_failures_solar9.txt", "a") do f
            for elt in inputs_solar9
                write(f, "$elt ")
            end
            write(f, "\n")
        end
        return Inf * ones(19)
    end

    return v
end

for seed in [0]
    for strategy in ["PB", "EB", "Penalty"]
        bbSolar9 = DMultiMadsPB.BBProblem(x -> solar9(x, seed),
                                          22, 19,
                                          [repeat([DMultiMadsPB.OBJ], 2); repeat([DMultiMadsPB.CSTR], 17)],
                                          lvar=[1.0, 1.0, 20.0, 1.0, 1.0, 1.0, 0.0, 1.0, 793.0, 1.0, 1.0, 0.01, 0.01, 495.0, 0.01, 0.0050, 0.006, 0.007, 0.5, 0.0050, 0.006, 0.15],
                                          uvar=[40.0, 40.0, 250.0, 30.0, 30.0, 89.0, 20.0, 20.0, 995.0, 50.0, 30.0, 5.00, 5.00, 650.0, 5.00, 0.1000, 0.100, 0.200, 10.0, 0.1000, 0.100, 0.40],
                                          name="Solar9")

        model = DMultiMadsPB.MadsModel(bbSolar9)
        model.options.neval_bb_max = 5000
        model.params.seed = seed
        model.options.display = true
        if strategy == "EB"
            model.params.h_init = 0
        elseif strategy == "Penalty"
            model.options.use_penalty_approach = true
        end

        DMultiMadsPB.solve!(model,
                            [[9.0, 9.0, 150.0, 6.0, 8.0, 45.0, 0.5, 5.0, 900.0, 9.0, 9.0, 0.30, 0.20, 560.0, 0.30, 0.0165, 0.018, 0.017, 10.0, 0.0155, 0.016, 0.20]])

        DMultiMadsPB.save_cache(model.cache, "SOLAR9_dmultimads_" * strategy * "_" * string(seed) * ".txt")
    end
end
