module AFMsimulations

    include("experiment_structs.jl")
    export Cantilever, Sample, AFM_LJ_experiment, AFM_DMT_experiment, AFM_eDMT_experiment, AFM_vDMT_experiment

    include("force_models.jl")
    export f_LJ!, f_DMT!, f_eDMT!, f_vDMT!

    include("utilities.jl")
    export eff_young_module, force_distance

    include("sweeps.jl")
    export freq_sweep

end
