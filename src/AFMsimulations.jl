module AFMsimulations

    include("experiment_structs.jl")
    export Cantilever, Sample, AFM_LJ_experiment, AFM_vLJ_experiment, AFM_DMT_experiment, AFM_eDMT_experiment, AFM_vDMT_experiment

    include("force_models.jl")
    export f_LJ!, f_vLJ!, f_vLJ, f_DMT!, f_DMT_auto!, f_eDMT!, f_vDMT!

    include("utilities.jl")
    export eff_young_module, force_distance, Φ_poincare, savitzky_golay_filter

    include("sweeps.jl")
    export freq_sweep, freq_sweep_stat, freq_sweep_auto, freq_sweep_Φ, phase_space, ampl_sweep

end
