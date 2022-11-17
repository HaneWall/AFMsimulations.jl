var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = AFMsimulations","category":"page"},{"location":"#AFMsimulations","page":"Home","title":"AFMsimulations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for AFMsimulations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [AFMsimulations]","category":"page"},{"location":"#AFMsimulations.AFM_DMT_experiment","page":"Home","title":"AFMsimulations.AFM_DMT_experiment","text":"DMT experiment, respects elastic behavior within the hertzian contact model \n\n\n\n\n\n","category":"type"},{"location":"#AFMsimulations.AFM_LJ_experiment","page":"Home","title":"AFMsimulations.AFM_LJ_experiment","text":"Lennard Jones experiment\n\n\n\n\n\n","category":"type"},{"location":"#AFMsimulations.AFM_eDMT_experiment","page":"Home","title":"AFMsimulations.AFM_eDMT_experiment","text":"DMT experiment with exponential damping term –> eDMT\n\n\n\n\n\n","category":"type"},{"location":"#AFMsimulations.AFM_vDMT_experiment","page":"Home","title":"AFMsimulations.AFM_vDMT_experiment","text":"DMT experiment with viscious damping term –> vDMT\n\n\n\n\n\n","category":"type"},{"location":"#AFMsimulations.AFM_vLJ_experiment","page":"Home","title":"AFMsimulations.AFM_vLJ_experiment","text":"Lennard Jones experiment with damping\n\n\n\n\n\n","category":"type"},{"location":"#AFMsimulations.Cantilever","page":"Home","title":"AFMsimulations.Cantilever","text":"R .. radius of tip, Q .. quality factor of tip, k .. spring constant,  E .. Young's module, ν .. poisson coefficient, f_0 .. eigenfrequency. \n\n\n\n\n\n","category":"type"},{"location":"#AFMsimulations.Sample","page":"Home","title":"AFMsimulations.Sample","text":"H .. Hamaker constant, E .. Young's module, ν .. poisson coefficient.\n\n\n\n\n\n","category":"type"},{"location":"#AFMsimulations.ampl_sweep-Tuple{Array{Float64}, Float64, AFM_DMT_experiment, Float64}","page":"Home","title":"AFMsimulations.ampl_sweep","text":"Amplitude sweeps for DMT-models.\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.ampl_sweep-Tuple{Array{Float64}, Float64, AFM_vLJ_experiment, Float64}","page":"Home","title":"AFMsimulations.ampl_sweep","text":"Amplitude-Sweeps. Keep the frequency constant, instead vary the amplitude of the forcing. This allows us to plot input amplitude against response amplitude. This curve represents a slice in the bifurcation diagramm of the frequncy against response amplitude. \n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.eff_young_module-Tuple{Sample, Cantilever}","page":"Home","title":"AFMsimulations.eff_young_module","text":"returns effective Young's module of two interacting materials\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.f_DMT!-NTuple{4, Any}","page":"Home","title":"AFMsimulations.f_DMT!","text":"DMT model. Respects elastic properties via Hertzian contact model. \n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.f_DMT_auto!-NTuple{4, Any}","page":"Home","title":"AFMsimulations.f_DMT_auto!","text":"Autonomous DMT model. Respects elastic properties via Hertzian contact model. \n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.f_LJ!-NTuple{4, Any}","page":"Home","title":"AFMsimulations.f_LJ!","text":"Lennard-Jones model, smooth interaction potential.\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.f_eDMT!-NTuple{4, Any}","page":"Home","title":"AFMsimulations.f_eDMT!","text":"DMT model with additional exponential decaying damping term, inspired by Haviland.\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.f_vDMT!-NTuple{4, Any}","page":"Home","title":"AFMsimulations.f_vDMT!","text":"DMR model with additional viscious damping term, only acting for x > a_0\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.f_vLJ!-NTuple{4, Any}","page":"Home","title":"AFMsimulations.f_vLJ!","text":"Lennard-Jones model, smooth interaction, soft sphere model to avoid singularities at the surface of the sample. Introduces damping via (x - d)^{-3} term. See I. Barke implementation. \n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.force_distance-Tuple{Float64, AFM_DMT_experiment}","page":"Home","title":"AFMsimulations.force_distance","text":"returns force-distance relation for DMT experiments\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.force_distance-Tuple{Float64, AFM_vLJ_experiment}","page":"Home","title":"AFMsimulations.force_distance","text":"returns force distance relation for LJ experiments\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.freq_sweep-Tuple{Float64, Array{Float64}, AFM_DMT_experiment, Float64}","page":"Home","title":"AFMsimulations.freq_sweep","text":"Forward or backward frequency sweep with constant excitation amplitude for DMT-model.\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.freq_sweep-Tuple{Float64, Array{Float64}, AFM_LJ_experiment, Float64}","page":"Home","title":"AFMsimulations.freq_sweep","text":"Forward or backward frequency sweep with constant excitation amplitude for LJ-model.\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.freq_sweep-Tuple{Float64, Array{Float64}, AFM_eDMT_experiment, Float64}","page":"Home","title":"AFMsimulations.freq_sweep","text":"same for eDMT models\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.freq_sweep-Tuple{Float64, Array{Float64}, AFM_vDMT_experiment, Float64}","page":"Home","title":"AFMsimulations.freq_sweep","text":"same for vDMT models\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.freq_sweep-Tuple{Float64, Array{Float64}, AFM_vLJ_experiment, Float64}","page":"Home","title":"AFMsimulations.freq_sweep","text":"Forward or backward frequency sweep with constant excitation amplitude for vLJ-model.\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.freq_sweep_auto-Tuple{Float64, Array{Float64}, AFM_DMT_experiment, Float64}","page":"Home","title":"AFMsimulations.freq_sweep_auto","text":"Forward / backward frequency sweep with constant excitation amplitude for DMT-model.\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.freq_sweep_Φ-Tuple{Float64, Array{Float64}, AFM_vLJ_experiment, Float64}","page":"Home","title":"AFMsimulations.freq_sweep_Φ","text":"Frequency sweep of LJ-potential with fixpoints x^{⋆} of Poincare map Φ: x(t) –> x(t + T) with T being the period of the forcing term. Should give even more accurate results and might be more faster (due to less iterations).\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.lock_in_amplifier-Tuple{Array{Float64}, Float64}","page":"Home","title":"AFMsimulations.lock_in_amplifier","text":"Implements digital lock in amplifier. Provide a reference signal (a local osciallator with ω=ωref) and gain information of the phase and amplitude of the input signal.  Input:: sig:Array of interest, ωref:reference frequency \n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.phase_space-Tuple{AFM_vLJ_experiment}","page":"Home","title":"AFMsimulations.phase_space","text":"Creates phase space representation for the damped Lennard-Jones model.\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.steady_state_check-Tuple{Array{Float64}, Float64}","page":"Home","title":"AFMsimulations.steady_state_check","text":"Steady state checker, stop simulation if std of amplitude for n periods stays below some tolerance: tol  \n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.Φ_poincare-Tuple{SciMLBase.ODEProblem, Array{Float64}, Array{Float64}, Float64}","page":"Home","title":"AFMsimulations.Φ_poincare","text":"Poincare map Φ: x(t) –> x(t + T), with T being the period of the forcing term. For us T = 2π * Ω. Input is a ODE-Problem, initial Vector x, problem specific parameters and Δt.\n\n\n\n\n\n","category":"method"},{"location":"#AFMsimulations.Φ_poincare_fixpoint-Tuple{SciMLBase.ODEProblem, Array{Float64}, Array{Float64}, Float64, Float64}","page":"Home","title":"AFMsimulations.Φ_poincare_fixpoint","text":"Poincare map Φ: x(t) –> x(t + T), with T being the period of the forcing term. For us T = 2π * Ω. Input is a ODE-Problem, initial Vector x, problem specific parameters and Δt. We found a periodic orbit iff |Φ(x) - x|_2 <= ϵ with ϵ being the absolute tolerance. This function can tzhan be used to determine a callback function in an ODEProblem. This way we Should be able to get way more efficient algorithms. \n\n\n\n\n\n","category":"method"}]
}