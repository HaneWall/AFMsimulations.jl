using DataFrames
using CSV
using CairoMakie
using DSP 
using LinearAlgebra
CairoMakie.activate!()

# data: 
data = DataFrame(CSV.File("./TimeSeries_1_17_11_22_JumpFlag.csv", header=true))

Δt = (data[2, :s] - data[1, :s])
Δf = 1/Δt




N = length(data[1:end, :s])

#N = _000_000
#Do a STFT
fig = Figure(resolution = (800, 800), fontsize=24)
spec = spectrogram(data[1:N, :mV], 70_000, 4_000; fs=Δf, window=blackman)
ax = Axis(fig[1, 1], ylabel=L"f[Hz]", xlabel=L"t[s]")
hmap = heatmap!(ax, spec.time, spec.freq, transpose(pow2db.(spec.power)), colormap = :turbo, colorrange=(-25, 25))
Colorbar(fig[2, 1], hmap; label = "Power in dB", ticksize = 15, tickalign = 1, vertical = false)
ylims!(ax, [4.8e4, 6.1e4])
fig



"""
Batches are pieces of the overall data at constant frequency. All batches from
one signal appended at each other result in the original signal again. 
"""
struct batches
    idx_min::Int
    idx_max::Int
    data_t::Array{Float64}
    data_x::Array{Float64}
    function batches(signal::DataFrame, idx_min, idx_max)
        data_t = signal[idx_min:idx_max, 1]
        data_x = signal[idx_min:idx_max, 2]
        new(idx_min, idx_max, data_t, data_x)
    end
end


"""
idea: create Bit-array out of rect array
bit_array = rect_array >. 0
save the 
"""
function create_batches_array(rect::Array{Float64})
    return 0
end

"""
We only have the position of the cantilever. In order to get the derivative, we
should either use the TV-diff or a central differences stencil scheme. This can
then be used on the smaller batches we created beforehand. 
"""
function tv_diff(noisy_arr::Array{Float64})
    return 0
end 