using DataFrames
using CSV
using CairoMakie
using DSP
using LinearAlgebra
CairoMakie.activate!()

function create_batches(rect::Array{Float64})
    bit_array = rect .> 0
    patch_ranges = []
    idx_min = 0
    idx_max = 0
    for idx ∈ eachindex(bit_array[1:end-1])
        idx_min = idx_max + 1 
        if bit_array[idx] == bit_array[idx + 1]
            continue
        else
            idx_max = idx
            push!(patch_ranges, idx_min:idx_max)
        end
    end
    return patch_ranges
end


# data: 
data = DataFrame(CSV.File("./TimeSeries_1_17_11_22_JumpFlag.csv", header=true))

Δt = (data[2, :s] - data[1, :s])
Δf = 1/Δt




N = length(data[1:end, 1])
rects = data[1:end, 3]
data_x = data[1:end, 2]
data_t = data[1:end, 1]

"""
Batches are pieces of the overall data at constant frequency. All batches from
one signal appended at each other result in the original signal again. 
"""
batches_idx = create_batches(rects)
# number of batches 
n_batches = length(batches)
# length of each individual batch 
l_batches = [length(batches[i]) for i in eachindex(batches)]
batches = [zeros(Float64, L) for L in l_batches]
t_batches = [zeros(Float64, L) for L in l_batches] 
for idx in eachindex(batches_idx) 
    batches[idx] .= data_x[batches_idx[idx]]
    t_batches[idx] .= data_t[batches_idx[idx]] 
end 


begin_t = [t_batches[idx][1] for idx in eachindex(t_batches)]



#function stft_batch(sig::Array{Float64}, fs::Float6)


# fig = Figure(resolution = (800, 800), fontsize=24)
# spec = spectrogram(batches[1000], 1_100, 200; fs=Δf, window=blackman)
# ax = Axis(fig[1, 1], ylabel=L"f[Hz]", xlabel=L"t[s]")
# hmap = heatmap!(ax, spec.time, spec.freq, transpose(pow2db.(spec.power)), colormap = :turbo, colorrange=(-50, 25))
# Colorbar(fig[2, 1], hmap; label = "Power in dB", ticksize = 15, tickalign = 1, vertical = false)
# ylims!(ax, [0.e4, 6.1e4])
# fig



N = 8_000_000
# Do a STFT
fig = Figure(resolution = (800, 800), fontsize=24)
spec = spectrogram(data[1:N, :mV], 70_000, 4_000; fs=Δf, window=blackman)
ax = Axis(fig[1, 1], ylabel=L"f[Hz]", xlabel=L"t[s]")
hmap = heatmap!(ax, spec.time, spec.freq, transpose(pow2db.(spec.power)), colormap = :turbo, colorrange=(-25, 25))
Colorbar(fig[2, 1], hmap; label = "Power in dB", ticksize = 15, tickalign = 1, vertical = false)
ylims!(ax, [4.8e4, 6.1e4])
fig



# """
# Batches are pieces of the overall data at constant frequency. All batches from
# one signal appended at each other result in the original signal again. 
# """
# struct batches
#     idx_min::Int
#     idx_max::Int
#     data_t::Array{Float64}
#     data_x::Array{Float64}
#     function batches(signal::DataFrame, idx_min, idx_max)
#         data_t = signal[idx_min:idx_max, 1]
#         data_x = signal[idx_min:idx_max, 2]
#         new(idx_min, idx_max, data_t, data_x)
#     end
# end




"""
We only have the position of the cantilever. In order to get the derivative, we
should either use the TV-diff or a central differences stencil scheme. This can
then be used on the smaller batches we created beforehand. 
"""
function tv_diff(noisy_arr::Array{Float64})
    return 0
end 