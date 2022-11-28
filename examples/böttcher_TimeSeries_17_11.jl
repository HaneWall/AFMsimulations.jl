using DataFrames
using CSV
using CairoMakie
using DSP
using LinearAlgebra
using AFMsimulations
CairoMakie.activate!(type="svg")

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
n_batches = length(batches_idx)
# length of each individual batch 
l_batches = [length(batches_idx[i]) for i in eachindex(batches_idx)]
batches = [zeros(Float64, L) for L in l_batches]
t_batches = [zeros(Float64, L) for L in l_batches] 
for idx in eachindex(batches_idx) 
    batches[idx] .= data_x[batches_idx[idx]]
    t_batches[idx] .= data_t[batches_idx[idx]] 
end 


begin_t = [t_batches[idx][1] for idx in eachindex(t_batches)]

"""
Apply savitzky_golay_filter to the batches to archieve denoised Signal and derivatives --> SINDy compatible
"""

batch_to_look_for = 650
window_size = 11
poly_fit = 3
fig = Figure(resolution=(1000, 800))
ax = Axis(fig[1, 1], title="Last 600 Points of batch")
lines!(ax, t_batches[batch_to_look_for][end-600:end],batches[batch_to_look_for][end-600:end], linewidth=1.5)
lines!(ax, t_batches[batch_to_look_for][end-600:end], savitzky_golay_filter(batches[batch_to_look_for][end-600:end], window_size, poly_fit), color=:red, linewidth=1.5, linestyle=:dash)
lines!(ax, t_batches[batch_to_look_for][end-600:end], savitzky_golay_filter(batches[batch_to_look_for][end-600:end], window_size, poly_fit; deriv_order=1), color=:green)

ax_d = Axis(fig[2, 1], title="First 600 points of batch")
lines!(ax_d, t_batches[batch_to_look_for][1:600],batches[batch_to_look_for][1:600], linewidth=1.5)
lines!(ax_d, t_batches[batch_to_look_for][1:600], savitzky_golay_filter(batches[batch_to_look_for][1:600], window_size, poly_fit), color=:red, linewidth=1.5)
lines!(ax_d, t_batches[batch_to_look_for][1:600], savitzky_golay_filter(batches[batch_to_look_for][1:600], window_size, poly_fit; deriv_order=1), color=:green)
fig





#function stft_batch(sig::Array{Float64}, fs::Float6)


# fig = Figure(resolution = (800, 800), fontsize=24)
# spec = spectrogram(batches[1000], 1_100, 200; fs=Δf, window=blackman)
# ax = Axis(fig[1, 1], ylabel=L"f[Hz]", xlabel=L"t[s]")
# hmap = heatmap!(ax, spec.time, spec.freq, transpose(pow2db.(spec.power)), colormap = :turbo, colorrange=(-50, 25))
# Colorbar(fig[2, 1], hmap; label = "Power in dB", ticksize = 15, tickalign = 1, vertical = false)
# ylims!(ax, [0.e4, 6.1e4])
# fig



# N = 8_000_000
# # Do a STFT
# fig = Figure(resolution = (800, 800), fontsize=24)
# spec = spectrogram(data[1:N, :mV], 70_000, 4_000; fs=Δf, window=blackman)
# ax = Axis(fig[1, 1], ylabel=L"f[Hz]", xlabel=L"t[s]")
# hmap = heatmap!(ax, spec.time, spec.freq, transpose(pow2db.(spec.power)), colormap = :turbo, colorrange=(-25, 25))
# Colorbar(fig[2, 1], hmap; label = "Power in dB", ticksize = 15, tickalign = 1, vertical = false)
# ylims!(ax, [4.8e4, 6.1e4])
# fig



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
