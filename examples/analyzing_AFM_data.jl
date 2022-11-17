using DataFrames
using CSV
using CairoMakie
using DSP 
using LinearAlgebra
CairoMakie.activate!()

data = DataFrame(CSV.File("./TimeSeries_3_14_11.csv", header=true))

"""
BIG PROBLEM: How to determine the end of one frequency. 
Todos:
    - create batches for each individual frequency --> SINDy determines coefficients for one model
"""

# around 80_000_000 rows, lets first concentrate on the first 20_000_000

# let us first get ther sample time (s) and the sample frequncy (Hz) 

Δt = (data[2, :ms] - data[1, :ms])*1e-3 
Δf = 1/Δt

# CairoMakie.activate!()
# #N = 70_000_000
# spec = spectrogram(data[1:end, :mV], 100_000, 4_000; fs=Δf, window=blackman)

# # Big Picture
# fig = Figure(resolution = (1000, 1200), fontsize=24)
# ax_t = Axis(fig[1, 1], xlabel=L"U[mV]", ylabel=L"t[s]")
# lines!(ax_t, data[1:end, :mV], data[1:end, :ms].*1.e-3)
# ax = Axis(fig[1, 2], xlabel=L"f[Hz]", ylabel=L"t[s]")
# hmap = heatmap!(ax, spec.freq, spec.time, pow2db.(spec.power), colormap = :turbo, colorrange=(-25, 20))
# Colorbar(fig[1, 3], hmap; label = "Power in dB", width = 15, ticksize = 15, tickalign = 1)
# ylims!(ax_t, [0, data[end, :ms]*1.e-3])
# xlims!(ax, [4.8e4, 6.1e4])
# fig

# # Only forward sweep 
# N = 35_000_000
# fig = Figure(resolution = (1200, 1000), fontsize=24)
# spec = spectrogram(data[1:N, :mV], 100_000, 4_000; fs=Δf, window=blackman)
# ax_t = Axis(fig[1, 1], xlabel=L"U[mV]", ylabel=L"t[s]")
# lines!(ax_t, data[1:N, :mV], data[1:N, :ms] .*1e-3)
# ax = Axis(fig[2, 1], xlabel=L"f[Hz]", ylabel=L"t[s]")
# hmap = heatmap!(ax, spec.freq, spec.time, pow2db.(spec.power), colormap = :turbo, colorrange=(-25, 20))
# Colorbar(fig[3, 1], hmap; label = "Power in dB", width = 15, ticksize = 15, tickalign = 1, vertical = false)
# xlims!(ax, [4.8e4, 6.1e4])
# ylims!(ax_t, [0., spec.time[end]])
# fig


# Only forward sweep trasnpose
N = 35_000_000 # consider the first N rows of the data (whole is around 75_000_000)
fig = Figure(resolution = (1200, 1000), fontsize=24)
spec = spectrogram(data[1:N, :mV], 100_000, 4_000; fs=Δf, window=blackman)
ax_t = Axis(fig[1, 1], xlabel=L"t[s]", ylabel=L"U[mV]")
lines!(ax_t, data[1:N, :ms] .*1e-3, data[1:N, :mV])
ax = Axis(fig[2, 1], ylabel=L"f[Hz]", xlabel=L"t[s]")
hmap = heatmap!(ax, spec.time, spec.freq, transpose(pow2db.(spec.power)), colormap = :turbo, colorrange=(-25, 20))
Colorbar(fig[3, 1], hmap; label = "Power in dB", ticksize = 15, tickalign = 1, vertical = false)
ylims!(ax, [4.8e4, 6.1e4])
xlims!(ax_t, [0., spec.time[end]])
fig

# Check if chirp is linear. 
spec = spectrogram(data[1:33_000_000, :mV], 100_000, 4_000; fs=Δf, window=blackman)
freq_border_min = argmin(norm.(spec.freq .- 50_000))
freq_border_max = argmin(norm.(spec.freq .- 60_000))
maximum_at_0 = argmax(spec.power[freq_border_min:freq_border_max, 1]) + freq_border_min
maximum_at_N = argmax(spec.power[freq_border_min:freq_border_max, end]) + freq_border_min 
fig = Figure(resolution = (800, 800), fontsize=24)
ax = Axis(fig[1, 1], xlabel=L"t[s]", ylabel=L"f[Hz]")
hmap = heatmap!(ax, spec.time, spec.freq, transpose(pow2db.(spec.power)), colormap = :turbo, colorrange=(-25, 20))
lines!(ax, [0, spec.time[end]], [spec.freq[maximum_at_0], spec.freq[maximum_at_N]], linestyle=:dash, linewidth=3, color=:white)
Colorbar(fig[2, 1], hmap; label = "Power in dB", ticksize = 15, tickalign = 1, vertical = false)
ylims!(ax, [5.e4, 6.e4])
fig


# # Let's compare the chirp with a linear chirp
# spec = spectrogram(data[1:33_000_000, :mV], 100_000, 4_000; fs=Δf, window=blackman)
# freq_border_min = argmin(norm.(spec.freq .- 50_000))
# freq_border_max = argmin(norm.(spec.freq .- 60_000))

# maximum_at_0 = argmax(spec.power[freq_border_min:freq_border_max, 1]) + freq_border_min
# maximum_at_N = argmax(spec.power[freq_border_min:freq_border_max, end]) + freq_border_min 

# fig = Figure(resolution = (800, 400), fontsize=24)
# ax = Axis(fig[1, 1], xlabel=L"f[Hz]", ylabel=L"t[s]")
# hmap = heatmap!(ax, spec.freq, spec.time, pow2db.(spec.power), colormap = :turbo, colorrange=(-25, 20))
# lines!(ax, [spec.freq[maximum_at_0], spec.freq[maximum_at_N]], [0, spec.time[end]], linestyle=:dash, linewidth=3, color=:white)
# Colorbar(fig[1, 2], hmap; label = "Power in dB", width = 15, ticksize = 15, tickalign = 1)
# xlims!(ax, [5.e4, 6.e4])
# fig

# # zoomed in forward sweep
# spec = spectrogram(data[1:33_000_000, :mV], 100_000, 0; fs=Δf, window=blackman)
# freq_border_min = argmin(norm.(spec.freq .- 50_000))
# freq_border_max = argmin(norm.(spec.freq .- 60_000))
# fig = Figure(resolution = (800, 400), fontsize=24)
# ax = Axis(fig[1, 1], xlabel=L"f[Hz]", ylabel=L"t[s]")
# hmap = heatmap!(ax, spec.freq, spec.time, pow2db.(spec.power), colormap = :turbo, colorrange=(-25, 20))
# lines!(ax, [spec.freq[maximum_at_0], spec.freq[maximum_at_N]], [0, spec.time[end]], linestyle=:dash, linewidth=3, color=:white)
# Colorbar(fig[1, 2], hmap; label = "Power in dB", width = 15, ticksize = 15, tickalign = 1)
# xlims!(ax, [5.45e4, 5.55e4])
# ylims!(ax, [6, 10])
# fig

# # plot timeseries of the data 
# fig = Figure(resolution=(800, 400), fontsize=24)
# ax = Axis(fig[1, 1], xlabel=L"t[s]", ylabel=L"V[mV]")
# lines!(x)




