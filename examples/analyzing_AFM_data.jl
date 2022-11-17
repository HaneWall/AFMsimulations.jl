using DataFrames
using CSV
using CairoMakie
using DSP 

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

CairoMakie.activate!()
#N = 70_000_000
spec = spectrogram(data[1:end, :mV], 100_000, 4_000; fs=Δf, window=blackman)

fig = Figure(resolution = (800, 1200), fontsize=24)
ax = Axis(fig[1, 1], xlabel=L"f[Hz]", ylabel=L"t[s]")
hmap = heatmap!(ax, spec.freq, spec.time, pow2db.(spec.power), colormap = :turbo, colorrange=(-25, 20))
Colorbar(fig[1, 2], hmap; label = "Power in dB", width = 15, ticksize = 15, tickalign = 1)
xlims!(ax, [4.8e4, 6.1e4])
fig