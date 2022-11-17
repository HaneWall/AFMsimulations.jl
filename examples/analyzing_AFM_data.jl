using DataFrames
using CSV
using CairoMakie
using DSP 

data = DataFrame(CSV.File("./TimeSeries_3_14_11.csv", header=true))

"""
BIG PROBLEM: How to determine the end of one frequency. 
Todos:
    - Implement Gabor Transform to determine endpoints of one constant frequency 
    - create batches for each individual frequency --> SINDy determines coefficients for one model
"""

# around 80_000_000 rows, lets first concentrate on the first 20_000_000

# let us first get ther sample time (s) and the sample frequncy (Hz) 

Δt = (data[2, :ms] - data[1, :ms])*1e-3 
Δf = 1/Δt

#N = 70_000_000
spec = spectrogram(data[1:end, :mV], 50_000, 1_000; fs=Δf)



fig = Figure(resolution = (500, 500))
ax = Axis(fig[1, 1], xlabel="f[Hz]")
heatmap!(ax, spec.freq, spec.time, pow2db.(spec.power))
xlims!(ax, [5e4, 6e4])
fig