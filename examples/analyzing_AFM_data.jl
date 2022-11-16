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


N = 1_000_000
number_of_peices = 10_000
samples_per_fft = N÷number_of_peices
spec = spectrogram(data[1:N, :mV], N, samples_per_fft; fs=Δf)

fig = Figure(resolution = (500, 500))
ax = Axis(fig[1, 1])
heatmap!(ax, collect(0:1:N).*Δt, spec.freq, spec.power)
fig