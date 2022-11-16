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



# fs = 44_100  # Hz
# ts = range(0, stop=5, step=1/fs)  # seconds
# signal = @. sin(2π*1_000*ts^2)   # sweep over f = (1000*ts) Hz
# #  plot(ts, signal)   # > 200K points, better to use InspectDR
# n = length(signal)
# nw = n÷50
# spec = spectrogram(signal, nw, nw÷2; fs=fs)


N = 1_000_000
number_of_peices = 10
samples_per_fft = N÷number_of_peices
spec = spectrogram(data[1:N, :mV], 300; fs=Δf)



fig = Figure(resolution = (500, 500))
ax = Axis(fig[1, 1])
heatmap!(ax, pow2db.(spec.power))
fig