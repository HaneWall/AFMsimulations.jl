using DataFrames
using CSV

data = DataFrame(CSV.File("./TimeSeries_3_14_11.csv", header=true))

