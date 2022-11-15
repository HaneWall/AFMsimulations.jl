using DataFrames
using CSV

data = DataFrame(CSV.File("./TimeSeries_3_14_11.csv", header=true))

"""
BIG PROBLEM: How to determine the end of one frequency. 
Todos:
    - Implement Gabor Transform to determine endpoints of one constant frequency 
    - create batches for each individual frequency --> SINDy determines coefficients for one model
"""