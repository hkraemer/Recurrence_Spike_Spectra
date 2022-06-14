using JLD2
using FileIO
using tau_recurrence
using DelayEmbeddings
using RecurrenceAnalysis
using DelimitedFiles
using PyPlot
using Distances
pygui(true)

f = JLD2.jldopen("./Applications/Kuramoto model/data/res_transition.jld2", "r")
phases = Dataset(read(f,"phases")')
frequencies = Dataset(read(f,"frequencies")')

# discard transients and downsampling
frequencies = frequencies[501:end,:]
phases = phases[501:end,:]

t = 1:length(phases)

## Computations have been carried out on the cluster (windowed spectrum determination)
# in the script `analyze_system_with_transition_cl.jl`

# load the results:

# λ1 = 1e-4:
spectrum1 = readdlm("./Applications/Kuramoto model/system with switch/results.spectrum1_trial.csv")
# λ2 = 5e-5:
