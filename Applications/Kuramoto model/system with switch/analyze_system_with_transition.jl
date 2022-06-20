using JLD2
using FileIO
using InterSpikeSpectra
using Statistics
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
N, M = size(phases)

t = 1:N # time vector

# define a artificial measurement function
function measurement_function(x, a, b, c, σ=0)
    return a*x + b*x^2 + c*x^3 + σ*randn()
end

# transform phases
a = 0.001
b = 0.0043
c = 0.002
n = 0.1
phases_transform = measurement_function.(Matrix(phases), a, b, c, n)

## Compute order parameter:
rs = zeros(N)
rs_transform = zeros(N)
for k in 1:N
    rs[k] = abs.(mean((exp.(im .* (phases[k,:])))))
    rs_transform[k] = abs.(mean((exp.(im .* (phases_transform[k,:])))))
end

figure()
subplot(211)
plot(t,rs)
title("Order parameter r for 50-dim. Kuramoto with switch")
ylabel("order parameter r")
grid()
subplot(212)
plot(t,rs_transform)
title("Order parameter r for 50-dim. transformed Kuramoto with switch")
xlabel("time [a.u.]")
ylabel("order parameter r")
grid()


figure()
subplot(211)
plot(t,Matrix(phases))
title("Phases for 50-dim. Kuramoto with switch")
grid()
subplot(212)
plot(t,phases_transform)
title("Transformed phases for 50-dim. Kuramoto with switch")
grid()

