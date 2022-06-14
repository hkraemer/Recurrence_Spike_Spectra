using JLD2
using FileIO
using InterSpikeSpectra
using DelayEmbeddings
using RecurrenceAnalysis
using DelimitedFiles
using PyPlot
using Distances
pygui(true)

f = JLD2.jldopen("./Applications/Kuramoto model/data/res_sigma001.jld2", "r")
phases = Dataset(read(f,"phases"))
frequencies = Dataset(read(f,"frequencies"))

# discard transients and downsampling
frequencies = frequencies[501:end,:]
phases = phases[501:end,:]

t = 1:length(phases)

# plot all phases
figure()
for i = 1:50
    plot(phases[:,i])
end
grid()

# plot all frequencies
figure()
for i = 1:50
    plot(frequencies[:,i])
end
grid()

# plot phase difference
figure()
plot(diff(Matrix(phases[:,1:2]); dims=2))
grid()

# plot phase difference
figure()
plot(diff(Matrix(frequencies[:,1:2]); dims=2))
grid()

phases2 = mod2pi.(Matrix(phases))

Y = Dataset(phases)
ε = 0.05
RP = RecurrenceAnalysis.RecurrenceMatrix(Y[1:3000,:], ε; metric=Distances.PeriodicEuclidean(2*pi*ones(50)), fixedrate = true)
τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
τ_rr ./= maximum(τ_rr)

plot(τ_rr)

ρ_thres = 0.8
method = lasso()
@time spectrum, rho = inter_spike_spectrum(τ_rr[1:1000]; ρ_thres, method, regression_type=normal())


Rg = grayscale(RP)
figure()
imshow(Rg, cmap = "binary_r", extent = (1, size(RP)[1], 1, size(RP)[2]))

@time begin
    figure()
    subplot(211)
    title("τ-RR spectrum of 50-dim. synchronized Kuramoto")
    plot(τ_rr)
    ylabel("τ-RR")
    grid()
    subplot(212)
    plot(spectrum)
    title("Inter spike spectrum of τ-RR of 50-dim. synchronized Kuramoto")
    xlabel("period [a.u.]")
    ylabel("rel. power per period")
    grid()
end
