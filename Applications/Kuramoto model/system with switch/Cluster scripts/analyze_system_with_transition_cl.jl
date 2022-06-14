using JLD2
using FileIO
using InterSpikeSpectra
using DelayEmbeddings
using RecurrenceAnalysis
using DelimitedFiles
using Distances


f = JLD2.jldopen("res_transition.jld2", "r")
phases = Dataset(read(f,"phases")')
frequencies = Dataset(read(f,"frequencies")')

# discard transients and downsampling
frequencies = frequencies[501:end,:]
phases = phases[501:end,:]

t = 1:length(phases)

## windowed Analysis
ε = 0.05
rho1 = 0.8
rho2 = 0.85
rho3 = 0.9
rho4 = 0.95
M = 500
window = 600

N = length(phases)

spectrum1 = zeros(length(1:N-window),Int(M/2))
spectrum2 = zeros(length(1:N-window),Int(M/2))
spectrum3 = zeros(length(1:N-window),Int(M/2))
spectrum4 = zeros(length(1:N-window),Int(M/2))
rhos1 = zeros(length(1:N-window))
rhos2 = zeros(length(1:N-window))
rhos3 = zeros(length(1:N-window))
rhos4 = zeros(length(1:N-window))

method = lasso()

for i = 1:N-window
    #println("i is: $i")
    Y = phases[i:i+window,:]

    RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε; metric=Distances.PeriodicEuclidean(2*pi*ones(size(Y,2))), fixedrate = true)
    τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
    τ_rr = τ_rr[1:M] ./ maximum(τ_rr[1:M])
    spectrum1[i,:], rhos1[i] = inter_spike_spectrum(τ_rr; ρ_thres = rho1, method, regression_type=normal())
    spectrum2[i,:], rhos2[i] = inter_spike_spectrum(τ_rr; ρ_thres = rho2, method, regression_type=normal())
    spectrum3[i,:], rhos3[i] = inter_spike_spectrum(τ_rr; ρ_thres = rho3, method, regression_type=normal())
    spectrum4[i,:], rhos4[i] = inter_spike_spectrum(τ_rr; ρ_thres = rho4, method, regression_type=normal())

end

writedlm("spectrum_0_8.csv",spectrum1)
writedlm("spectrum_0_85.csv",spectrum2)
writedlm("spectrum_0_9.csv",spectrum3)
writedlm("spectrum_0_95.csv",spectrum4)
writedlm("rhos_0_8.csv", rhos1)
writedlm("rhos_0_85.csv", rhos2)
writedlm("rhos_0_9.csv", rhos3)
writedlm("rhos_0_95.csv", rhos4)