using DelimitedFiles
using Distances
using OptimalTransport
using Distributions

## We compute the pairwise earth mover's distance of the computed spectra of 
## the 50-dimensional Kuramoto-model with switch (computations have been carried 
## out on the cluster -windowed spectrum determination- in the script `analyze_system_with_transition_cl.jl`

# load the results:

spectrum1 = readdlm("./Applications/Kuramoto model/system with switch/results/spectrum_0_8.csv")
spectrum2 = readdlm("./Applications/Kuramoto model/system with switch/results/spectrum_0_85.csv")
spectrum3 = readdlm("./Applications/Kuramoto model/system with switch/results/spectrum_0_9.csv")
spectrum4 = readdlm("./Applications/Kuramoto model/system with switch/results/spectrum_0_95.csv")


N, M = size(spectrum1)

emd1 = zeros(N,N)
emd2 = zeros(N,N)
emd3 = zeros(N,N)
emd4 = zeros(N,N)

Threads.@threads for i = 1:N
    for j = 1:N
        emd1[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum1[i,:]), DiscreteNonParametric(1:M,spectrum1[j,:]); p=2)
        emd2[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum2[i,:]), DiscreteNonParametric(1:M,spectrum2[j,:]); p=2)
        emd3[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum3[i,:]), DiscreteNonParametric(1:M,spectrum3[j,:]); p=2)
        emd4[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum4[i,:]), DiscreteNonParametric(1:M,spectrum4[j,:]); p=2)
    end
end

writedlm("./Applications/Kuramoto model/system with switch/results/emd1.csv",emd1)
writedlm("./Applications/Kuramoto model/system with switch/results/emd2.csv",emd2)
writedlm("./Applications/Kuramoto model/system with switch/results/emd3.csv",emd3)
writedlm("./Applications/Kuramoto model/system with switch/results/emd4.csv",emd4)

