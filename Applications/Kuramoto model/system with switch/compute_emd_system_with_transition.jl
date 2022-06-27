using DelimitedFiles
using Distances
using OptimalTransport
using Distributions

## We compute the pairwise earth mover's distance of the computed spectra of 
## the 50-dimensional Kuramoto-model with switch (computations have been carried 
## out on the cluster -windowed spectrum determination- in the script `analyze_system_with_transition_cl.jl`

# load the results:

spectrum1 = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_8.csv")
spectrum2 = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_85.csv")
spectrum3 = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_9.csv")
spectrum4 = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_95.csv")

spectrum1t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_8_t.csv")
spectrum2t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_85_t.csv")
spectrum3t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_9_t.csv")
spectrum4t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_95_t.csv")

spectrum1t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_8_t_no_noise.csv")
spectrum2t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_85_t_no_noise.csv")
spectrum3t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_9_t_no_noise.csv")
spectrum4t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_95_t_no_noise.csv")

spectrum1t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_8_t.csv")
spectrum2t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_85_t.csv")
spectrum3t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_9_t.csv")
spectrum4t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_spectrum_0_95_t.csv")

spectrum_fft = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_fft_spectrum.csv")
spectrum_fft_t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_fft_spectrum_t.csv")

spectrum_fft = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_fft_spectrum_no_noise.csv")
spectrum_fft_t = readdlm("./Applications/Kuramoto model/system with switch/results/results_kuramoto_switch_fft_spectrum_t_no_noise.csv")


N, M = size(spectrum1)

# standardize FFT-spectra to probabilities
for j = 1:N
    spectrum_fft[j,:] = spectrum_fft[j,:] ./ sum(spectrum_fft[j,:])
    spectrum_fft_t[j,:] = spectrum_fft_t[j,:] ./ sum(spectrum_fft_t[j,:])
end

# emd1 = zeros(N,N)
# emd2 = zeros(N,N)
# emd3 = zeros(N,N)
# emd4 = zeros(N,N)

emd1t = zeros(N,N)
emd2t = zeros(N,N)
emd3t = zeros(N,N)
emd4t = zeros(N,N)

emd_fft = zeros(N,N)
emd_fft_t = zeros(N,N)


Threads.@threads for i = 1:N
    for j = 1:N
        # emd1[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum1[i,:]), DiscreteNonParametric(1:M,spectrum1[j,:]); p=2)
        # emd2[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum2[i,:]), DiscreteNonParametric(1:M,spectrum2[j,:]); p=2)
        # emd3[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum3[i,:]), DiscreteNonParametric(1:M,spectrum3[j,:]); p=2)
        # emd4[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum4[i,:]), DiscreteNonParametric(1:M,spectrum4[j,:]); p=2)

        emd1t[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum1t[i,:]), DiscreteNonParametric(1:M,spectrum1t[j,:]); p=2)
        emd2t[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum2t[i,:]), DiscreteNonParametric(1:M,spectrum2t[j,:]); p=2)
        emd3t[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum3t[i,:]), DiscreteNonParametric(1:M,spectrum3t[j,:]); p=2)
        emd4t[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum4t[i,:]), DiscreteNonParametric(1:M,spectrum4t[j,:]); p=2)

        emd_fft[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum_fft[i,:]), DiscreteNonParametric(1:M,spectrum_fft[j,:]); p=2)
        emd_fft_t[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum_fft_t[i,:]), DiscreteNonParametric(1:M,spectrum_fft_t[j,:]); p=2)
    end
end

# writedlm("./Applications/Kuramoto model/system with switch/results/emd1.csv",emd1)
# writedlm("./Applications/Kuramoto model/system with switch/results/emd2.csv",emd2)
# writedlm("./Applications/Kuramoto model/system with switch/results/emd3.csv",emd3)
# writedlm("./Applications/Kuramoto model/system with switch/results/emd4.csv",emd4)

# writedlm("./Applications/Kuramoto model/system with switch/results/emd1t.csv",emd1t)
# writedlm("./Applications/Kuramoto model/system with switch/results/emd2t.csv",emd2t)
# writedlm("./Applications/Kuramoto model/system with switch/results/emd3t.csv",emd3t)
# writedlm("./Applications/Kuramoto model/system with switch/results/emd4t.csv",emd4t)

# writedlm("./Applications/Kuramoto model/system with switch/results/emd_fft.csv",emd_fft)
# writedlm("./Applications/Kuramoto model/system with switch/results/emd_fft_t.csv",emd_fft_t)

writedlm("./Applications/Kuramoto model/system with switch/results/emd1t_no_noise.csv",emd1t)
writedlm("./Applications/Kuramoto model/system with switch/results/emd2t_no_noise.csv",emd2t)
writedlm("./Applications/Kuramoto model/system with switch/results/emd3t_no_noise.csv",emd3t)
writedlm("./Applications/Kuramoto model/system with switch/results/emd4t_no_noise.csv",emd4t)

writedlm("./Applications/Kuramoto model/system with switch/results/emd_fft_no_noise.csv",emd_fft)
writedlm("./Applications/Kuramoto model/system with switch/results/emd_fft_t_no_noise.csv",emd_fft_t)
