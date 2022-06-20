using DelimitedFiles
using Distances
using OptimalTransport
using Distributions
using DelayEmbeddings

## We compute the pairwise earth mover's distance of the computed spectra of 
## the Roessler model with transitions (computations have been carried 
## out on the cluster -windowed spectrum determination- in the script `comm_roessler_cl.jl`

# load the results:
NN = 1000
spectrum1 = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_ISS_0_85.csv")
spectrum2 = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_ISS_0_9.csv")
spectrum3 = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_ISS_0_95.csv")
spectrum4 = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_ISS_0_99.csv")

spectrum1t = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_ISS_0_85_true.csv")
spectrum2t = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_ISS_0_9_true.csv")
spectrum3t = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_ISS_0_95_true.csv")
spectrum4t = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_ISS_0_99_true.csv")

FFT1 = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_FFT_tau_rr_recon.csv")
FFT2 = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_FFT_tau_rr_true.csv")
FFT3 = readdlm("./Applications/Roessler Transition/results/results_Roessler_N_$(NN)_FFT_time_series.csv")

N, M = size(spectrum1)

# standardize FFT-spectra to probabilities
for j = 1:N
    FFT1[j,:] = FFT1[j,:] ./ sum(FFT1[j,:])
    FFT2[j,:] = FFT2[j,:] ./ sum(FFT2[j,:])
    FFT3[j,:] = FFT3[j,:] ./ sum(FFT3[j,:])
end


emd1 = zeros(N,N)
emd2 = zeros(N,N)
emd3 = zeros(N,N)
emd4 = zeros(N,N)
emd1t = zeros(N,N)
emd2t = zeros(N,N)
emd3t = zeros(N,N)
emd4t = zeros(N,N)
emd_fft1 = zeros(N,N)
emd_fft2 = zeros(N,N)
emd_fft3 = zeros(N,N)
@time begin
    Threads.@threads for i = 1:N
        for j = 1:N
            emd1[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum1[i,:]), DiscreteNonParametric(1:M,spectrum1[j,:]); p=2)
            emd2[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum2[i,:]), DiscreteNonParametric(1:M,spectrum2[j,:]); p=2)
            emd3[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum3[i,:]), DiscreteNonParametric(1:M,spectrum3[j,:]); p=2)
            emd4[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum4[i,:]), DiscreteNonParametric(1:M,spectrum4[j,:]); p=2)

            emd1t[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum1t[i,:]), DiscreteNonParametric(1:M,spectrum1t[j,:]); p=2)
            emd2t[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum2t[i,:]), DiscreteNonParametric(1:M,spectrum2t[j,:]); p=2)
            emd3t[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum3t[i,:]), DiscreteNonParametric(1:M,spectrum3t[j,:]); p=2)
            emd4t[i,j] = wasserstein(DiscreteNonParametric(1:M,spectrum4t[i,:]), DiscreteNonParametric(1:M,spectrum4t[j,:]); p=2)

            emd_fft1[i,j] = wasserstein(DiscreteNonParametric(1:M,FFT1[i,:]), DiscreteNonParametric(1:M,FFT1[j,:]); p=2)
            emd_fft2[i,j] = wasserstein(DiscreteNonParametric(1:M,FFT2[i,:]), DiscreteNonParametric(1:M,FFT2[j,:]); p=2)
            emd_fft3[i,j] = wasserstein(DiscreteNonParametric(1:M,FFT3[i,:]), DiscreteNonParametric(1:M,FFT3[j,:]); p=2)
        end
    end
end

writedlm("./Applications/Roessler Transition/results/emd1.csv",emd1)
writedlm("./Applications/Roessler Transition/results/emd2.csv",emd2)
writedlm("./Applications/Roessler Transition/results/emd3.csv",emd3)
writedlm("./Applications/Roessler Transition/results/emd4.csv",emd4)

writedlm("./Applications/Roessler Transition/results/emd1t.csv",emd1t)
writedlm("./Applications/Roessler Transition/results/emd2t.csv",emd2t)
writedlm("./Applications/Roessler Transition/results/emd3t.csv",emd3t)
writedlm("./Applications/Roessler Transition/results/emd4t.csv",emd4t)

writedlm("./Applications/Roessler Transition/results/emd_FFT1.csv",emd_fft1)
writedlm("./Applications/Roessler Transition/results/emd_FFT2.csv",emd_fft2)
writedlm("./Applications/Roessler Transition/results/emd_FFT3.csv",emd_fft3)
