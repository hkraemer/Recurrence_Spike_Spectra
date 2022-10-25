using OptimalTransport
using Distances
using Distributions
using DelimitedFiles


num_levels = 51; # number of noise levels
thresholds = [0.9; 0.95; 0.99]; # regularization thresholds
N = 500; # length of the spectra 
# the sampling times for the different setups (limit-2, limit-3, chaos)
dt1 = 0.1
dt2 = 0.1
dt3 = 0.05
supports = [dt1:dt1:(N*dt1) dt2:dt2:(N*dt2) dt3:dt3:(N*dt3)];

# preallocate matrix storing all distances
dist_kl_div = zeros(length(thresholds), num_levels-1, 3)
dist_cosine = zeros(length(thresholds), num_levels-1, 3)
dist_hellinger = zeros(length(thresholds), num_levels-1, 3)
dist_wasserstein = zeros(length(thresholds), num_levels-1, 3)

# compute distances to noise-free spectrum
for (i, rho) in enumerate(thresholds)
    println("rho is $rho:")
    for system = 1:3
        println("system is $system:")
        # noise-free spectrum
        # load the computed data (computed in `compute_roessler_spike_spectra_for_different_noise_levels.jl`)
        a = vec(readdlm("./Applications/Roessler Example/computed data/tauRR_spectrum$(system)_roessler_$(rho)_1.csv"));
        # normalization
        a = a./sum(a);
        # for wasserstein distance make a distribution
        a_ = DiscreteNonParametric(supports[:, system], a)

        for k = 2:num_levels

            # load the computed data (computed in `compute_roessler_spike_spectra_for_different_noise_levels.jl`)
            b = vec(readdlm("./Applications/Roessler Example/computed data/tauRR_spectrum$(system)_roessler_$(rho)_$k.csv"));
            # normalization
            b = b./sum(b);
            # for wasserstein distance make a distribution
            b_ = DiscreteNonParametric(supports[:, system], b)

            # compute the distances
            dist_kl_div[i, k-1, system] = kl_divergence(a, b)
            dist_cosine[i, k-1, system] = cosine_dist(a, b)
            dist_hellinger[i, k-1, system] = hellinger(a, b)
            dist_wasserstein[i, k-1, system] = OptimalTransport.wasserstein(a_, b_)
        end
    end
end

# save the results
for (i, rho) in enumerate(thresholds)
    for system = 1:3
        writedlm("./Applications/Roessler Example/computed data/dist_kl_div_rho_$(rho)_system_$(system).csv", vec(dist_kl_div[i, :, system]))
        writedlm("./Applications/Roessler Example/computed data/dist_cosine_rho_$(rho)_system_$(system).csv", vec(dist_cosine[i, :, system]))
        writedlm("./Applications/Roessler Example/computed data/dist_hellinger_rho_$(rho)_system_$(system).csv", vec(dist_hellinger[i, :, system]))
        writedlm("./Applications/Roessler Example/computed data/dist_wasserstein_rho_$(rho)_system_$(system).csv", vec(dist_wasserstein[i, :, system]))
    end
end