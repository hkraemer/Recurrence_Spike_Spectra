using OptimalTransport
using Distances
using Distributions
using DelimitedFiles


num_levels = 51; # number of noise levels
sigmas = range(0, 0.5, num_levels); # noise levels
thresholds = [0.9; 0.95; 0.99]; # regularization thresholds
N = 500; # length of the spectra 
# the sampling times for the different setups (limit-2, limit-3, chaos)
dt1 = 0.1
dt2 = 0.1
dt3 = 0.05
support1 = dt1:dt1:(N*dt1)
support2 = dt2:dt2:(N*dt2)
support3 = dt3:dt3:(N*dt3)

k = 1;
rho = thresholds[1]
a = vec(readdlm("./Applications/Roessler Example/computed data/tauRR_spectrum1_roessler_($(rho))_$k.csv"));
k = 2;
b = vec(readdlm("./Applications/Roessler Example/computed data/tauRR_spectrum1_roessler_($(rho))_$k.csv"));

# normalization
a = DiscreteNonParametric(sigmas, a./sum(a))
b = DiscreteNonParametric(sigmas, b./sum(b))

dist = OptimalTransport.wasserstein(a, b)