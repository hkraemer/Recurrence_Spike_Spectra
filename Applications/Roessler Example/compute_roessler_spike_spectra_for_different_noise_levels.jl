import Pkg
Pkg.activate(".")
using InterSpikeSpectra
using DelimitedFiles
using DSP

# load the computed data (computed in `compute_roessler_trajectories_and_tau_rr_for_different_noise_levels.m`)
τ_rrs_1 = readdlm("./Applications/Roessler Example/computed data/tau_rr_1_all.csv")
τ_rrs_2 = readdlm("./Applications/Roessler Example/computed data/tau_rr_2_all.csv")
τ_rrs_3 = readdlm("./Applications/Roessler Example/computed data/tau_rr_3_all.csv")

N, M = size(τ_rrs_1)

# the sampling times for the different setups (limit-2, limit-3, chaos)
dt1 = 0.1
dt2 = 0.1
dt3 = 0.05

# compute power spectra
perio_1_freqs = zeros(N, Int(floor(M/2)+1))
perio_2_freqs = zeros(N, Int(floor(M/2)+1))
perio_3_freqs = zeros(N, Int(floor(M/2)+1))
perio_1_power = zeros(N, Int(floor(M/2)+1))
perio_2_power = zeros(N, Int(floor(M/2)+1))
perio_3_power = zeros(N, Int(floor(M/2)+1))

for i = 1:N
    print(i)
    perio1 = DSP.Periodograms.periodogram(vec(τ_rrs_1[i,:]); fs = 1/dt1)
    perio2 = DSP.Periodograms.periodogram(vec(τ_rrs_2[i,:]); fs = 1/dt2)
    perio3 = DSP.Periodograms.periodogram(vec(τ_rrs_3[i,:]); fs = 1/dt3)

    perio_1_freqs[i,:] = perio1.freq
    perio_1_power[i,:] = perio1.power
    perio_2_freqs[i,:] = perio2.freq
    perio_2_power[i,:] = perio2.power
    perio_3_freqs[i,:] = perio3.freq
    perio_3_power[i,:] = perio3.power
end


# compute spike spectra for different thresholds
thresholds = [0.9 0.95 0.99]
spectrum1 = zeros(length(thresholds), N, Int(M/2))
spectrum2 = zeros(length(thresholds), N, Int(M/2))
spectrum3 = zeros(length(thresholds), N, Int(M/2))

for (i, ρ) in enumerate(thresholds)
    println("This is for ρ: $ρ")
    for k = 1:N  
        println(k)
        spectrum1[i, k, :], _ = inter_spike_spectrum(vec(τ_rrs_1[i,:]); ρ_thres = ρ)
        spectrum2[i, k, :], _ = inter_spike_spectrum(vec(τ_rrs_2[i,:]); ρ_thres = ρ)
        spectrum3[i, k, :], _ = inter_spike_spectrum(vec(τ_rrs_3[i,:]); ρ_thres = ρ)
    end
end

# Save data
writedlm("./Applications/Roessler Example/computed data/perio_1_freqs.csv", perio_1_freqs)
writedlm("./Applications/Roessler Example/computed data/perio_2_freqs.csv", perio_2_freqs)
writedlm("./Applications/Roessler Example/computed data/perio_3_freqs.csv", perio_3_freqs)
writedlm("./Applications/Roessler Example/computed data/perio_1_power.csv", perio_1_power)
writedlm("./Applications/Roessler Example/computed data/perio_2_power.csv", perio_2_power)
writedlm("./Applications/Roessler Example/computed data/perio_3_power.csv", perio_3_power)


for (i, ρ) in enumerate(thresholds)
    for k = 1:N
        writedlm("./Applications/Roessler Example/computed data/tauRR_spectrum1_roessler_($ρ)_$k.csv", spectrum1[i, k, :])
        writedlm("./Applications/Roessler Example/computed data/tauRR_spectrum2_roessler_($ρ)_$k.csv", spectrum2[i, k, :])
        writedlm("./Applications/Roessler Example/computed data/tauRR_spectrum3_roessler_($ρ)_$k.csv", spectrum3[i, k, :])
    end
end
