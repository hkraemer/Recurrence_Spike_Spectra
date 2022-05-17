import Pkg
Pkg.activate(".")
using InterSpikeSpectra
using DelimitedFiles
using DSP

# load the computed data (computed in `compute_roessler_trajectories_and_tau_rr.jl`)
τ_rr1 = vec(readdlm("./Applications/Roessler Example/computed data/tau_rr_1.csv"))
τ_rr1_n = vec(readdlm("./Applications/Roessler Example/computed data/tau_rr_1_n.csv"))
τ_rr2 = vec(readdlm("./Applications/Roessler Example/computed data/tau_rr_2.csv"))
τ_rr2_n = vec(readdlm("./Applications/Roessler Example/computed data/tau_rr_2_n.csv"))
τ_rr3 = vec(readdlm("./Applications/Roessler Example/computed data/tau_rr_3.csv"))
τ_rr3_n = vec(readdlm("./Applications/Roessler Example/computed data/tau_rr_3_n.csv"))

# compute power spectra
dt = 0.1
perio1 = DSP.Periodograms.periodogram(τ_rr1; fs = 1/dt)
perio1_n = DSP.Periodograms.periodogram(τ_rr1_n; fs = 1/dt)
dt = 0.1
perio2 = DSP.Periodograms.periodogram(τ_rr2; fs = 1/dt)
perio2_n = DSP.Periodograms.periodogram(τ_rr2_n; fs = 1/dt)
dt = 0.05
perio3 = DSP.Periodograms.periodogram(τ_rr3; fs = 1/dt)
perio3_n = DSP.Periodograms.periodogram(τ_rr3_n; fs = 1/dt)

# compute spike spectra for different thresholds
N = length(τ_rr1)
thresholds = [0.9 0.95 0.99]
spectrum1 = zeros(length(thresholds), Int(N/2))
spectrum2 = zeros(length(thresholds), Int(N/2))
spectrum3 = zeros(length(thresholds), Int(N/2))
spectrum1_n = zeros(length(thresholds), Int(N/2))
spectrum2_n = zeros(length(thresholds), Int(N/2))
spectrum3_n = zeros(length(thresholds), Int(N/2))

for ρ in thresholds 
    println("This is for ρ: $ρ")
    println("s1")
    spectrum1, _ = inter_spike_spectrum(τ_rr1; ρ_thres = ρ)
    spectrum1_n, _ = inter_spike_spectrum(τ_rr1_n; ρ_thres = ρ)
    println("s2")
    spectrum2, _ = inter_spike_spectrum(τ_rr2; ρ_thres = ρ)
    spectrum2_n, _ = inter_spike_spectrum(τ_rr2_n; ρ_thres = ρ)
    println("s3")
    spectrum3, _ = inter_spike_spectrum(τ_rr3; ρ_thres = ρ)
    spectrum3_n, _ = inter_spike_spectrum(τ_rr3_n; ρ_thres = ρ)
end

# Save data
writedlm("./Applications/Roessler Example/computed data/perio1_roessler_freq.csv", perio1.freq)
writedlm("./Applications/Roessler Example/computed data/perio1_roessler_power.csv", perio1.power)
writedlm("./Applications/Roessler Example/computed data/perio1_n_roessler_freq.csv", perio1_n.freq)
writedlm("./Applications/Roessler Example/computed data/perio1_n_roessler_power.csv", perio1_n.power)

writedlm("./Applications/Roessler Example/computed data/perio2_roessler_freq.csv", perio2.freq)
writedlm("./Applications/Roessler Example/computed data/perio2_roessler_power.csv", perio2.power)
writedlm("./Applications/Roessler Example/computed data/perio2_n_roessler_freq.csv", perio2_n.freq)
writedlm("./Applications/Roessler Example/computed data/perio2_n_roessler_power.csv", perio2_n.power)

writedlm("./Applications/Roessler Example/computed data/perio3_roessler_freq.csv", perio3.freq)
writedlm("./Applications/Roessler Example/computed data/perio3_roessler_power.csv", perio3.power)
writedlm("./Applications/Roessler Example/computed data/perio3_n_roessler_freq.csv", perio3_n.freq)
writedlm("./Applications/Roessler Example/computed data/perio3_n_roessler_power.csv", perio3_n.power)

for (i, ρ) in enumerate(thresholds)
    writedlm("./Applications/Roessler Example/computed data/tauRR_spectrum1_roessler_$ρ.csv", spectrum1[i,:])
    writedlm("./Applications/Roessler Example/computed data/tauRR_spectrum1_n_roessler_$ρ.csv", spectrum1_n[i,:])
    writedlm("./Applications/Roessler Example/computed data/tauRR_spectrum2_roessler_$ρ.csv", spectrum2[i,:])
    writedlm("./Applications/Roessler Example/computed data/tauRR_spectrum2_n_roessler_$ρ.csv", spectrum2_n[i,:])
    writedlm("./Applications/Roessler Example/computed data/tauRR_spectrum3_roessler_$ρ.csv", spectrum3[i,:])
    writedlm("./Applications/Roessler Example/computed data/tauRR_spectrum3_n_roessler_$ρ.csv", spectrum3_n[i,:])
end
