
using DelayEmbeddings
using RecurrenceAnalysis
using DelimitedFiles
using PyPlot
using Statistics
pygui(true)

# Parameter-values for Logistic map
rs = 3.4:0.001:4

# delay for embedding
τ = 1
# recurrence threshold
ε = 0.05
# length of time series
N = 201

# load results
nsp = vec(readdlm("./Applications/Logistic Map Example/Results/results_Logistic_N_201_nsp.csv"))   # H₀: Random
nsp2 = vec(readdlm("./Applications/Logistic Map Example/Results/results_Logistic_N_201_nsp2.csv")) # H₀: iAAFT surrogates

for (i,r) in enumerate(rs)
    println(i)
    # compute τ-RR
    s = tau_recurrence.logistic_map_time_series(N, r)
    Y = DelayEmbeddings.hcat_lagged_values(s, s, τ)
    RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε)
    τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
    τ_rr /= maximum(τ_rr)

    # compute iAAFT-surrogates
    surros2 = tau_recurrence.compute_surrogate_τ_RR_iAAFT(100, s, ε, [0, τ])
    spec_surros2 = tau_recurrence.compute_surrogate_spectra(surros2) # compute spectra of these surrogates τ-RR-time series
    upper2, _ = tau_recurrence.compute_percentiles_of_surrogate_spectra(spec_surros)

    # spectrum of the τ-RR
    spectrum = tau_recurrence.get_spike_spectrum(τ_rr)
    # Maxima of the spectrum
    _, idx = tau_recurrence.get_maxima(spectrum)
    # compute number of significant peaks
    significant_peaks = spectrum[idx] .> upper[idx]
    significant_peaks2 = spectrum[idx] .> upper2[idx]
    number_of_significant_peaks[i] = sum(significant_peaks)
    number_of_significant_peaks2[i] = sum(significant_peaks2)
end

# compute Lyapunov exponent
lyap = zeros(length(rs))
for i = 1:length(rs)
    s = tau_recurrence.logistic_map_time_series(N, rs[i])
    lyap[i] = mean(log.(abs.(rs[i] .- 2*rs[i].*s)))
end

figure()
subplot(211)
plot(rs,lyap)
ylabel("λ")
title("Lyapunov exponent and τ-recurrence rate spectrum of Logistic map")
grid()
subplot(212)
plot(rs, number_of_significant_peaks)
xlabel("control parameter r")
ylabel("no. of significant peaks in the spectrum")
grid()


r = 3.57
s = tau_recurrence.logistic_map_time_series(N, r)
Y = DelayEmbeddings.hcat_lagged_values(s, s, τ)
RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε)
τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
τ_rr /= maximum(τ_rr)

# spectrum of the τ-RR
spectrum = tau_recurrence.get_spike_spectrum(τ_rr)

figure()
subplot(211)
title("τ-recurrence rate (r = $r)")
ylabel("τ-RR")
plot(τ_rr)
grid()
subplot(212)
title("Spectrum of τ_rr")
ylabel("Power")
xlabel("time [a.u.]")
plot(spectrum)
plot(upper)
grid()
subplots_adjust(hspace=.8)
