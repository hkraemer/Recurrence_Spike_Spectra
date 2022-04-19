begin
    import Pkg
    Pkg.activate(".")
    using InterSpikeSpectra
    using DelayEmbeddings
    using RecurrenceAnalysis
    using DelimitedFiles
    using Random

    include("../../../Methods/handy_functions.jl")

    # Parameter-values for Logistic map
    rs = 3.4:0.001:4

    # delay for embedding
    τ = 1
    # recurrence threshold
    ε = 0.05
    # length of time series
    N = 201
    # threshold for ISS
    thresholds = [0.9, 0.95, 0.99]

    # compute surrogates and their spectra
    surros = compute_surrogate_τ_RR(1000, .05, N-τ) # surrogate τ-RR-time series
    spec_surros1 = compute_surrogate_spectra(surros; ρ_thres=thresholds[1], regression_type=InterSpikeSpectra.normal()) # compute spectra of these surrogates τ-RR-time series
    spec_surros2 = compute_surrogate_spectra(surros; ρ_thres=thresholds[2], regression_type=InterSpikeSpectra.normal())
    spec_surros3 = compute_surrogate_spectra(surros; ρ_thres=thresholds[3], regression_type=InterSpikeSpectra.normal())
    upper1, _ = compute_percentiles_of_surrogate_spectra(spec_surros1)
    upper2, _ = compute_percentiles_of_surrogate_spectra(spec_surros2)
    upper3, _ = compute_percentiles_of_surrogate_spectra(spec_surros3)

    params = tuple(N, ε, rs, thresholds)
end

number_of_significant_peaks1 = zeros(length(rs))
number_of_significant_peaks2 = zeros(length(rs))
number_of_significant_peaks3 = zeros(length(rs))

number_of_significant_peaks_iafft1 = zeros(length(rs))
number_of_significant_peaks_iafft2 = zeros(length(rs))
number_of_significant_peaks_iafft3 = zeros(length(rs))

# loop over different r's
for i in eachindex(rs)
    print("Run $i")
    r = rs[i]
    # compute τ-RR
    s = logistic_map_time_series(N, r)
    Y = DelayEmbeddings.hcat_lagged_values(s, s, τ)
    RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε)
    τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
    τ_rr /= maximum(τ_rr)

    # compute iAAFT-surrogates
    Random.seed!(123)
    surros_iafft = compute_surrogate_τ_RR_iAAFT(100, s, ε, [0, τ])

    spec_surros_iafft1 = compute_surrogate_spectra(surros_iafft; ρ_thres=thresholds[1], regression_type=InterSpikeSpectra.normal()) # compute spectra of these surrogates τ-RR-time series
    spec_surros_iafft2 = compute_surrogate_spectra(surros_iafft; ρ_thres=thresholds[2], regression_type=InterSpikeSpectra.normal())
    spec_surros_iafft3 = compute_surrogate_spectra(surros_iafft; ρ_thres=thresholds[3], regression_type=InterSpikeSpectra.normal())
    
    upper_iafft1, _ = compute_percentiles_of_surrogate_spectra(spec_surros_iafft1)
    upper_iafft2, _ = compute_percentiles_of_surrogate_spectra(spec_surros_iafft2)
    upper_iafft3, _ = compute_percentiles_of_surrogate_spectra(spec_surros_iafft3)

    # spectrum of the τ-RR
    tol = 1e-3
    spectrum1, _ = inter_spike_spectrum(τ_rr; ρ_thres = thresholds[1], tol, regression_type=InterSpikeSpectra.normal(), verbose = false)
    spectrum2, _ = inter_spike_spectrum(τ_rr; ρ_thres = thresholds[2], tol, regression_type=InterSpikeSpectra.normal(), verbose = false)
    spectrum3, _ = inter_spike_spectrum(τ_rr; ρ_thres = thresholds[3], tol, regression_type=InterSpikeSpectra.normal(), verbose = false)
    # Maxima of the spectrum
    _, idx1 = get_maxima(spectrum1)
    _, idx2 = get_maxima(spectrum2)
    _, idx3 = get_maxima(spectrum3)
    # compute number of significant peaks
    # simple H0
    significant_peaks1 = spectrum1[idx1] .> upper1[idx1]
    significant_peaks2 = spectrum2[idx2] .> upper2[idx2]
    significant_peaks3 = spectrum3[idx3] .> upper3[idx3]
    # simple H0 - iafft's
    significant_peaks_iafft1 = spectrum1[idx1] .> upper_iafft1[idx1]
    significant_peaks_iafft2 = spectrum2[idx2] .> upper_iafft2[idx2]
    significant_peaks_iafft3 = spectrum3[idx3] .> upper_iafft3[idx3]


    number_of_significant_peaks1[i] = sum(significant_peaks1)
    number_of_significant_peaks2[i] = sum(significant_peaks2)
    number_of_significant_peaks3[i] = sum(significant_peaks3)

    number_of_significant_peaks_iafft1[i] = sum(significant_peaks_iafft1)
    number_of_significant_peaks_iafft2[i] = sum(significant_peaks_iafft2)
    number_of_significant_peaks_iafft3[i] = sum(significant_peaks_iafft3)

end

varnames = ["nsp", "nsp_iafft"]

i = 1
writestr1 = "results_Logistic_N_$(N)_thres_$(thresholds[1])_"*varnames[i]*".csv"
writestr2 = "results_Logistic_N_$(N)_thres_$(thresholds[2])_"*varnames[i]*".csv"
writestr3 = "results_Logistic_N_$(N)_thres_$(thresholds[3])_"*varnames[i]*".csv"
i = 2
writestr4 = "results_Logistic_N_$(N)_thres_$(thresholds[1])_"*varnames[i]*".csv"
writestr5 = "results_Logistic_N_$(N)_thres_$(thresholds[2])_"*varnames[i]*".csv"
writestr6 = "results_Logistic_N_$(N)_thres_$(thresholds[3])_"*varnames[i]*".csv"

writedlm(writestr1, number_of_significant_peaks1)
writedlm(writestr2, number_of_significant_peaks2)
writedlm(writestr3, number_of_significant_peaks3)

writedlm(writestr4, number_of_significant_peaks_iafft1)
writedlm(writestr5, number_of_significant_peaks_iafft2)
writedlm(writestr6, number_of_significant_peaks_iafft3)