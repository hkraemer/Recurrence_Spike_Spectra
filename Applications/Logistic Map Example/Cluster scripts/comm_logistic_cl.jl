using ClusterManagers
using Distributed
@everywhere N_tasks = parse(Int, ARGS[1])
@everywhere N_worker = N_tasks
addprocs(SlurmManager(N_worker))

@everywhere begin
    using InterSpikeSpectra
    using DelayEmbeddings
    using RecurrenceAnalysis
    using DelimitedFiles

    include("handy_functions.jl")

    # Parameter-values for Logistic map
    rs = 3.4:0.001:4

    # delay for embedding
    τ = 1
    # recurrence threshold
    ε = 0.05
    # length of time series
    N = 201

    # compute surrogates and their spectra
    surros = compute_surrogate_τ_RR(1000, .05, N-τ) # surrogate τ-RR-time series
    spec_surros = compute_surrogate_spectra(surros) # compute spectra of these surrogates τ-RR-time series
    upper, lower = compute_percentiles_of_surrogate_spectra(spec_surros)

    params = tuple(N, ε, rs)

end

@time begin
# loop over different r's
results = @distributed (vcat) for i in eachindex(rs)

    r = rs[i]
    # compute τ-RR
    s = logistic_map_time_series(N, r)
    Y = DelayEmbeddings.hcat_lagged_values(s, s, τ)
    RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε)
    τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
    τ_rr /= maximum(τ_rr)

    # compute iAAFT-surrogates
    surros2 = compute_surrogate_τ_RR_iAAFT(100, s, ε, [0, τ])
    spec_surros2 = compute_surrogate_spectra(surros2) # compute spectra of these surrogates τ-RR-time series
    upper2, _ = compute_percentiles_of_surrogate_spectra(spec_surros2)

    # spectrum of the τ-RR
    spectrum, _ = inter_spike_spectrum(τ_rr)
    # Maxima of the spectrum
    _, idx = get_maxima(spectrum)
    # compute number of significant peaks
    significant_peaks = spectrum[idx] .> upper[idx]
    significant_peaks2 = spectrum[idx] .> upper2[idx]
    number_of_significant_peaks = sum(significant_peaks)
    number_of_significant_peaks2 = sum(significant_peaks2)

    # Output
    tuple(number_of_significant_peaks, number_of_significant_peaks2, upper2)
end

end

writedlm("results_Logistic_N_$(N)_params.csv", params)
writedlm("spec_surros.csv", spec_surros)
writedlm("upper.csv", upper)

varnames = ["nsp", "nsp2", "upper2"]

for i = 1:length(varnames)
    writestr = "results_Logistic_N_$(N)_"*varnames[i]*".csv"
    data = []
    for j = 1:length(results)
        push!(data,results[j][i])
        writedlm(writestr, data)
    end
end
