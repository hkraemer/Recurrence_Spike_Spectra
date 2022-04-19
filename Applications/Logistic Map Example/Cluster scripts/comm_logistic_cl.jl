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
    using Random

    include("handy_functions.jl")

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


    number_of_significant_peaks1 = sum(significant_peaks1)
    number_of_significant_peaks2 = sum(significant_peaks2)
    number_of_significant_peaks3 = sum(significant_peaks3)

    number_of_significant_peaks_iafft1 = sum(significant_peaks_iafft1)
    number_of_significant_peaks_iafft2 = sum(significant_peaks_iafft2)
    number_of_significant_peaks_iafft3 = sum(significant_peaks_iafft3)

    # Output
    tuple(number_of_significant_peaks1, number_of_significant_peaks2, number_of_significant_peaks3, number_of_significant_peaks_iafft1, 
    number_of_significant_peaks_iafft2, number_of_significant_peaks_iafft3, upper_iafft1, upper_iafft2, upper_iafft3)
end

end

writedlm("results_Logistic_N_$(N)_params.csv", params)
writedlm("results_Logistic_N_$(N)_thres_$(thresholds[1])_spec_surros.csv", spec_surros1)
writedlm("results_Logistic_N_$(N)_thres_$(thresholds[2])_spec_surros.csv", spec_surros2)
writedlm("results_Logistic_N_$(N)_thres_$(thresholds[3])_spec_surros.csv", spec_surros3)
writedlm("results_Logistic_N_$(N)_thres_$(thresholds[1])_upper.csv", upper1)
writedlm("results_Logistic_N_$(N)_thres_$(thresholds[2])_upper.csv", upper2)
writedlm("results_Logistic_N_$(N)_thres_$(thresholds[3])_upper.csv", upper3)


varnames = ["nsp", "nsp_iafft"]

i = 1
writestr1 = "results_Logistic_N_$(N)_thres_$(thresholds[1])_"*varnames[i]*".csv"
writestr2 = "results_Logistic_N_$(N)_thres_$(thresholds[2])_"*varnames[i]*".csv"
writestr3 = "results_Logistic_N_$(N)_thres_$(thresholds[3])_"*varnames[i]*".csv"
i = 2
writestr4 = "results_Logistic_N_$(N)_thres_$(thresholds[1])_"*varnames[i]*".csv"
writestr5 = "results_Logistic_N_$(N)_thres_$(thresholds[2])_"*varnames[i]*".csv"
writestr6 = "results_Logistic_N_$(N)_thres_$(thresholds[3])_"*varnames[i]*".csv"

writestr7 = "results_Logistic_N_$(N)_thres_$(thresholds[1])_upper_iafft.csv"
writestr8 = "results_Logistic_N_$(N)_thres_$(thresholds[2])_upper_iafft.csv"
writestr9 = "results_Logistic_N_$(N)_thres_$(thresholds[3])_upper_iafft.csv"


data1 = []
data2 = []
data3 = []
data4 = []
data5 = []
data6 = []
data7 = []
data8 = []
data9 = []
for j = 1:length(results)
    push!(data1,results[j][1])
    writedlm(writestr1, data1)
    push!(data2,results[j][2])
    writedlm(writestr2, data2)
    push!(data3,results[j][3])
    writedlm(writestr3, data3)
    push!(data4,results[j][4])
    writedlm(writestr4, data4)
    push!(data5,results[j][5])
    writedlm(writestr5, data5)
    push!(data6,results[j][6])
    writedlm(writestr6, data6)
    push!(data7,results[j][7])
    writedlm(writestr7, data7)
    push!(data8,results[j][8])
    writedlm(writestr8, data8)
    push!(data9,results[j][9])
    writedlm(writestr9, data9)
end

