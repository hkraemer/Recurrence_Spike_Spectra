
begin
    using InterSpikeSpectra
    using Revise
    using DelayEmbeddings
    using RecurrenceAnalysis
    using DelimitedFiles
    using PyPlot
    pygui(true)

    include("./Methods/handy_functions.jl")

    # Parameter-values for Logistic map
    rs = 3.4:0.001:4

    # delay for embedding
    τ = 1
    # recurrence threshold
    ε = 0.05
    # length of time series
    N = 201
    # threshold for ISS
    threshold = 0.95

    # compute surrogates and their spectra
    surros = compute_surrogate_τ_RR(1000, .05, N-τ) # surrogate τ-RR-time series
    spec_surros = compute_surrogate_spectra(surros; ρ_thres=threshold, regression_type=InterSpikeSpectra.normal()) # compute spectra of these surrogates τ-RR-time series
    upper, lower = compute_percentiles_of_surrogate_spectra(spec_surros)

    params = tuple(N, ε, rs, threshold)

end



include("./Methods/handy_functions.jl")
i = 6
begin 
    r = rs[i]
    # compute τ-RR
    s = logistic_map_time_series(N, r)
    Y = DelayEmbeddings.hcat_lagged_values(s, s, τ)
    RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε)
    τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
    τ_rr /= maximum(τ_rr)

    # compute iAAFT-surrogates
    surros2 = compute_surrogate_τ_RR_iAAFT(100, s, ε, [0, τ])
    spec_surros2 = compute_surrogate_spectra(surros2; ρ_thres=threshold, regression_type=InterSpikeSpectra.normal())

    for row = 1:size(spec_surros2,1)
        if any(findall(x->isnan(x), spec_surros2[row,:]))
            print("THIS IS NAN:")
            print(row)
        end
    end
end


i = 57
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

for i = 1:100
    print(i)
    if findall(x->isnan(x), spec_surros_iafft3[i,:]) == []
        continue
    else
        nanas = findall(x->isnan(x), spec_surros_iafft3[i,:])
        break
    end
end

writedlm("test_series.csv",surros_iafft[59,:])


##

data = readdlm("test_series.csv")
