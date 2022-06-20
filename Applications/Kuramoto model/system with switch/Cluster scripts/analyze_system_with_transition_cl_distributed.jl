using ClusterManagers
using Distributed
@everywhere N_tasks = parse(Int, ARGS[1])
@everywhere N_worker = N_tasks
addprocs(SlurmManager(N_worker))

@everywhere begin
    using InterSpikeSpectra
    using DynamicalSystemsBase
    using RecurrenceAnalysis
    using DelimitedFiles
    using Random
    using FFTW
    using JLD2
    using FileIO
    using DelayEmbeddings
    using Distances

    f = JLD2.jldopen("res_transition.jld2", "r")
    phases = Dataset(read(f,"phases")')
    frequencies = Dataset(read(f,"frequencies")')

    # discard transients and downsampling
    frequencies = frequencies[501:end,:]
    phases = phases[501:end,:]
    N, _ = size(phases)

    t = 1:N # time vector

    # define a artificial measurement function
    function measurement_function(x, a, b, c, σ=0)
        return a*x + b*x^2 + c*x^3 + σ*randn()
    end

    # transform phases
    a = 0.001
    b = 0.0043
    c = 0.002
    σ = 0.05
    phases_transform = Dataset(measurement_function.(Matrix(phases), a, b, c, σ))

    # windowed Analysis params
    ε = 0.05
    rho1 = 0.8
    rho2 = 0.85
    rho3 = 0.9
    rho4 = 0.95
    M = 500
    window = 600
    # regression method for the inter spike spectrum
    method = InterSpikeSpectra.lasso()

end

@time begin
# loop over different windows
results = @distributed (vcat) for i = 1:N-window
#results = @distributed (vcat) for i = 1:10

    Y = phases[i:i+window,:]
    Y2 = phases_transform[i:i+window,:]

    RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε; metric=Distances.PeriodicEuclidean(2*pi*ones(size(Y,2))), fixedrate = true)
    τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
    τ_rr = τ_rr[1:M] ./ maximum(τ_rr[1:M])

    RP2 = RecurrenceAnalysis.RecurrenceMatrix(Y2, ε; metric=Distances.PeriodicEuclidean(2*pi*ones(size(Y2,2))), fixedrate = true)
    τ_rr2 = RecurrenceAnalysis.tau_recurrence(RP2)
    τ_rr2 = τ_rr2[1:M] ./ maximum(τ_rr2[1:M])

    # Inter Spike Spectrum
    spectrum1, rhos1 = inter_spike_spectrum(τ_rr; ρ_thres = rho1, method, regression_type=normal())
    spectrum2, rhos2 = inter_spike_spectrum(τ_rr; ρ_thres = rho2, method, regression_type=normal())
    spectrum3, rhos3 = inter_spike_spectrum(τ_rr; ρ_thres = rho3, method, regression_type=normal())
    spectrum4, rhos4 = inter_spike_spectrum(τ_rr; ρ_thres = rho4, method, regression_type=normal())
    # transformed data
    spectrum1t, rhos1t = inter_spike_spectrum(τ_rr2; ρ_thres = rho1, method, regression_type=normal())
    spectrum2t, rhos2t = inter_spike_spectrum(τ_rr2; ρ_thres = rho2, method, regression_type=normal())
    spectrum3t, rhos3t = inter_spike_spectrum(τ_rr2; ρ_thres = rho3, method, regression_type=normal())
    spectrum4t, rhos4t = inter_spike_spectrum(τ_rr2; ρ_thres = rho4, method, regression_type=normal())

    # Fourier spectrum
    F1 = fftshift(fft(τ_rr)) 
    F2 = fftshift(fft(τ_rr2)) 
    FF1 = abs.(F1[Int(M/2)+1:end])
    FF2 = abs.(F2[Int(M/2)+1:end])

    # Output
    tuple(spectrum1, spectrum2, spectrum3, spectrum4, spectrum1t, spectrum2t, spectrum3t, spectrum4t, FF1, FF2, rhos1, rhos2, rhos3, rhos4, rhos1t, rhos2t, rhos3t, rhos4t)
end

end

# save the computed data

varnames = ["spectrum_0_8", "spectrum_0_85", "spectrum_0_9", "spectrum_0_95", "spectrum_0_8_t", "spectrum_0_85_t", "spectrum_0_9_t", "spectrum_0_95_t", "fft_spectrum", "fft_spectrum_t", "rhos1", "rhos2", "rhos3", "rhos4", "rhos1t", "rhos2t", "rhos3t", "rhos4t"]

for i = 1:length(varnames)
    writestr = "results_kuramoto_switch_"*varnames[i]*".csv"
    data = []
    for j = 1:length(results)
        push!(data,results[j][i])
        writedlm(writestr, data)
    end
end