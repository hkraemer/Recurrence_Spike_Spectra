using ClusterManagers
using Distributed
@everywhere N_tasks = parse(Int, ARGS[1])
@everywhere N_worker = N_tasks
addprocs(SlurmManager(N_worker))

@everywhere begin
    using InterSpikeSpectra
    using DelayEmbeddings
    using DynamicalSystemsBase
    using RecurrenceAnalysis
    using DelimitedFiles
    using Random
    using FFTW

    # Parameter-values for Roessler transitions
    as = 0.36:0.0001:0.43

    # recurrence threshold
    ε = 0.08
    # length of time series
    N = 1000
    # time step
    Δt = 0.2
    # according sampling frequency
    fs = 1/Δt
    # transients
    transients = 2000
    # size of the considered τ-RR
    tau_window = 500
    w2 = Int(tau_window/2)
    # make a time vector for the Fourier analysis
    tt = 0:Δt:(tau_window-1)*Δt
    # method for the inter spike spectrum regression_type
    method = InterSpikeSpectra.lasso()
    # rho-thresholds for the inter spike spectrum
    rho1 = 0.85
    rho2 = 0.9
    rho3 = 0.95
    rho4 = 0.99
    # init roessler
    roe = Systems.roessler(a = 0.36, b = 2, c = 4)

    freqs = fftshift(fftfreq(length(tt), fs)) # according frequencies
    f_freqs = collect(freqs[w2+1:end])
    t_freqs = [1/i for i in f_freqs]

    params = tuple(freqs, f_freqs, t_freqs, as, ε, Δt, N, rho1, rho2, rho3, rho4, )

end

@time begin
# loop over different r's
results = @distributed (vcat) for i in eachindex(as)

    Random.seed!(1234)
    a = as[i]
    # compute τ-RR
    set_parameter!(roe, 1, a)
    data = trajectory(roe, Δt*N; Δt, Ttr = transients*Δt)

    # reconstruction
    Y, tau_pec, _, _, _ = pecuzal_embedding(data[:,1]; τs = 0:50, w = estimate_delay(data[:,1], "mi_min"), verbose=false)
    RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε; fixedrate = true)
    τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
    τ_rr /= maximum(τ_rr)

    # true trajectory
    RP_true = RecurrenceAnalysis.RecurrenceMatrix(data, ε; fixedrate = true)
    τ_rr_true = RecurrenceAnalysis.tau_recurrence(RP_true)
    τ_rr_true /= maximum(τ_rr_true)

    # spectrum of the τ-RR
    spectrum1, _ = inter_spike_spectrum(τ_rr[1:tau_window]; ρ_thres = rho1, method, regression_type=normal())
    spectrum2, _ = inter_spike_spectrum(τ_rr[1:tau_window]; ρ_thres = rho2, method, regression_type=normal())
    spectrum3, _ = inter_spike_spectrum(τ_rr[1:tau_window]; ρ_thres = rho3, method, regression_type=normal())
    spectrum4, _ = inter_spike_spectrum(τ_rr[1:tau_window]; ρ_thres = rho4, method, regression_type=normal())

    spectrum1t, _ = inter_spike_spectrum(τ_rr_true[1:tau_window]; ρ_thres = rho1, method, regression_type=normal())
    spectrum2t, _ = inter_spike_spectrum(τ_rr_true[1:tau_window]; ρ_thres = rho2, method, regression_type=normal())
    spectrum3t, _ = inter_spike_spectrum(τ_rr_true[1:tau_window]; ρ_thres = rho3, method, regression_type=normal())
    spectrum4t, _ = inter_spike_spectrum(τ_rr_true[1:tau_window]; ρ_thres = rho4, method, regression_type=normal())

    # Fourier spectrum
    F1 = fftshift(fft(τ_rr[1:tau_window]))  # embedded RP
    F2 = fftshift(fft(τ_rr_true[1:tau_window])) # RP based on true trajectory
    F3 = fftshift(fft(data[1:tau_window,1])) # spectrum based on x-time series

    FF1 = abs.(F1[w2+1:end])
    FF2 = abs.(F2[w2+1:end])
    FF3 = abs.(F3[w2+1:end])

    # Output
    tuple(spectrum1, spectrum2, spectrum3, spectrum4, spectrum1t, spectrum2t, spectrum3t, spectrum4t, FF1, FF2, FF3)
end

end

# save the computed data
writedlm("results_Roessler_N_$(N)_params.csv", params)

varnames = ["ISS_0_85", "ISS_0_9", "ISS_0_95", "ISS_0_99", "ISS_0_85_true", "ISS_0_9_true", "ISS_0_95_true", "ISS_0_99_true", "FFT_tau_rr_recon", "FFT_tau_rr_true", "FFT_time_series"]

for i = 1:length(varnames)
    writestr = "results_Roessler_N_$(N)_"*varnames[i]*".csv"
    data = []
    for j = 1:length(results)
        push!(data,results[j][i])
        writedlm(writestr, data)
    end
end
