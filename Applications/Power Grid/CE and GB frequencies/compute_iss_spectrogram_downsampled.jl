using DelayEmbeddings
using TreeEmbedding  # this will be merged into DelayEmbeddings soon
using DelimitedFiles
using Random
using InterSpikeSpectra
using RecurrenceAnalysis
using Statistics

## Downsampled frequency data
# f_ce = vec(readdlm("./Applications/Power Grid/CE and GB frequencies/data/frequencies_ce_downsampled.csv"))
# f_gb = vec(readdlm("./Applications/Power Grid/CE and GB frequencies/data/frequencies_gb_downsampled.csv"))

# f_ce = vec(readdlm("./Applications/Power Grid/CE and GB frequencies/data/frequencies_ce_downsampled_30s.csv"))
# f_gb = vec(readdlm("./Applications/Power Grid/CE and GB frequencies/data/frequencies_gb_downsampled_30s.csv"))

f_ce = vec(readdlm("./Applications/Power Grid/CE and GB frequencies/data/frequencies_ce_downsampled_20s.csv"))
f_gb = vec(readdlm("./Applications/Power Grid/CE and GB frequencies/data/frequencies_gb_downsampled_20s.csv"))

# T_ce = 60 # sampling time in secs
# T_gb = 60

# T_ce = 30 # sampling time in secs
# T_gb = 30

T_ce = 20 # sampling time in secs
T_gb = 20

# L_block_ce = 1440 # block length covering one full day (24 hours)
# L_block_gb = 1440

# L_block_ce = 2880 # block length covering one full day (24 hours)
# L_block_gb = 2880

L_block_ce = 4320 # block length covering one full day (24 hours)
L_block_gb = 4320

rho = 0.95
epsilon = 0.08
τs = 0:150
tau_window = 600

spec_ce, d_ce, taus_ce = iss_spec(f_ce, L_block_ce, τs, epsilon, tau_window, rho; embed_method="afnn")
avrg_spec_ce = [mean(spec_ce[i,:]) for i = 1:size(spec_ce,1)]

spec_gb, d_gb, taus_gb = iss_spec(f_gb, L_block_gb, τs, epsilon, tau_window, rho; embed_method="afnn")
avrg_spec_gb = [mean(spec_gb[i,:]) for i = 1:size(spec_gb,1)]

writedlm("./Applications/Power Grid/CE and GB frequencies/data/iss_spec_ce_standard_20_0_95_0_8.csv", avrg_spec_ce)
writedlm("./Applications/Power Grid/CE and GB frequencies/data/iss_spec_gb_standard_20_0_95_0_8.csv", avrg_spec_gb)

using PyPlot
pygui(true)

figure()
subplot(211)
plot(1:300, avrg_spec_ce[1:300])
grid()
subplot(212)
plot(1:300, avrg_spec_gb[1:300])
grid()


# ##
# method = InterSpikeSpectra.lasso()
# blocksize = 1440
# x = f_ce

# i = 10
# time = 1+(i-1)*blocksize:(i*blocksize)
# # embedding
# w = estimate_delay(x[time], "mi_min", τs)
# Y, τ, _ = optimal_traditional_de(x[time], "afnn", "mi_min"; τs, w)
# delays = τ
# dimensions = size(Y,2)
# # RP and tauRR
# RP = RecurrenceAnalysis.RecurrenceMatrix(Y, epsilon; fixedrate = true)
# τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
# τ_rr /= maximum(τ_rr)
# τ_rr[1:w].=0

# spec, _ = inter_spike_spectrum(τ_rr[1:tau_window]; ρ_thres = rho, method, regression_type=normal())


# figure()
# imshow(grayscale(RP))

# figure()
# subplot(211)
# plot(τ_rr[1:tau_window])
# grid()
# subplot(212)
# plot(spec)
# grid()



##

function iss_spec(x::Vector, blocksize::Integer, τs::AbstractRange, ε::Real, tau_window::Integer, rho::Real; method = InterSpikeSpectra.lasso(), embed_method::String="ifnn")
    
    @assert tau_window ≤ blocksize
    T = length(x)
    blocks = Int(floor(T/blocksize))
    
    spectrum = zeros(Int(tau_window/2), blocks)
    delays = zeros(Int, blocks)
    dimensions = zeros(Int, blocks)

    for i = 1:blocks
        time = 1+(i-1)*blocksize:(i*blocksize)
        # embedding
        w = estimate_delay(x[time], "mi_min", τs)
        #Y, τ, _ = optimal_traditional_de(x[time], embed_method, "mi_min"; τs, w)
        τ = w
        Y = embed(x[time], 5, τ)
        delays[i] = τ
        dimensions[i] = size(Y,2)
        # RP and tauRR
        RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε; fixedrate = true)
        τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
        τ_rr /= maximum(τ_rr)
        #τ_rr[1:w].=0 # theiler window
        # spectrum of the τ-RR
        spectrum[:,i], _ = inter_spike_spectrum(τ_rr[1:tau_window]; ρ_thres = rho, method, regression_type=normal())
    end
    return spectrum, delays, dimensions
end