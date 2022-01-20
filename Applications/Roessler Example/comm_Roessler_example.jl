using tau_recurrence
using DelimitedFiles
using DynamicalSystemsBase
using PyPlot
using Statistics
using DSP
pygui(true)

# recurrence threshold
ε = 0.08
# length of time series
N = 6000
# time step
dt = 0.2
# transients
transients = 2000

a = 0.38
roe = Systems.roessler(a = a, b = 2, c = 4)
data = trajectory(roe, dt*N; dt = dt, Ttr = transients*dt)
data = regularize(data)

writedlm("data_roessler.csv", data)

τ_rr = vec(readdlm("./application/simple example/tauRR_limit_roessler.csv"))
τ_rr_d = vec(readdlm("./application/simple example/tauRR_skeleton_limit_roessler.csv"))

τ_rr = τ_rr./maximum(τ_rr)
τ_rr_d = τ_rr_d./maximum(τ_rr_d)

perio = DSP.Periodograms.periodogram(τ_rr[1:500]; fs = 1/dt)
perio_d = DSP.Periodograms.periodogram(τ_rr_d[1:500]; fs = 1/dt)

τ_rr_spectrum = tau_recurrence.get_spike_spectrum(τ_rr[1:500])
τ_rr_spectrum_d = tau_recurrence.get_spike_spectrum(τ_rr_d[1:500])

writedlm("perio_roessler_freq.csv", perio.freq)
writedlm("perio_roessler_power.csv", perio.power)
writedlm("perio_d_roessler_freq.csv", perio_d.freq)
writedlm("perio_d_roessler_power.csv", perio_d.power)

writedlm("tauRR_spectrum_roessler.csv", τ_rr_spectrum)
writedlm("tauRR_spectrum_d_roessler.csv", τ_rr_spectrum_d)
