using tau_recurrence
using DelimitedFiles
using DynamicalSystemsBase
using DelayEmbeddings
using RecurrenceAnalysis
using PyPlot
using Statistics
using DSP
pygui(true)

dt = 0.01
N = 6000
transients = 2000
lo = Systems.lorenz()
data = trajectory(lo, dt*N; dt = dt, Ttr = transients*dt)
data = regularize(data)

# compute RP and τ_rr
ε = 0.08

Y, tau_pec, _, _, _ = pecuzal_embedding(data[:,1]; τs = 0:50, w = estimate_delay(data[:,1], "mi_min"))

RP = RecurrenceAnalysis.RecurrenceMatrix(Y[1:1000,:], ε; fixedrate = true)
τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
τ_rr /= maximum(τ_rr)

Rg = grayscale(RP)

# FT of x and of Y
perio = DSP.Periodograms.periodogram(data[:,1]; fs = 1/dt)
perio2 = DSP.Periodograms.periodogram(τ_rr; fs = 1/dt)


# save variables and plot everything in Matlab
writedlm("./application/simple example/tau_rr.csv", τ_rr)
writedlm("./application/simple example/data.csv", data)
writedlm("./application/simple example/Y.csv", Y)
writedlm("./application/simple example/tau_pec.csv", tau_pec)
writedlm("./application/simple example/RP.csv", Rg)
writedlm("./application/simple example/perio_freq.csv", perio.freq)
writedlm("./application/simple example/perio_power.csv", perio.power)
writedlm("./application/simple example/perio_freq2.csv", perio2.freq)
writedlm("./application/simple example/perio_power2.csv", perio2.power)
