using DynamicalSystemsBase
using InterSpikeSpectra
using RecurrenceAnalysis
using DelimitedFiles
using DelayEmbeddings
using Random

using PyPlot
pygui(true)

## We illustrate the Inter Spike Spectrum for the paradigmatic Roessler system

dt = 0.05
N = 1000
transients = 5000
ε = 0.05
NN = 5000
σ = 0.05

# Period 2
roe = Systems.roessler(a = 0.36, b = 2, c = 4)
data1 = trajectory(roe, dt*NN; dt = dt, Ttr = transients*dt)
data1 = standardize(Dataset(data1))

RP1 = RecurrenceAnalysis.RecurrenceMatrix(data1, ε; fixedrate = true)
τ_rr1 = RecurrenceAnalysis.tau_recurrence(RP1)
τ_rr1 = τ_rr1[1:N] ./ maximum(τ_rr1[1:N])


Rg = grayscale(RP1)
imshow(Rg, cmap = "binary_r", extent = (1, size(RP1)[1], 1, size(RP1)[2]))
fig.tight_layout(pad=0.3)

figure()
scatter(1:N,τ_rr1)
grid()

spectrum1 = inter_spike_spectrum(τ_rr1)

# Period 3
u0=[.7, -1, 0.4]
roe = Systems.roessler(u0; a = 0.41, b = 2, c = 4)
data2 = trajectory(roe, dt*NN; dt = dt, Ttr = transients*dt)
data2 = standardize(Dataset(data2))

RP2 = RecurrenceAnalysis.RecurrenceMatrix(data2, ε; fixedrate = true)
Rg2 = grayscale(RP2)
τ_rr2 = RecurrenceAnalysis.tau_recurrence(RP2)
τ_rr2 = τ_rr2[1:N] ./ maximum(τ_rr2[1:N])

spectrum2 = inter_spike_spectrum(τ_rr2)

# Chaotic
roe = Systems.roessler(a = 0.428, b = 2, c = 4)
data3 = trajectory(roe, dt*NN; dt = dt, Ttr = transients*dt)
data3 = standardize(Dataset(data3))

RP3 = RecurrenceAnalysis.RecurrenceMatrix(data3, ε; fixedrate = true)
τ_rr3 = RecurrenceAnalysis.tau_recurrence(RP3)
τ_rr3 = τ_rr3[1:N] ./ maximum(τ_rr3[1:N])

spectrum3 = inter_spike_spectrum(τ_rr3)

## noise

# Period 2
data1_n = Matrix(data1) .+ σ.*randn(size(data1))

RP1_n = RecurrenceAnalysis.RecurrenceMatrix(data1_n, ε; fixedrate = true)
τ_rr1_n = RecurrenceAnalysis.tau_recurrence(RP1_n)
τ_rr1_n = τ_rr1_n[1:N] ./ maximum(τ_rr1_n[1:N])

spectrum1_n = inter_spike_spectrum(τ_rr1_n)

# Period 3
data2_n = Matrix(data2) .+ σ.*randn(size(data2))

RP2_n = RecurrenceAnalysis.RecurrenceMatrix(data2_n, ε; fixedrate = true)
τ_rr2_n = RecurrenceAnalysis.tau_recurrence(RP2_n)
τ_rr2_n = τ_rr2_n[1:N] ./ maximum(τ_rr2_n[1:N])

spectrum2_n = inter_spike_spectrum(τ_rr2_n)

# Chaotic
data3_n = Matrix(data3) .+ σ.*randn(size(data3))

RP3_n = RecurrenceAnalysis.RecurrenceMatrix(data3_n, ε; fixedrate = true)
τ_rr3_n = RecurrenceAnalysis.tau_recurrence(RP3_n)
τ_rr3_n = τ_rr3_n[1:N] ./ maximum(τ_rr3_n[1:N])

spectrum3_n = inter_spike_spectrum(τ_rr3_n)

Rg1 = grayscale(RP1[1:1000,1:1000])
Rg2 = grayscale(RP2[1:1000,1:1000])
Rg3 = grayscale(RP3[1:1000,1:1000])
Rg1_n = grayscale(RP1_n[1:1000,1:1000])
Rg2_n = grayscale(RP2_n[1:1000,1:1000])
Rg3_n = grayscale(RP3_n[1:1000,1:1000])

writedlm("./Applications/Roessler example/computed data/RP_1.csv", Rg1)
writedlm("./Applications/Roessler example/computed data/RP_2.csv", Rg2)
writedlm("./Applications/Roessler example/computed data/RP_3.csv", Rg3)
writedlm("./Applications/Roessler example/computed data/RP_1_n.csv", Rg1_n)
writedlm("./Applications/Roessler example/computed data/RP_2_n.csv", Rg2_n)
writedlm("./Applications/Roessler example/computed data/RP_3_n.csv", Rg3_n)

writedlm("./Applications/Roessler example/computed data/Y_1.csv", data1)
writedlm("./Applications/Roessler example/computed data/Y_2.csv", data2)
writedlm("./Applications/Roessler example/computed data/Y_3.csv", data3)
writedlm("./Applications/Roessler example/computed data/Y_1_n.csv", data1_n)
writedlm("./Applications/Roessler example/computed data/Y_2_n.csv", data2_n)
writedlm("./Applications/Roessler example/computed data/Y_3_n.csv", data3_n)

writedlm("./Applications/Roessler example/computed data/spectrum1.csv", spectrum1)
writedlm("./Applications/Roessler example/computed data/spectrum2.csv", spectrum2)
writedlm("./Applications/Roessler example/computed data/spectrum3.csv", spectrum3)
writedlm("./Applications/Roessler example/computed data/spectrum1_n.csv", spectrum1_n)
writedlm("./Applications/Roessler example/computed data/spectrum2_n.csv", spectrum2_n)
writedlm("./Applications/Roessler example/computed data/spectrum3_n.csv", spectrum3_n)

writedlm("./Applications/Roessler example/computed data/tau_rr_1.csv", τ_rr1)
writedlm("./Applications/Roessler example/computed data/tau_rr_2.csv", τ_rr2)
writedlm("./Applications/Roessler example/computed data/tau_rr_3.csv", τ_rr3)
writedlm("./Applications/Roessler example/computed data/tau_rr_1_n.csv", τ_rr1_n)
writedlm("./Applications/Roessler example/computed data/tau_rr_2_n.csv", τ_rr2_n)
writedlm("./Applications/Roessler example/computed data/tau_rr_3_n.csv", τ_rr3_n)


## Plot results for checking

figure(figsize=(15,10))
plot3D(data2[:,1],data2[:,2],data2[:,3])
title("DATA 2")
grid()

figure()
subplot(211)
plot(τ_rr2)
title("tauRR DATA 2")
grid()
subplot(212)
plot(spectrum2)
title("spectrum DATA 2")
grid()


figure(figsize=(15,10))
plot3D(data1[:,1],data1[:,2],data1[:,3])
title("DATA 1")
grid()

figure()
subplot(211)
plot(τ_rr1)
title("tauRR DATA 1")
grid()
subplot(212)
plot(spectrum1)
title("spectrum DATA 1")
grid()


##
# plot
figure(figsize=(15,10))

subplot(331)
Rg1 = grayscale(RP1)
imshow(Rg1, cmap = "binary_r", extent = (1, size(RP1)[1], 1, size(RP1)[2]))
subplot(332)
Rg2 = grayscale(RP2)
imshow(Rg2, cmap = "binary_r", extent = (1, size(RP2)[1], 1, size(RP2)[2]))
subplot(333)
Rg3 = grayscale(RP3)
imshow(Rg3, cmap = "binary_r", extent = (1, size(RP3)[1], 1, size(RP3)[2]))

subplot(334)
plot((1:length(spectrum1)),τ_rr1)
grid()
subplot(335)
plot((1:length(spectrum2)),τ_rr2)
grid()
subplot(336)
plot((1:length(spectrum3)),τ_rr3)
grid()
subplot(337)
plot((1:length(spectrum1)),spectrum1.^2)
yscale("log")
grid()
subplot(338)
plot((1:length(spectrum2)),spectrum2.^2)
yscale("log")
grid()
subplot(339)
plot((1:length(spectrum3)),spectrum3.^2)
yscale("log")
grid()


figure(figsize=(20,8))
ax1 = subplot()
ax1.plot(as, lyap, linewidth=2)
hlines(0, as[1], as[end], linestyle="dashed")
xlabel("control parameter r")
ylabel("λ₁")
grid()
xlim([as[1], as[end]])
ax2 = ax1.twinx()
grid()
roe = Systems.roessler(a = 0.36, b = 2, c = 4)
for a in as
  set_parameter!(roe, 1, a)
  #data = trajectory(roe, dt*Nplot; dt = dt, Ttr = transients*dt)
  plane = (1, 0.0)
  psos = poincaresos(roe, plane, 2000; Ttr = 1000)
  #plot(a*ones(length(psos)), psos[:,2], "k.", markersize=.1)
  ax2.scatter(a*ones(length(psos)), psos[:,2], s=.1, c="gray", alpha=0.2)
end
PyPlot.yticks([])
title("Bifurcation diagram and λ₁ of Roessler system")
grid()
tight_layout()
