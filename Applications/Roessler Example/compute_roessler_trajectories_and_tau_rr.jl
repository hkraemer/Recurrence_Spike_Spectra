import Pkg
Pkg.activate(".")

using DynamicalSystemsBase
using InterSpikeSpectra
using RecurrenceAnalysis
using DelimitedFiles
using DelayEmbeddings
using Random

using PyPlot
pygui(true)

## We illustrate the Inter Spike Spectrum for the paradigmatic Roessler system
N = 1000
transients = 5000
NN = 5000
σ = 0.05

# Period 2
u0=[.7, -1, 0.4]
dt = 0.05
ε = 0.1
roe = Systems.roessler(u0; a = 0.36, b = 2, c = 4)
data1 = trajectory(roe, dt*NN; Δt = dt, Ttr = transients*dt)
data1 = standardize(Dataset(data1))

RP1 = RecurrenceAnalysis.RecurrenceMatrix(data1, ε; fixedrate = true)
τ_rr1 = RecurrenceAnalysis.tau_recurrence(RP1)
τ_rr1 = τ_rr1[1:N] ./ maximum(τ_rr1[1:N])

# # plot data for checking
# begin
#     Rg = grayscale(RP1[1:N,1:N])
#     fig = imshow(Rg, cmap = "binary_r", extent = (1, size(RP1[1:N,1:N])[1], 1, size(RP1[1:N,1:N])[2]))
#     figure()
#     plot(1:N,τ_rr1)
#     grid()
#     figure(figsize=(15,10))
#     plot3D(data1[:,1],data1[:,2],data1[:,3])
#     title("DATA 2")
#     grid()
# end

# Period 3
u0=[.7, -1, 0.4]
dt = 0.05
ε = 0.1
roe = Systems.roessler(u0; a = 0.41, b = 2, c = 4)
data2 = trajectory(roe, dt*NN; Δt = dt, Ttr = transients*dt)
data2 = standardize(Dataset(data2))

RP2 = RecurrenceAnalysis.RecurrenceMatrix(data2, ε; fixedrate = true)
τ_rr2 = RecurrenceAnalysis.tau_recurrence(RP2)
τ_rr2 = τ_rr2[1:N] ./ maximum(τ_rr2[1:N])

# # plot data for checking
# begin
#     Rg = grayscale(RP2[1:N,1:N])
#     fig = imshow(Rg, cmap = "binary_r", extent = (1, size(RP2[1:N,1:N])[1], 1, size(RP2[1:N,1:N])[2]))
#     figure()
#     plot(1:N,τ_rr2)
#     grid()
#     figure(figsize=(15,10))
#     plot3D(data2[:,1],data2[:,2],data2[:,3])
#     title("DATA 2")
#     grid()
# end


# Chaos
u0=[-0.1242,-2.5415,0.2772]
dt = 0.1
ε = 0.05
roe = Systems.roessler(a = 0.428, b = 2, c = 4)
data3 = trajectory(roe, dt*NN; Δt = dt, Ttr = transients*dt)
data3 = standardize(Dataset(data3))

RP3 = RecurrenceAnalysis.RecurrenceMatrix(data3, ε; fixedrate = true)
τ_rr3 = RecurrenceAnalysis.tau_recurrence(RP3)
τ_rr3 = τ_rr3[1:N] ./ maximum(τ_rr3[1:N])

# # plot data for checking
# begin
#     Rg = grayscale(RP3[1:N,1:N])
#     fig = imshow(Rg, cmap = "binary_r", extent = (1, size(RP3[1:N,1:N])[1], 1, size(RP3[1:N,1:N])[2]))
#     figure()
#     plot(1:N,τ_rr3)
#     grid()
#     figure(figsize=(15,10))
#     plot3D(data3[:,1],data3[:,2],data3[:,3])
#     title("DATA 3")
#     grid()
# end


## Additive noise

# Period 2
ε = 0.1
data1_n = Matrix(data1) .+ σ.*randn(size(data1))
RP1_n = RecurrenceAnalysis.RecurrenceMatrix(data1_n, ε; fixedrate = true)
τ_rr1_n = RecurrenceAnalysis.tau_recurrence(RP1_n)
τ_rr1_n = τ_rr1_n[1:N] ./ maximum(τ_rr1_n[1:N])

# # plot data for checking
# begin
#     Rg = grayscale(RP1_n[1:N,1:N])
#     fig = imshow(Rg, cmap = "binary_r", extent = (1, size(RP3[1:N,1:N])[1], 1, size(RP3[1:N,1:N])[2]))
#     figure()
#     plot(1:N,τ_rr1_n)
#     grid()
#     figure(figsize=(15,10))
#     plot3D(data1_n[:,1],data1_n[:,2],data1_n[:,3])
#     grid()
# end

# Period 3
ε = 0.1
data2_n = Matrix(data2) .+ σ.*randn(size(data2))
RP2_n = RecurrenceAnalysis.RecurrenceMatrix(data2_n, ε; fixedrate = true)
τ_rr2_n = RecurrenceAnalysis.tau_recurrence(RP2_n)
τ_rr2_n = τ_rr2_n[1:N] ./ maximum(τ_rr2_n[1:N])

# # plot data for checking
# begin
#     Rg = grayscale(RP2_n[1:N,1:N])
#     fig = imshow(Rg, cmap = "binary_r", extent = (1, size(RP3[1:N,1:N])[1], 1, size(RP3[1:N,1:N])[2]))
#     figure()
#     plot(1:N,τ_rr2_n)
#     grid()
#     figure(figsize=(15,10))
#     plot3D(data2_n[:,1],data2_n[:,2],data2_n[:,3])
#     grid()
# end

# Chaotic
ε = 0.05
data3_n = Matrix(data3) .+ σ.*randn(size(data3))
RP3_n = RecurrenceAnalysis.RecurrenceMatrix(data3_n, ε; fixedrate = true)
τ_rr3_n = RecurrenceAnalysis.tau_recurrence(RP3_n)
τ_rr3_n = τ_rr3_n[1:N] ./ maximum(τ_rr3_n[1:N])

# # plot data for checking
# begin
#     Rg = grayscale(RP3_n[1:N,1:N])
#     fig = imshow(Rg, cmap = "binary_r", extent = (1, size(RP3[1:N,1:N])[1], 1, size(RP3[1:N,1:N])[2]))
#     figure()
#     plot(1:N,τ_rr3_n)
#     grid()
#     figure(figsize=(15,10))
#     plot3D(data3_n[:,1],data3_n[:,2],data3_n[:,3])
#     grid()
# end


## Save data
Rg1 = grayscale(RP1[1:N,1:N])
Rg2 = grayscale(RP2[1:N,1:N])
Rg3 = grayscale(RP3[1:N,1:N])
Rg1_n = grayscale(RP1_n[1:N,1:N])
Rg2_n = grayscale(RP2_n[1:N,1:N])
Rg3_n = grayscale(RP3_n[1:N,1:N])

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

writedlm("./Applications/Roessler example/computed data/tau_rr_1.csv", τ_rr1)
writedlm("./Applications/Roessler example/computed data/tau_rr_2.csv", τ_rr2)
writedlm("./Applications/Roessler example/computed data/tau_rr_3.csv", τ_rr3)
writedlm("./Applications/Roessler example/computed data/tau_rr_1_n.csv", τ_rr1_n)
writedlm("./Applications/Roessler example/computed data/tau_rr_2_n.csv", τ_rr2_n)
writedlm("./Applications/Roessler example/computed data/tau_rr_3_n.csv", τ_rr3_n)

# # bifurcation diagram

# figure(figsize=(20,8))
# ax1 = subplot()
# ax1.plot(as, lyap, linewidth=2)
# hlines(0, as[1], as[end], linestyle="dashed")
# xlabel("control parameter r")
# ylabel("λ₁")
# grid()
# xlim([as[1], as[end]])
# ax2 = ax1.twinx()
# grid()
# roe = Systems.roessler(a = 0.36, b = 2, c = 4)
# for a in as
#   set_parameter!(roe, 1, a)
#   #data = trajectory(roe, dt*Nplot; dt = dt, Ttr = transients*dt)
#   plane = (1, 0.0)
#   psos = poincaresos(roe, plane, 2000; Ttr = 1000)
#   #plot(a*ones(length(psos)), psos[:,2], "k.", markersize=.1)
#   ax2.scatter(a*ones(length(psos)), psos[:,2], s=.1, c="gray", alpha=0.2)
# end
# PyPlot.yticks([])
# title("Bifurcation diagram and λ₁ of Roessler system")
# grid()
# tight_layout()
