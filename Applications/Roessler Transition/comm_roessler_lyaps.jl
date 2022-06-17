using DelayEmbeddings
using DynamicalSystemsBase
using ChaosTools
using DelimitedFiles
using PyPlot
using Statistics
pygui(true)

# Parameter-values for Roessler transitions
as = 0.36:0.0001:0.43

roe = Systems.roessler(a = 0.36, b = 2, c = 4) # init Roessler
# check Lyapunox spectrum
λ1 = zeros(length(as))
λs = zeros(3,length(as))
for (i,a) in enumerate(as)
  println(i)
  set_parameter!(roe, 1, a)
  λ1[i] = lyapunov(roe, 100000; Ttr = 10000)
  λs[i,:] = lyapunovspectrum(roe, 100000; Ttr = 10000)

end

writedlm("./Applications/Roessler Transition/results/Lyaps_Roessler_0_36_to_0_43.csv", λ1)

λ1 = readdlm("./Applications/Roessler Transition/results/Lyaps_Roessler_0_36_to_0_43.csv")
λ1 = vec(λ1[:])

figure()
plot(as, λ1)
title("Largest Lyapunv exponent for Roessler")
grid()

fig, ax1 = subplots()
ax1.plot(as, λ1, linewidth=2)
xlabel("control parameter a")
ylabel("λ₁")
ax2 = ax1.twinx()
Nplot = 100
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
