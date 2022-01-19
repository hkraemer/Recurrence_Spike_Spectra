using DelimitedFiles
using PyPlot
using Statistics
pygui(true)

include("../../Methods/handy_functions.jl")
# Parameter-values for Logistic map
rs = 3.4:0.001:4

# delay for embedding
τ = 1
# recurrence threshold
ε = 0.05
# length of time series
N = 201
# threshold for ISS
thres = 0.99

# spectra and upper CI of random surrogates
upper = vec(readdlm("./Applications/Logistic Map Example/Results/results_Logistic_N_$(N)_thres_$(thres)_upper.csv"))
# upper CI of iAAFT surrogates
upper2 = vec(readdlm("./Applications/Logistic Map Example/Results/results_Logistic_N_$(N)_thres_$(thres)_upper2.csv"))
# number of significant peaks vs random H0
nsp = vec(readdlm("./Applications/Logistic Map Example/Results/results_Logistic_N_$(N)_thres_$(thres)_nsp.csv"))
# number of significant peaks vs iAAFT H0
nsp2 = vec(readdlm("./Applications/Logistic Map Example/Results/results_Logistic_N_$(N)_thres_$(thres)_nsp2.csv"))

# compute Lyapunov exponent
lyap = zeros(length(rs))
lyap_pos = zeros(length(rs))
for i = 1:length(rs)
    s = logistic_map_time_series(N, rs[i])
    lyap[i] = mean(log.(abs.(rs[i] .- 2*rs[i].*s)))
    if lyap[i]>0
        lyap_pos[i]=lyap[i]
    end
end

begin
    ## Corr coeff
    coeff1 = round(Statistics.cor(lyap_pos,nsp); digits = 2)
    coeff2 = round(Statistics.cor(lyap_pos,nsp2); digits = 2)

    figure(figsize=(15,10))

    ax1 = subplot(311)
    ax1.plot(rs, lyap, linewidth=2)
    hlines(0, rs[1], rs[end], linestyle="dashed")
    #xlabel("control parameter r")
    ylabel("λ₁")
    grid()
    xlim([rs[1], rs[end]])
    ax2 = ax1.twinx()
    grid()
    Nplot = 100
    for r in rs
      s = logistic_map_time_series(Nplot, r)
      ax2.scatter(r*ones(Nplot), s, s=.1, c="gray", alpha=0.2)
    end
    PyPlot.yticks([])
    title("Lyapunov exponent of the Logistic map")
    grid()


    ax3 = subplot(312)
    ax3.plot(rs, nsp, linewidth = 2)
    textstr = "ρ Pearson = $coeff1"
    props = Dict("boxstyle"=>"round", "facecolor"=>"wheat", "alpha"=>0.5)
    ax3.text(0.03, 0.93, textstr, transform=ax3.transAxes, fontsize=10,
            verticalalignment="top", bbox=props)
    ylabel("no. of significant peaks")
    title("no. of significant peaks in the inter spike spectrum of τ-RR (H0 random)")
    xlim([rs[1], rs[end]])
    ylim([-2, 60])
    grid()

    ax4 = subplot(313)
    ax4.plot(rs, nsp2, linewidth = 2)
    textstr = "ρ Pearson = $coeff2"
    props = Dict("boxstyle"=>"round", "facecolor"=>"wheat", "alpha"=>0.5)
    ax4.text(0.03, 0.93, textstr, transform=ax4.transAxes, fontsize=10,
            verticalalignment="top", bbox=props)
    xlabel("control parameter r")
    ylabel("no. of significant peaks")
    title("no. of significant peaks in the spike powerspectrum of τ-RR (H0 iAAFT)")
    xlim([rs[1], rs[end]])
    ylim([-2, 60])
    grid()
    subplots_adjust(hspace=.4)
end
