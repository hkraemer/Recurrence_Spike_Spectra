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
thresholds = [0.9, 0.95, 0.99]

nsp = zeros(length(thresholds), length(rs))
nsp_iafft = zeros(length(thresholds), length(rs))


begin
    # compute Lyapunov exponent
    lyaps = zeros(length(rs))
    lyap_pos = zeros(length(rs))
    for i = 1:length(rs)
        s = logistic_map_time_series(N, rs[i])
        lyaps[i] = mean(log.(abs.(rs[i] .- 2*rs[i].*s)))
        if lyaps[i]>0
            lyap_pos[i]=lyaps[i]
        end
    end

    ## Corr coeff preallocation
    coeff1 = zeros(length(thresholds))
    coeff2 = zeros(length(thresholds))

    for i = 1:lenght(thresholds)

        # number of significant peaks vs random H0 (normal LASSO)
        nsp[i,:] = vec(readdlm("./Applications/Logistic Map Example/Results/results_Logistic_N_$(N)_thres_$(thresholds[i])_nsp.csv"))
        # number of significant peaks vs iAAFT H0 (normal LASSO)
        nsp_iafft[i,:] = vec(readdlm("./Applications/Logistic Map Example/Results/results_Logistic_N_$(N)_thres_$(thresholds[i])_nsp_iafft.csv"))


        coeff1[i] = round(Statistics.cor(lyap_pos,nsp[i,:]); digits = 2)
        coeff2[i] = round(Statistics.cor(lyap_pos,nsp_iafft[i,:]); digits = 2)

        ## results for normal LASSO regression
        figure(figsize=(15,10))

        ax1 = subplot(311)
        ax1.plot(rs, lyaps, linewidth=2)
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
        ax3.plot(rs, nsp[i,:], linewidth = 2)
        textstr = "ρ Pearson = $(coeff1[i])"
        props = Dict("boxstyle"=>"round", "facecolor"=>"wheat", "alpha"=>0.5)
        ax3.text(0.03, 0.93, textstr, transform=ax3.transAxes, fontsize=10,
                verticalalignment="top", bbox=props)
        ylabel("no. of significant peaks")
        title("Inter spike spectrum of τ-RR (H0 random, normal LASSO, thres=$(thresholds[i]))")
        xlim([rs[1], rs[end]])
        ylim([-2, 60])
        grid()

        ax4 = subplot(313)
        ax4.plot(rs, nsp_iafft[i,:], linewidth = 2)
        textstr = "ρ Pearson = $(coeff2[i])"
        props = Dict("boxstyle"=>"round", "facecolor"=>"wheat", "alpha"=>0.5)
        ax4.text(0.03, 0.93, textstr, transform=ax4.transAxes, fontsize=10,
                verticalalignment="top", bbox=props)
        xlabel("control parameter r")
        ylabel("no. of significant peaks")
        title("Inter spike spectrum of τ-RR  (H0 iAAFT, normal LASSO, thres=$(thresholds[i]))")
        xlim([rs[1], rs[end]])
        ylim([-2, 60])
        grid()
        subplots_adjust(hspace=.4)

    end
end
