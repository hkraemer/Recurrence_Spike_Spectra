using InterSpikeSpectra
using LinearAlgebra
using Random
using Statistics
using StatsBase
using TimeseriesSurrogates
using DelayEmbeddings
using RecurrenceAnalysis
using Distributions

"""
    tau_spectrum(Y::Vector, ϵ::Real=0.05; kwargs...)  → ISS_spectrum

    Compute the inter spike spectrum of the τ-recurrence rate `ISS_spectrum`
    from a vector or trajectory `Y`, the recurrence threshold `ϵ` and a fixed
    recurrence rate.

    Keyword arguments:
    `ρ_thres = 0.99`: The agreement of the regenerated decomposed signal with the
                      true signal. This depends on the LASSO regularization parameter
                      `λ`. `λ` gets adjusted automatically with respect to `ρ_thres`.
    `tol = 1e-3`: Allowed tolerance between `ρ_thres` and `ρ`.
    `max_iter = 20`: Determines after how many tried Lambdas the algorithm stopps.
    `method::AbstractRegressionMethod = lasso()` : The method for sparse regression. 
                        Pick either lasso() or STLS() (sequential thresholded least squares)
    `regression_type::AbstractRegressionType = auto()` : Regression type. For probability 
                        input (like the τ-recurrence rate), a logit regression is needed 
                        'regression_type=logit()'. If this is not the case, select 
                        'regression_type=normal()'. By default it is automatically set with 
                        respect to the distribution of `s`.
"""
function tau_spectrum(Y::Union{Vector,Dataset}, ϵ::Real=0.05; fixedrate::Bool=true, ρ_thres::Real = 0.99,
        tol::Real=1e-3, max_iter::Integer=20, regression_type=InterSpikeSpectra.normal(), 
        method=InterSpikeSpectra.lasso())

    @assert 0.8 <= ρ_thres <= 1 "Optional input `ρ_thres` must be a value in the interval [0.8, 1]"
    @assert 1e-5 <= tol <= 1 "Optional input `tol` must be a value in the interval [1e-5, 1]"
    @assert 0 < maxλ <= 100 "Optional input `maxλ` must be an integer in the interval (0, 100]."

    # create RP and τ-recurrence rate from trajectory Y
    RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ϵ; fixedrate = fixedrate)
    τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
    τ_rr /= maximum(τ_rr)

    spectrum, _ = inter_spike_spectrum(τ_rr; ρ_thres, tol, max_iter, regression_type, method)
    return spectrum
end
 

"""
    compute_surrogate_τ_RR(NumTrials::Int, p::Real, N::Int) → τ_RR_surrogates

    Compute `NumTrials` surrogate τ-recurrence rate time series of length `N` under
    H₀ that the data is randomly distributed. That way the τ-recurrence rate time
    series can be obtained from a binomial distribution with parameter `p`, which
    corresponds to the fixed recurrence rate of the underlying RP.
"""
function compute_surrogate_τ_RR(NumTrials::Int, p::Real, N::Int)
    @assert NumTrials > 0 "Number of trials must be an positive integer value."
    @assert p > 0 "Provide a positive probability."
    @assert N > 0 "Length of sample spectra must be positive."
    if NumTrials > 1000
        println("You attempt to produce a sample with size larger than 1,000.")
    end
    if N > 10000
        println("You attempt to produce a samples each of which is longer than 10,000.")
    end

    NN = Int(N + floor(N*(1/3)))   # increase length, in order to avoid artifacts

    τ_RR_surrogates = ones(NumTrials,NN)
    # iterate over the different trials for surrogates time series
    for i = 1:NumTrials
        # draw from binomial distribution for each tau value
        for j = 2:NN
            pos_num = NN-j+1
            τ_RR_surrogates[i,j] = rand(Binomial(pos_num,p))/pos_num
        end
    end
    return τ_RR_surrogates[:, 1:N]
end

"""
    compute_surrogate_τ_RR_iAAFT(NumTrials::Int, x::Vector, ε::Real, τs::Vector) → τ_RR_surrogates

    Compute `NumTrials` surrogate τ-recurrence rate time series of the input time
    series `x`, which has length `M`, under the recurrence threshold `ε`. `τs` stores
    the delay values for the embedding of each iAAFT-surrogate, and this also determines
    the length `N` of the embedded trajectory and, thus, the length of each surrogate
    τ-recurrence rate time series. The output `τ_RR_surrogates` is a
    (`NumTrials`-by-`N`) Dataset.
"""
function compute_surrogate_τ_RR_iAAFT(NumTrials::Int, x::Vector{T}, ε::Real,
    τs::Vector{Int}; fixedrate = false) where {T}

    @assert NumTrials > 0 "Number of trials must be an positive integer value."
    @assert 1 > ε > 0 "Recurrence threshold must be ∈ [0,1]."
    if NumTrials > 1000
        println("You attempt to produce a sample with size larger than 1,000.")
    end
    Y_test = DelayEmbeddings.genembed(x, τs)
    N = length(Y_test)

    τ_RR_surrogates = ones(NumTrials,N)
    # iterate over the different trials for surrogates time series
    for i = 1:NumTrials
        s = surrogate(x, IAAFT())
        Y = DelayEmbeddings.genembed(s, τs)
        RP = RecurrenceAnalysis.RecurrenceMatrix(Y, ε; fixedrate = fixedrate)
        τ_rr = RecurrenceAnalysis.tau_recurrence(RP)
        τ_rr = τ_rr ./ maximum(τ_rr)
        τ_RR_surrogates[i,:] = τ_rr
    end
    return τ_RR_surrogates
end


"""
    compute_surrogate_spectra(τ_RR_surrogates::AbstractMatrix; kwargs...) → τ_surrogate_spectra

    Compute `τ-recurrence-spectra` for all surrogate τ-recurrence rate time series in
    `τ_RR_surrogates`.

    Keyword arguments:
    `ρ_thres = 0.99`: The agreement of the regenerated decomposed signal with the
                      true signal. This depends on the LASSO regularization parameter
                      `λ`. `λ` gets adjusted automatically with respect to `ρ_thres`.
    `tol = 1e-3`: Allowed tolerance between `ρ_thres` and `ρ`.
    `max_iter = 20`: Determines after how many tried Lambdas the algorithm stopps.
"""
function compute_surrogate_spectra(τ_RR_surrogates::Union{AbstractMatrix,AbstractVector}; ρ_thres::Real = 0.99,
        tol::Real=1e-3, max_iter::Integer=20, regression_type=InterSpikeSpectra.normal(), 
        method=InterSpikeSpectra.lasso())

    if typeof(τ_RR_surrogates)<:AbstractVector
        NumTrials = 1
        N = length(τ_RR_surrogates)
        τ_surrogate_spectra = zeros(Int(ceil(N/2)))
    else
        NumTrials, N = size(τ_RR_surrogates)
        τ_surrogate_spectra = zeros(NumTrials, Int(ceil(N/2)))
    end
    for i = 1:NumTrials
        if typeof(τ_RR_surrogates)<:AbstractVector
            τ_surrogate_spectra[:], _ = inter_spike_spectrum(τ_RR_surrogates;
                                        ρ_thres, tol, max_iter, regression_type, method, verbose=false)
        else
            τ_surrogate_spectra[i,:], _ = inter_spike_spectrum(τ_RR_surrogates[i,:];
            ρ_thres, tol, max_iter, regression_type, method, verbose=false)
        end
    end
    return τ_surrogate_spectra
end

"""
    compute_percentiles_of_surrogate_spectra(τ_surrogate_spectra::AbstractMatrix, α::Real = 0.05) → upper_bound, lower_bound

    For a significance-level `α` compute the `upper_bound` and `lower_bound` of the
    confidence interval obtained from the surrogate τ-recurrence spectra stored in
    `τ_surrogate_spectra`.
"""
function compute_percentiles_of_surrogate_spectra(τ_surrogate_spectra::Union{AbstractMatrix,AbstractVector}, α::Real = 0.05)
    @assert 1 > α > 0 "Significance level α must be ∈ [0, 1]"

    NumTrials, N = size(τ_surrogate_spectra)
    upper_bound = zeros(N)
    lower_bound = zeros(N)

    # adjust α for multiple comparison problem (Sidak correction)
    α_m = 1-(1-α)^(1/N)

    for i = 1:N
        upper_bound[i] = Statistics.quantile(τ_surrogate_spectra[:,i],1-α)
        lower_bound[i] = Statistics.quantile(τ_surrogate_spectra[:,i],α)
    end

    return upper_bound, lower_bound
end

"""
    logistic_map_time_series(N::Int, r::Real; u0 = rand()) → s

Compute a time series `s` of length `N` of the logistic map with control parameter
`r` and a initial condition `u0`.
"""
function logistic_map_time_series(N::Int, r::Real; u0 = rand())
    transients = 1000
    s = zeros(N+transients)
    s[1] = u0
    for i = 2:N+transients
        s[i] = (r * s[i-1]) - (r * s[i-1]^2)
    end
    return s[transients+1:end]
end


"""
Return the maxima of the given time series s and its indices
"""
function get_maxima(s::Vector{T}) where {T}
    maximas = T[]
    maximas_idx = Int[]
    N = length(s)
    flag = false
    first_point = 0
    for i = 2:N-1
        if s[i-1] < s[i] && s[i+1] < s[i]
            flag = false
            push!(maximas, s[i])
            push!(maximas_idx, i)
        end
        # handling constant values
        if flag
            if s[i+1] < s[first_point]
                flag = false
                push!(maximas, s[first_point])
                push!(maximas_idx, first_point)
            elseif s[i+1] > s[first_point]
                flag = false
            end
        end
        if s[i-1] < s[i] && s[i+1] == s[i]
            flag = true
            first_point = i
        end
    end
    # make sure there is no empty vector returned
    if isempty(maximas)
        maximas, maximas_idx = findmax(s)
    end
    return maximas, maximas_idx
end
