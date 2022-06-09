# we compute the optimal embedding parameters for the Laskar eccentricity data and for the 
# CENOGRID data in order to compute RPs, τ-RRs and the according inter spike spectra

using DelayEmbeddings
using TreeEmbedding  # this will be merged into DelayEmbeddings soon
using DelimitedFiles

## Eccentricity
ecc = vec(readdlm("./Applications/CENOGRID/data/Laskar/ecc_Laskar.txt"))

# PECUZAL
w = estimate_delay(ecc, "mi_min", 0:200)
L_threshold = 0.005
KNN = 5
τs = 0:150
Y, τ_pec, _, Ls , _ = pecuzal_embedding(ecc; w, L_threshold, KNN, τs, econ=true)

# MCDTS
runs = 50
tws = 2:2:τs[end]
pecuzal = TreeEmbedding.MCDTSOptimGoal(TreeEmbedding.L_statistic(-L_threshold,KNN,tws), TreeEmbedding.Continuity_function())
tree = mcdts_embedding(ecc, pecuzal, w, τs, runs; verbose=true)
best_node = TreeEmbedding.best_embedding(tree)
τ_mcdts = best_node.τs

#out
# τ_pec = [0, 8, 4, 21, 15, 33, 27, 128, 74, 135, 64]
# τ_mcdts = []

## O18

data_O18 = 