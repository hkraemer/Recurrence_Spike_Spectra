# Here we compute the optimal embedding parameters for the Laskar eccentricity data and for the 
# CENOGRID data in order to compute RPs, τ-RRs and the according inter spike spectra. We use 
# PECUZAL algorithm and the corresponding MCDTS-routine.

using DelayEmbeddings
using TreeEmbedding  # this will be merged into DelayEmbeddings soon
using DelimitedFiles
using Random

## Eccentricity
ecc = vec(readdlm("./Applications/Laskar/data/Laskar/ecc_Laskar.txt"))

# PECUZAL
w = estimate_delay(ecc, "mi_min", 0:200)
L_threshold = 0.005
KNN = 5
τs = 0:150 # sampling time is 5 kyrs, so this covers 750 kyrs
Y, τ_pec_ecc, _, Ls , _ = pecuzal_embedding(ecc; w, L_threshold, KNN, τs, econ=true)

# MCDTS -> FNN + Continuity_function
runs = 50
fnn_threshold = 0.01 # 1% false neighbors
optimgoal = TreeEmbedding.MCDTSOptimGoal(TreeEmbedding.FNN_statistic(fnn_threshold), TreeEmbedding.Continuity_function())
Random.seed!(2)
tree = mcdts_embedding(ecc, optimgoal, w, τs, runs; verbose=true)
best_node = TreeEmbedding.best_embedding(tree)
τ_mcdts_ecc_fnn = best_node.τs

#out
# τ_pec_ecc = [0, 8, 4, 21, 15, 33, 27, 128, 74, 135, 64]
# τ_mcdts_ecc_fnn = [0  8  101  62  114  36  18  4]
