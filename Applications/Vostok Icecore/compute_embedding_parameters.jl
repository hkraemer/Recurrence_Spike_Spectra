# Here we compute the optimal embedding parameters for the Vostok data

using DelayEmbeddings
using TreeEmbedding  # this will be merged into DelayEmbeddings soon
using DelimitedFiles
using Random

## Vostok temperature
temp = vec(readdlm("./Applications/Vostok Icecore/data/vostok_interpolated_temperature.csv"))
temp = temp .+ 0.000000001 .* randn()

# PECUZAL
w = estimate_delay(temp, "mi_min", 0:150)
L_threshold = 0.005
KNN = 5
τs = 0:150 # sampling time is 81 yrs, so this covers 12150 years
Y, τ_pec_temp, _, Ls , _ = pecuzal_embedding(temp; w, L_threshold, KNN, τs, econ=true)

# MCDTS -> FNN + Continuity_function
runs = 50
fnn_threshold = 0.01 # 1% false neighbors
optimgoal = TreeEmbedding.MCDTSOptimGoal(TreeEmbedding.FNN_statistic(fnn_threshold), TreeEmbedding.Continuity_function())
Random.seed!(2)
tree = mcdts_embedding(temp, optimgoal, w, τs, runs; verbose=true)
best_node = TreeEmbedding.best_embedding(tree)
τ_mcdts_temp_fnn = best_node.τs

#out
# τ_pec_temp = []
# τ_mcdts_temp_fnn = [0, 140, 46, 84, 112, 66, 57]
