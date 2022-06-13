# Here we compute the optimal embedding parameters for the Laskar eccentricity data and for the 
# CENOGRID data in order to compute RPs, τ-RRs and the according inter spike spectra. We use 
# PECUZAL algorithm and the corresponding MCDTS-routine.

using DelayEmbeddings
using TreeEmbedding  # this will be merged into DelayEmbeddings soon
using DelimitedFiles
using Random

## Eccentricity
ecc = vec(readdlm("./Applications/CENOGRID/data/Laskar/ecc_Laskar.txt"))

# PECUZAL
w = estimate_delay(ecc, "mi_min", 0:200)
L_threshold = 0.005
KNN = 5
τs = 0:150 # sampling time is 5 kyrs, so this covers 750 kyrs
Y, τ_pec_ecc, _, Ls , _ = pecuzal_embedding(ecc; w, L_threshold, KNN, τs, econ=true)

# MCDTS -> PECUZAL
runs = 50
tws = 2:2:τs[end]
pecuzal = TreeEmbedding.MCDTSOptimGoal(TreeEmbedding.L_statistic(-L_threshold,KNN,tws), TreeEmbedding.Continuity_function())
Random.seed!(1)
tree = mcdts_embedding(ecc, pecuzal, w, τs, runs; verbose=true)
best_node = TreeEmbedding.best_embedding(tree)
τ_mcdts_ecc_pec = best_node.τs
# MCDTS -> FNN + Continuity_function
optimgoal = TreeEmbedding.MCDTSOptimGoal(TreeEmbedding.FNN_statistic(), TreeEmbedding.Continuity_function())
Random.seed!(1)
tree = mcdts_embedding(ecc, optimgoal, w, τs, runs; verbose=true)
best_node = TreeEmbedding.best_embedding(tree)
τ_mcdts_ecc_fnn = best_node.τs

#out
# τ_pec_ecc = [0, 8, 4, 21, 15, 33, 27, 128, 74, 135, 64]
# τ_mcdts_ecc_pec = [0, 8, 4, 21, 15, 33, 27, 128, 74, 135, 64]
# τ_mcdts_ecc_fnn = [0  8  45  4  82  18  34  13  10]
# τ_mcdts_ecc_fnn = [0, 8, 62, 115, 4, 3]
# τ_mcdts_ecc_fnn = [0, 91, 130, 8, 47, 4, 2]

## CENOGRID filtered
O18 = vec(readdlm("./Applications/CENOGRID/data/O18_filtered_time_reverse.txt"))
C13 = vec(readdlm("./Applications/CENOGRID/data/C13_filtered_time_reverse.txt"))

w_O18 = estimate_delay(O18, "mi_min", 0:200)
w_C13 = estimate_delay(C13, "mi_min", 0:200)

# MCDTS -> FNN + Continuity_function
runs = 50
optimgoal = TreeEmbedding.MCDTSOptimGoal(TreeEmbedding.FNN_statistic(0.01), TreeEmbedding.Continuity_function())
Random.seed!(1)
tree_O18 = mcdts_embedding(O18, optimgoal, w_O18, τs, runs; verbose=true)
best_node_O18 = TreeEmbedding.best_embedding(tree_O18)
τ_mcdts_O18_fnn = best_node_O18.τs
Random.seed!(1)
tree_C13 = mcdts_embedding(C13, optimgoal, w_C13, τs, runs; verbose=true)
best_node_C13 = TreeEmbedding.best_embedding(tree_C13)
τ_mcdts_C13_fnn = best_node_C13.τs

#out
# τ_mcdts_O18_fnn = [0  50  104  35  43  149  39  47  45], threshold 0
# τ_mcdts_O18_fnn = [0  50  142  105  130  136  133], threshold 0.01
# τ_mcdts_C13_fnn = [0  136  32  123  130  126  133  128] , threshold 0
# τ_mcdts_C13_fnn = [0  136  32  101  114  53  88  42  125  37], threshold 0.01


## CENOGRID raw
data = vec(readdlm("./Applications/CENOGRID/data/detrended.txt"))
len = Int(length(data)/3)
C13 = reverse(data[len+1:2*len])
O18 = reverse(data[2*len+1:end])

τs = 0:150 # sampling time is 5 kyrs, so this covers 750 kyrs
w_O18 = estimate_delay(O18, "mi_min", τs)
w_C13 = estimate_delay(C13, "mi_min", τs)

# PECUZAL
L_threshold = 0
KNN = 3
Y, τ_pec_O18, _, Ls , _ = pecuzal_embedding(O18; w=w_O18, L_threshold, KNN, τs, econ=true)
Y, τ_pec_C13, _, Ls , _ = pecuzal_embedding(C13; w=w_C13, L_threshold, KNN, τs, econ=true)

# MCDTS -> FNN + Continuity_function
runs = 50
optimgoal = TreeEmbedding.MCDTSOptimGoal(TreeEmbedding.FNN_statistic(), TreeEmbedding.Continuity_function())
Random.seed!(1)
tree_O18 = mcdts_embedding(O18, optimgoal, w_O18, τs, runs; verbose=true)
best_node_O18 = TreeEmbedding.best_embedding(tree_O18)
τ_mcdts_O18_raw_fnn = best_node_O18.τs
Random.seed!(1)
tree_C13 = mcdts_embedding(C13, optimgoal, w_C13, τs, runs; verbose=true)
best_node_C13 = TreeEmbedding.best_embedding(tree_C13)
τ_mcdts_C13_raw_fnn = best_node_C13.τs

#out
# τ_mcdts_O18_fnn_raw = [0  85  43  78  72  65  21  75  77  76], threshold 0
# τ_mcdts_O18_fnn_raw = [], threshold 0.01
# τ_mcdts_C13_fnn_raw = [0  113  53  27  68  60  56  58  57  59] , threshold 0
# τ_mcdts_C13_fnn_raw = [], threshold 0.01

