# Here we compute the optimal embedding parameters for the Vostok data

using DelayEmbeddings
using TreeEmbedding  # this will be merged into DelayEmbeddings soon
using DelimitedFiles
using Random

## Downsampled frequency data
f_ce = vec(readdlm("./Applications/Power Grid/CE and GB frequencies/data/frequencies_ce_downsampled.csv"))
f_gb = vec(readdlm("./Applications/Power Grid/CE and GB frequencies/data/frequencies_gb_downsampled.csv"))

# Cao Embedding
w = estimate_delay(f_ce[1:14400], "mi_min", 0:150)
D_ce, τ_ce, E_ce = optimal_traditional_de(f_ce[1:14400], "ifnn", "mi_min"; τs, w)

w = estimate_delay(f_gb[1:14400], "mi_min", 0:150)
D_gb, τ_gb, E_gb = optimal_traditional_de(f_gb[1:14400], "ifnn", "mi_min"; τs, w)

# out Cao
# D_ce = 7, τ_ce = 19 
# D_gb = 7, τ_gb = 17

# changing embedding params for different signal lengths