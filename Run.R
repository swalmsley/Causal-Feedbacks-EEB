

library(targets)

###############
# Run analysis
tar_make() 
###############

# Visualize analysis
tar_visnetwork() # visualizes analysis pipeline


# Examine object from analysis
#tar_read(printsByPopulation, branches=1)

# Diagnostics
View(tar_meta()) # useful tool for diagnostics
View(tar_meta(targets_only = TRUE)) # simplified


# Notes 

# Feedbacks seen as barrier “Second, when causal relationships between the microbiome and behaviour do exist, their directionality may be hard to infer as effects in both directions, and feedback, may be expected” -- Identifying Microbiome-Mediated Behaviour in Wild Vertebrates















