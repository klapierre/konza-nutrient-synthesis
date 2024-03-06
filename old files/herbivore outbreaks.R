herbivores <- read.csv('konza_herbivore_site mean.csv')

hist(herbivores$grasshopper_abund, breaks=10, xlab='Grasshopper Abundance')
#export at 600x600

hist(herbivores$small_mammal_abund, breaks=10, xlab='Small Mammal Abundance')
#export at 600x600
