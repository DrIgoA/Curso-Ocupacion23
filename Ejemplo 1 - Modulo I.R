############################################################
###########          Ejemplo Muestreo            ###########
############################################################
########              Basado en                     ########
########      Kery, M., & Royle, J. A. (2016).      ########
########       Applied Hierarchical Modeling        ######## 
########              in Ecology                    ######## 
########                Volume 1                    ########
############################################################

set.seed(82)
tmp <-sim.fn(quad.size=10, cell.size=1, intensity =1)
install.packages("TMB",type = "source")
