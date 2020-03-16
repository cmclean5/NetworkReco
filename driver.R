# References:
# [1] M. E. J. Newman, Network sturcture from rich but noisy data, nature physics, 2018.
# [2] Supp. Material, https://doi.org/10.1038/s41567-018-0076-1
# [3] M.E. J. Newman, Network structure from rich but noisy data, arXiv:1703.07376v2, 2018.
# [4] example EM algo structure: http://www.di.fc.ul.pt/~jpn/r/EM/EM.html


# reset 
rm(list=ls())

#--- imports/headers
require(igraph)
library(methods)

# load EM functions
source('EM.R')
# load Random_Measurement functions
source('Generate_Rndm_Measurements.R')

# Generate set of random measurements from test graphs 
cat('Generate oberserved and measured datasets...\n')

#--- number of nodes
n = 200

#--- number of edges
m = 1000

#--- No: of measurements
Nmeas = c(4,8,16)

N     = 1

errRates <- seq(0,m,100)/m

#--- number of points
points <- length(errRates)


#--- generate our ground truth Aij's for each point
Aij <- list()

#--- generate Bernoulli random graph with n=200 and <k> = 10
adj <- sample_gnm(n=n,m=m,directed=F,loops=F)

for( p in 1:points ){
    Aij[[p]]      <- adj
    names(Aij)[p] <- sprintf("Aij_%d",p)
}


#--- generate measurements for each Aij
#Ecount <- list()
#
#for( p in 1:points ){
#
#    Ecount[[p]]  <- generateMeasurements(GG=Aij[[p]],GGname=names(Aij)[p],
#                                         Nrand=1, errorRate=errRates[p],
#                                         Nmeas=Nmeas[N])
#    names(Ecount)[p] <- sprintf("Ecount_%d",p)
#
#}
#
#cat('done!\n')

#epsilon=10^-9

# test
if ( 1 ){
    my.params <- list()
    my.params[[1]] = m/choose(n,2)
    names(my.params)[1] = "rho"
    my.params[[2]] = 0
    names(my.params)[2] = "alpha"
    my.params[[3]] = 0
    names(my.params)[3] = "beta"
}


recall     = matrix(0,ncol=(1+length(Nmeas)), nrow=points )
recall[,1] = errRates


for(N in 1:length(Nmeas) ){

    #--- generate measurements for each Aij
    cat("Generate measurments for N = ", Nmeas[N], ".\n")
    Ecount <- list()

    for( p in 1:points ){

        Ecount[[p]]  <- generateMeasurements(GG=Aij[[p]],GGname=names(Aij)[p],
                                             Nrand=1, errorRate=errRates[p],
                                             Nmeas=Nmeas[N])
        names(Ecount)[p] <- sprintf("Ecount_%d",p)

    }

    cat('done!\n')

    
    RES <- list()
    for( p in 1:points ){

        #RES[[p]] = run.em(Adj=Aij[[p]], Nmeas=Nmeas, obs=Ecount[[p]], study=N, errRATE=errRates[p], initPARAMS=NULL)
        RES[[p]] = run.em(Adj=Aij[[p]], Nmeas=Nmeas, obs=Ecount[[p]], study=N, errRATE=errRates[p], initPARAMS=my.params, fixPARAMS=c(1,0,0), max.steps=100)
        names(RES)[p] = sprintf("N_%d_errorRate_%.1f",Nmeas[N],errRates[p])
    
        val = RES[[p]]@Qij[get.edgelist(Aij[[p]])]
        recall[p,(N+1)] = sum(val > 0.5)/m
   
    }

}
    
write.table(recall,"recallv3.csv",sep="\t", row.names=F, col.names=T, quote=F)

