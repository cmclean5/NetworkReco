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

steps <- vector(length=points)
steps[1]    = 5
steps[2:5]  = 25
steps[5:11] = 250


#--- generate our ground truth Aij's for each point
Aij <- list()

#--- generate Bernoulli random graph with n=200 and <k> = 10
adj <- sample_gnm(n=n,m=m,directed=F,loops=F)

for( p in 1:points ){
    Aij[[p]]      <- adj
    names(Aij)[p] <- sprintf("Aij_%d",p)
}


# test
if ( 0 ){
    my.params <- list()
    my.params[[1]] = m/choose(n,2)
    names(my.params)[1] = "rho"
    my.params[[2]] = 0
    names(my.params)[2] = "alpha"
    my.params[[3]] = 0
    names(my.params)[3] = "beta"

    my.const <- list()
    my.const[[1]] = c(0.01, 0.01)
    names(my.const)[1] = "rho"
    my.const[[2]]  = c(0.1, 0.1)
    names(my.const)[2] = "alpha"
    my.const[[3]] = c(0.1, 0.1)
    names(my.const)[3] = "beta"
}

recall     = matrix(0,ncol=(1+length(Nmeas)), nrow=points )
recall[,1] = errRates

conv       = matrix(1,ncol=(1+length(Nmeas)), nrow=points )
conv[,1]   = errRates

diff       = matrix(0,ncol=(1+length(Nmeas)), nrow=points )
diff[,1]   = errRates


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
    restarts=3#10
    tol=1e-7#1e-5
    conv.params = TRUE
    for( p in 1:points ){

        
        RES[[p]] = run.em( Adj=Aij[[p]], meas=Nmeas[N], obs=Ecount[[p]],
                          max.steps=steps[p], restarts=restarts,
                          store.delta.N=ceiling(steps[p]/2),
                          conv.params=conv.params, tol=tol )
        names(RES)[p] = sprintf("N_%d_errorRate_%.1f",Nmeas[N],errRates[p])

        cat("> RES for errorRate : ", errRates[p], " converage: ", RES[[p]]@converage,".\n")
        
        val = RES[[p]]@Qij[get.edgelist(Aij[[p]])]
        recall[p,(N+1)] = sum(val > 0.5)/m
        conv[p,(N+1)]   = RES[[p]]@converage
        diff[p,(N+1)]   = RES[[p]]@delta.diff
        
        write.table(recall,"recall.csv",sep="\t", row.names=F, col.names=T, quote=F)
        write.table(conv,"converage.csv",sep="\t", row.names=F, col.names=T, quote=F)
        write.table(diff,"deltadiff.csv",sep="\t", row.names=F, col.names=T, quote=F)
        
    }

}
    


