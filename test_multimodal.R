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
Nmeas    <- c(4,4,4)

errRates <- c(0.1,0.5,0.9)

#--- No: of modes
modes    <- length(errRates)

#--- generate our ground truth Aij's for each point
Aij <- list()

#--- generate Bernoulli random graph with n=200 and <k> = 10
adj <- sample_gnm(n=n,m=m,directed=F,loops=F)

for( p in 1:modes ){
    Aij[[p]]      <- adj
    names(Aij)[p] <- sprintf("Aij_%d",p)
}


#--- generate measurements for each Aij
Ecount <- list()
Ncount <- list()

for( k in 1:modes ){

    N = sample(Nmeas,1)

    Ncount[[k]]  <- N    
    Ecount[[k]]  <- generateMeasurements(GG=Aij[[k]],GGname=names(Aij)[k],
                                         Nrand=1, errorRate=errRates[k],
                                         Nmeas=N)
    names(Ecount)[k] <- sprintf("Ecount_%d",k)

}

cat('done!\n')


# test
if ( 0 ){
    my.params <- list()
    my.params[[1]] = m/choose(n,2)
    names(my.params)[1] = "rho"
    my.params[[2]] = matrix(0,ncol=modes, nrow=1)
    names(my.params)[2] = "alpha"
    my.params[[3]] = matrix(0,ncol=modes, nrow=1)
    names(my.params)[3] = "beta"

    my.const <- list()
    my.const[[1]] = matrix(0.01, ncol=1, nrow=2)
    names(my.const)[1] = "rho"
    my.const[[2]]  = matrix(0.1, ncol=modes, nrow=2)
    names(my.const)[2] = "alpha"
    my.const[[3]] = matrix(0.1, ncol=modes, nrow=2)
    names(my.const)[3] = "beta"
}

   
RES <- list()
steps=10
restarts=10
epsilon=1e-100
tol=1e-5

RES[[1]] = run.em(Adj=Aij[[1]], meas=Ncount, obs=Ecount,
                  max.steps=steps, restarts=restarts, epsilon=epsilon,
                  store.delta.N=ceiling(steps/2), conv.params=TRUE, tol=tol)

val = RES[[1]]@Qij[get.edgelist(Aij[[1]])]
recall = sum(val > 0.5)/m
conv   = RES[[1]]@converage
diff   = RES[[1]]@delta.diff
        


