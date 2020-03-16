# References:
# M. E. J. Newman, Network sturcture from rich but noisy data, nature physics, 2018.
# example EM algo structure: http://www.di.fc.ul.pt/~jpn/r/EM/EM.html

#library(methods)

#--- return type for em algo
setClass(Class="runEM",representation(
                           Eij="list",
                           Nij="matrix",
                           params="list",
                           Qij="matrix",
                           steps="numeric",
                           ll="numeric") )			   


#--- return type for em algo
setClass(Class="EM",representation(
                        params="list",
                        Qij="matrix",
                        steps="numeric",
                        ll="numeric") )			   


print_params <- function( theta, modes ){

    rho   = theta[[which(names(theta)=="rho")]]
    alpha = theta[[which(names(theta)=="alpha")]]
    beta  = theta[[which(names(theta)=="beta")]]

    
    header = sprintf("\t modes %d \t",seq(1,modes,1))
    cat(header, "\n")
    cat("rho: ")
    for( i in 1:length(rho) ){
        cat("\t", rho[i])
    }
    cat("\n")

    cat("alpha: ")
    for( i in 1:length(alpha) ){
        cat("\t", alpha[i])
    }
    cat("\n")

    cat("beta: ")
    for( i in 1:length(beta) ){
        cat("\t", beta[i])
    }
    cat("\n")

}

#--- init model parameters to random values
init_params <- function( theta, fix ){

    N <- length(names(theta))
    for( i in 1:N ){
        if( !fix[i] ){
            for(m in 1:length(theta[[i]]) ){
                theta[[i]][1,m] <- runif(1)
            }
        }
    }
           
    return( theta )
}


init_qij <- function( QQ ){

    NROW = nrow(QQ)
    NCOL = ncol(QQ)
    
    return(matrix(runif(NROW*NCOL),NROW,NCOL))
    
}


Qij.param <- function( meas, obs, theta, modes ){

    Nij   = meas
    Eij   = obs
    
    rho   = theta[[which(names(theta)=="rho")]]
    alpha = theta[[which(names(theta)=="alpha")]]
    beta  = theta[[which(names(theta)=="beta")]]

    num.m = matrix(0, ncol=modes, nrow=1)
    dem.m = matrix(0, ncol=modes, nrow=1)

    for( k in 1:modes ){
        num.m[1,k] = rho * alpha[1,k]^(Eij[1,k]) * (1-alpha[1,k])^(Nij[1,k]-Eij[1,k])
    }

    for( k in 1:modes ){
        dem.m[1,k] = (1-rho) * beta[1,k]^(Eij[1,k]) * (1-beta[1,k])^(Nij[1,k]-Eij[1,k])
    }

    num = prod(num.m)
    dem = num + prod(dem.m)    
    
    Qij = 0

    if( dem != 0 ){ Qij = num / dem }
    
    return(Qij)
    
}

# The e-step calculates the prosterior distribution over networks, given 
# our model parameters, theta, from the m-step.
e.step <- function( n, QQ, meas, obs, theta){

    # the number of multimodes in study
    modes  <- length(obs) 

    nij.m = matrix(0, ncol=modes, nrow=1)
    eij.m = matrix(0, ncol=modes, nrow=1)
    
    for( i in 1:n ){
        for( j in 1:n ){

            if( i < j ){

                nij.m[1,] = 0
                eij.m[1,] = 0
                
                for( k in 1:modes ){                
                    nij.m[1,k] = meas[1,k]
                    eij.m[1,k] = obs[[k]][i,j]
                }
                
                QQ[i,j] = Qij.param( meas=nij.m, obs=eij.m, theta=theta, modes=modes ) 
            }

        }
    }

    return(QQ)
    
}


alpha.m <- function( qij, nij, eij ){

    return( c( eij*qij, nij*qij ) )
    
}


beta.m <- function( qij, nij, eij ){

    return( c( eij*(1-qij), nij*(1-qij) ) )
    
}


# At least obs must be a list of size m.
m.step <- function( n, QQ, meas, obs, theta ){

    # the number of multimodes in study
    modes  <- length(obs) 

    h_alpha = matrix(0,ncol=modes, nrow=1)
    h_beta  = matrix(0,ncol=modes, nrow=1)
    h_rho   = 0

    num_alpha = matrix(0,ncol=modes, nrow=1)
    num_beta  = matrix(0,ncol=modes, nrow=1)
    num_rho   = 0
    
    dem_alpha = matrix(0,ncol=modes, nrow=1)
    dem_beta  = matrix(0,ncol=modes, nrow=1)
    
    for( i in 1:n ){
        for( j in 1:n ){

            if( i < j ){

                qij = QQ[i,j]                

                num_rho = num_rho + qij
                
                for( k in 1:modes ){

                    nij.m = meas[1,k]#[i,j] #for the moment 
                    eij.m = obs[[k]][i,j]

                    temp.m = alpha.m(qij, nij.m, eij.m)
                    num_alpha[1,k] = num_alpha[1,k] + temp.m[1]
                    dem_alpha[1,k] = dem_alpha[1,k] + temp.m[2]

                    temp.m = beta.m(qij, nij.m, eij.m)
                    num_beta[1,k] = num_beta[1,k] + temp.m[1]
                    dem_beta[1,k] = dem_beta[1,k] + temp.m[2]
                    
                }

            }
            
        }
    }

    for( k in 1:modes ){
    
        if( dem_alpha[1,k] != 0 ){
            h_alpha[1,k] = num_alpha[1,k] / dem_alpha[1,k]
        }

        if( dem_beta[1,k] != 0 ){
            h_beta[1,k]  = num_beta[1,k] / dem_beta[1,k]
        }

    }
    
    dem_rho = choose( n, 2 ) 
    if( dem_rho != 0 ){
        h_rho   = num_rho / dem_rho
    }

    theta[[which(names(theta)=="rho")]]   = h_rho
    theta[[which(names(theta)=="alpha")]] = h_alpha
    theta[[which(names(theta)=="beta")]]  = h_beta
   
    #--- return result of running m.step
    return(theta)
    
}

adj.log <- function( arg.l=NULL, epsilon=NULL ){

    if( !is.null( arg.l ) ){

        if( is.null(epsilon) ){
            epsilon = .Machine$double.eps
        }

        if( arg.l == 0 ){ return( log( epsilon ) ) }
        else            { return( log( arg.l ) )   }
        
    }
    
}

# return the log-likelihood of the posterior, i.e. P(A,theta | data)
# See ref [2] eqn S4
posterior <- function( gg, QQ, meas, obs, theta, epsilon=NULL){

    n   <- length(V(gg))
    Aij <- get.adjacency(gg)

    modes <- length(obs)
    
    Nij   = meas
    Eij   = obs
    
    rho   = theta[[which(names(theta)=="rho")]]
    alpha = theta[[which(names(theta)=="alpha")]]
    beta  = theta[[which(names(theta)=="beta")]]
    
    #log likelihood
    ll   <- 0

    if( is.null(epsilon) ){
        epsilon = .Machine$double.eps
    }
    
     for( i in 1:n ){
        for( j in 1:n ){

            if( i < j ){
                                
                #qij = QQ[i,j]
                aij = Aij[i,j]               

                ll = ll + aij * adj.log(    rho,  epsilon=epsilon ) +
                      (1-aij) * adj.log( (1-rho), epsilon=epsilon )
                
                for( k in 1:modes ){

                    nij.m = Nij[1,k]
                    eij.m = Eij[[k]][i,j]
                    
                    alpha.m = alpha[1,k] 
                    beta.m  = beta[1,k]

                   ll = ll + aij * eij.m  * adj.log(    alpha.m,  epsilon=epsilon ) +
                    aij * (nij.m - eij.m) * adj.log( (1-alpha.m), epsilon=epsilon ) +
                         (1-aij) * eij.m  * adj.log(     beta.m,  epsilon=epsilon ) +
                (1-aij) * (nij.m - eij.m) * adj.log(  (1-beta.m), epsilon=epsilon )

                }                    
            }

        }
     }


    return(ll)
    
}



em <- function( Adj, Nij, Eij, Qij, params, modes, tol, max.steps, initPARAMS, fixPARAMS, epsilon ){

    #--- No: of nodes in graph    
    n <- length(V(Adj))
    
    #--- No: of edges in graph                  
    m <- length(E(Adj))

    modes <- modes
    
    #--- Init-step: assign random values to the model's parameters
    if( is.null(initPARAMS) ){
        Qij    <- init_qij( Qij )
        params <- m.step( n=n, QQ=Qij, meas=Nij, obs=Eij, theta=params )
    } else {
        params <- init_params(initPARAMS, fixPARAMS) 
        Qij    <- e.step( n=n, QQ=Qij, meas=Nij, obs=Eij, theta=params )
    }
    
    cat("init params: \n")
    print_params( params, modes )

    #--- em alg. steps    
    steps = maxdelta = 1

    #--- em test flags  
    emTEST1 = 0
    emTEST2 = 0
    emTEST3 = 0    

    diff = vector(length=2)

    #run the em algo.
    while( maxdelta > tol ) {
    
        #--- record old ll for run
        old.ll     <- posterior( gg=Adj, QQ=Qij, meas=Nij, obs=Eij, theta=params, epsilon=epsilon )

        #--- maximazation step
        old.params <- params
        params     <- m.step( n=n, QQ=Qij, meas=Nij, obs=Eij, theta=params )

        
        #--- expectation step
        Qij        <- e.step( n=n, QQ=Qij, meas=Nij, obs=Eij, theta=params)
        
                
        #--- record new ll for run
        new.ll     <- posterior( gg=Adj, QQ=Qij, meas=Nij, obs=Eij, theta=params, epsilon=epsilon )
        
        #--- test convergence of parameters        
        maxdelta = abs( abs(old.ll) - abs(new.ll) )
        cat(" LL_old ", old.ll ," LL_new ", new.ll, " diff.LL ", maxdelta ," params: \n")
        print_params( params, modes )

        if( steps == 1 ){
            diff[1] = maxdelta
            diff[2] = maxdelta
        } else {
            temp    = diff[1]
            diff[1] = maxdelta
            diff[2] = temp
        }

        #--- EM STOPPING CONDITIONS ---# 
        
        #--- test number of iterations taken
        steps = steps + 1
        if( steps > max.steps ){
            emTEST1 = 1
            cat("Stopping because steps >", max.steps, "!\n")
            break
        }

         #--- test if our current ll is nan?
        if( is.nan(new.ll) ){
            emTEST2 = 1
            cat("Stopping because ll = nan!\n")
            break
        }

        #--- test if old.ll < new.ll?
        if( old.ll > new.ll ){
            emTEST3 = 1
            cat("Stopping because old.ll ", old.ll , " > new.ll ", new.ll ,"!\n")
            break
        }
        
        #--- EM STOPPING CONDITIONS ---# 

        
    }#end em algo

    
        
    #--- return result of running em alg
    return(new("EM",
               params=params,
               Qij=Qij,
               steps=steps,
               ll=new.ll))
    
}
 
    
run.em <- function( Adj, obs, meas, tol=1e-5, max.steps=1e3, restarts=10, initPARAMS=NULL, fixPARAMS=NULL, epsilon=NULL ){

    #--- No: of nodes in graph    
    n <- length(V(Adj))
    
    #--- No: of edges in graph                  
    m <- length(E(Adj))

    modes <- length(obs)
    
    #nmeas <- Nmeas[study]
    
    #--- format edge measurements     
    #Nij <- list2matrix( meas[study], NN )
    #Nij = ifelse(Nij > 0, nmeas,0) ##test
    #Nij <- matrix(nmeas,NN,NN)
    #Nij <- nmeas
    Nij = matrix(0, ncol=modes, nrow=1)
    for( k in 1:modes ){
        Nij[1,k] = meas[[k]]
    }
    
    #--- format edge counts
    #Eij <- list2matrix( obs[study], n )
    #Eij <- list2matrix( obs, n )
    Eij = list()
    for( k in 1:modes ){
        Eij[[k]] <- matrix(as.vector(unlist(obs[[k]])),ncol=n,byrow=FALSE)
    }
    
    #--- model parameters
    # rho   ==> prior probablity, the existance of an edge in any position
    # m/choose(n,2) ==> most likely rho value for Bernoulli model given network 
    # alpha ==> true-positive rate, prob observing an edge where one truely exists
    # beta  ==> false-positive rate, prob observing an edge where none exists
    #---
    
    rho   = 0
    alpha = matrix(0, ncol=modes, nrow=1)
    beta  = matrix(0, ncol=modes, nrow=1)
        
    #--- the model parameters
    params <- list()

    params[[1]]      <- rho
    names(params)[1] <- "rho"
    params[[2]]      <- alpha
    names(params)[2] <- "alpha"
    params[[3]]      <- beta
    names(params)[3] <- "beta"

    if( !is.null(initPARAMS) ){
         params[[1]] = initPARAMS[[which(names(initPARAMS)=="rho")]]
         params[[2]] = initPARAMS[[which(names(initPARAMS)=="alpha")]]
         params[[3]] = initPARAMS[[which(names(initPARAMS)=="beta")]]
         if( is.null(fixPARAMS) ){ fixPARAMS = c(0,0,0) }
    }

    #--- Qij
    Q = matrix(0, n, n)
    
    res <- em( Adj=Adj, Nij=Nij, Eij=Eij, Qij=Q, params=params, modes=modes,
               tol=tol, max.steps=max.steps, initPARAMS=initPARAMS, fixPARAMS=fixPARAMS, epsilon=epsilon )
    
    #--- the parameters to retrun
    save.steps = res@steps;
    save.params= res@params;
    save.Q     = res@Qij;
    save.ll    = res@ll;    

    cat("init.ll ", save.ll, "\n")
    
    for( r in 1:restarts ){

        res <- em( Adj=Adj, Nij=Nij, Eij=Eij, Qij=Q, params=params, modes=modes, tol=tol, max.steps=max.steps, initPARAMS=initPARAMS, fixPARAMS=fixPARAMS, epsilon=epsilon )
        
        cat("saved.ll ", save.ll, " new.ll ", res@ll, "\n")
        
        if( abs(res@ll) < abs(save.ll) && is.finite(res@ll) ){
            #--- the parameters to retrun
            save.steps = res@steps;
            save.params= res@params;
            save.Q     = res@Qij;
            save.ll    = res@ll;
        }

    }

        
    #--- return result of running em alg
    return(new("runEM",
               Eij=Eij,
               Nij=Nij,
               params=save.params,
               Qij=save.Q,
               steps=save.steps,
               ll=save.ll))
    
}


    
