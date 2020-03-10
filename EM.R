# References:
# M. E. J. Newman, Network sturcture from rich but noisy data, nature physics, 2018.
# example EM algo structure: http://www.di.fc.ul.pt/~jpn/r/EM/EM.html

#library(methods)

#--- return type for em algo
setClass(Class="runEM",representation(
                           Eij="matrix",
                           Nij="numeric",
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


#--- init model parameters to random values
init_params <- function( theta, fix ){

    N <- length(names(theta))
    for( i in 1:N ){
        if( !fix[i] ){
            theta[[i]] <- runif(1)
        }
    }
           
    return( theta )
}


init_qij <- function( QQ ){

    NROW = nrow(QQ)
    NCOL = ncol(QQ)
    
    return(matrix(runif(NROW*NCOL),NROW,NCOL))
    
}


Qij.param <- function( meas, obs, theta ){

    N     = meas
    Eij   = obs
    
    rho   = theta[[which(names(theta)=="rho")]]
    alpha = theta[[which(names(theta)=="alpha")]]
    beta  = theta[[which(names(theta)=="beta")]]
    
    num = rho * alpha^(Eij) * (1-alpha)^(N-Eij)
    dem = num + (1-rho) * beta^(Eij) * (1-beta)^(N-Eij)

    Qij = 0

    if( dem != 0 ){ Qij = num / dem }
    
    return(Qij)
    
}

# The e-step calculates the prosterior distribution over networks, given 
# our model parameters, theta, from the m-step.
e.step <- function( gg, QQ, meas, obs, theta){

    n <- length(V(gg))

    for( i in 1:n ){
        for( j in 1:n ){

            if( i < j ){
                QQ[i,j] = Qij.param( meas=meas, obs=obs[i,j], theta=theta )                
            }

        }
    }

    return(QQ)
    
}


m.step <- function( gg, QQ, meas, obs, theta ){

    n  <- length(V(gg))

    N  <- meas
    
    h_alpha   = 0
    h_beta    = 0
    h_rho     = 0
        
    num_alpha = 0
    num_beta  = 0
    #num_rho   = 0
    
    dem_alpha = 0
    dem_beta  = 0
    dem_rho   = 0

    for( i in 1:n ){
        for( j in 1:n ){

            if( i < j ){
                                
                eij = obs[i,j]
                qij = QQ[i,j]                
                
                num_alpha = num_alpha + (eij * qij)
                dem_alpha = dem_alpha + qij

                num_beta  = num_beta  + (eij * (1-qij))                
                dem_beta  = dem_beta  + (1-qij)
                
                #num_rho   = num_rho   + qij

            }

            
        }
    }

    if( dem_alpha != 0 ){
        h_alpha = num_alpha / (N*dem_alpha)
    }

    if( dem_beta != 0 ){
        h_beta  = num_beta / (N*dem_beta)
    }        
    
    dem_rho = choose( n, 2 ) 
    #OR, i.e. just the number of elements in the upper triangle,
    #    excluding the diagonal 
    #dem_rho = (n*(n-1))/2
    
    if( dem_rho != 0 ){
        h_rho   = dem_alpha / dem_rho
        #h_rho   = num_rho / dem_rho
    }
        
    theta[[which(names(theta)=="rho")]]   = h_rho
    theta[[which(names(theta)=="alpha")]] = h_alpha
    theta[[which(names(theta)=="beta")]]  = h_beta
   
    #--- return result of running m.step
    return(theta)
    
}

#  return the log-likelihood of the posterior, i.e. q(A)
## ll  = ll + aij * log( qij ) + (1-aij) * log( (1-qij) )
posterior <- function( gg, QQ,  epsilon=NULL){

    n   <- length(V(gg))

    Aij <- get.adjacency(gg)
    
    #log likelihood
    ll   <- 0

    if( is.null(epsilon) ){
        epsilon = .Machine$double.eps
    }
    
     for( i in 1:n ){
        for( j in 1:n ){

            if( i < j ){
                                
                qij = QQ[i,j]
                aij = Aij[i,j]

                if( qij == 0 || qij == 1 ){
                    if( qij == 0 ){
                        ll = ll + aij * log( epsilon ) + (1-aij) * log( 1 )
                    }
                    if( qij == 1 ){
                        ll = ll + aij * log( 1 ) + (1-aij) * log( (1-epsilon) )
                    }
                } else {
                    ll  = ll + aij * log( qij ) + (1-aij) * log( (1-qij) )
                }
            }

        }
     }


    return(ll)
    
}



em <- function( Adj, Nij, Eij, Qij, params, tol, max.steps, initPARAMS, fixPARAMS ){

    #--- Init-step: assign random values to the model's parameters
    if( is.null(initPARAMS) ){
        Qij    <- init_qij( Qij )
        params <- m.step( gg=Adj, QQ=Qij, meas=Nij, obs=Eij, theta=params)
    } else {
        params <- init_params(initPARAMS, fixPARAMS) 
        Qij    <- e.step( gg=Adj, QQ=Qij, meas=Nij, obs=Eij, theta=params)
    }
    
    cat("init params ", unlist(params),"\n")
 
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
        old.ll     <- posterior( gg=Adj, QQ=Qij, epsilon=10^-9 )

        #--- maximazation step
        old.params <- params
        params     <- m.step( gg=Adj, QQ=Qij, meas=Nij, obs=Eij, theta=params)

        
        #--- expectation step
        Qij        <- e.step( gg=Adj, QQ=Qij, meas=Nij, obs=Eij, theta=params)
        
                
        #--- record new ll for run
        new.ll     <- posterior( gg=Adj, QQ=Qij, epsilon=10^-9 )
        
        #--- test convergence of parameters        
        maxdelta = abs( abs(old.ll) - abs(new.ll) )
        cat(" LL_old ", old.ll ," LL_new ", new.ll, " diff.LL ", maxdelta ," params ", unlist(params)  ,"\n")

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
 
    
run.em <- function( Adj, obs, Nmeas=c(4,8,16), study=1, errRATE, tol=1e-5, max.steps=1e3, restarts=10, initPARAMS=NULL, fixPARAMS=NULL ){

    #--- No: of nodes in graph    
    n <- length(V(Adj))
    
    #--- No: of edges in graph                  
    m <- length(E(Adj))

    nmeas <- Nmeas[study]
    
    #--- format edge measurements     
    #Nij <- list2matrix( meas[study], NN )
    #Nij = ifelse(Nij > 0, nmeas,0) ##test
    #Nij <- matrix(nmeas,NN,NN)
    Nij <- nmeas
    
    #--- format edge counts
    #Eij <- list2matrix( obs[study], n )
    #Eij <- list2matrix( obs, n )
    Eij <- matrix(as.vector(unlist(obs)),ncol=n,byrow=FALSE)

    #--- model parameters
    # rho   ==> prior probablity, the existance of an edge in any position
    # m/choose(n,2) ==> most likely rho value for Bernoulli model given network 
    # alpha ==> true-positive rate, prob observing an edge where one truely exists
    # beta  ==> false-positive rate, prob observing an edge where none exists
    #---
    
    rho   = 0
    alpha = 0
    beta  = 0
        
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
    
    res <- em( Adj=Adj, Nij=Nij, Eij=Eij, Qij=Q, params=params, 
               tol=tol, max.steps=max.steps, initPARAMS=initPARAMS, fixPARAMS=fixPARAMS )
    
    #--- the parameters to retrun
    save.steps = res@steps;
    save.params= res@params;
    save.Q     = res@Qij;
    save.ll    = res@ll;    

    cat("init.ll ", save.ll, "\n")
    
    for( r in 1:restarts ){

        res <- em( Adj=Adj, Nij=Nij, Eij=Eij, Qij=Q, params=params, 
                  tol=tol, max.steps=max.steps, initPARAMS=initPARAMS, fixPARAMS=fixPARAMS )
        
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

    
