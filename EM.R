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
                           ll="numeric",
                           delta.diff="numeric",
                           converage="numeric") )			   


#--- return type for em algo
setClass(Class="EM",representation(
                        params="list",
                        Qij="matrix",
                        steps="numeric",
                        ll="numeric",
                        delta.diff="numeric",
                        converage="numeric",
                        emTESTS="vector") )			   

#--- return edge odds ratio
setClass(Class="OR",representation(
                        norm="list",
                        term1="list",
                        term2="list") )


#--- return em stopping conditions
setClass(Class="CND",representation(
                        emTESTS="vector",
                        delta.diff="vector") )


print_params <- function( theta, modes ){

    rho   = theta[[which(names(theta)=="rho")]]
    alpha = theta[[which(names(theta)=="alpha")]]
    beta  = theta[[which(names(theta)=="beta")]]

    
    header = sprintf("\t mode %d \t",seq(1,modes,1))
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

fix_all_params <- function( theta ){

    fix <- list()    
    N   <- length(names(theta))
    for( i in 1:N ){
        m = length(theta[[i]])
        fix[[i]] = rep(1,length=m) 
    }
            
    return(fix)
    
}

use_param_constraints <- function( theta, const_min, const_max ){

    if( is.na(const_min) ){ const_min = 0 }
    if( is.na(const_max) ){ const_max = 1 }

    if( is.na(theta) )    { return( runif(1, const_min, const_max) ) }
    else{

        if( (theta - const_min) < 0 ){ const_min = 0 }
        else {                         const_min = theta - const_min }
        if( (theta + const_max) > 1 ){ const_max = 1 }
        else {                         const_max = theta + const_max }

        return( runif(1, const_min, const_max) ) 

    }

    return(runif(1))
    
}

#--- init model parameters to random values
init_params <- function( theta, fix, const ){

    N <- length(names(theta))
    for( i in 1:N ){
        for( m in 1:length(theta[[i]]) ){
            if( !fix[[i]][m] ){
                theta[[i]][1,m] <- use_param_constraints(as.numeric(theta[[i]][1,m]),
                                                         as.numeric(const[[i]][1,m]),
                                                         as.numeric(const[[i]][2,m]) )
            }
        }
    }
           
    return( theta )
}


init_qij <- function( gg, QQ ){

    n       = length(V(gg))
    m       = length(E(gg))
    ed      = get.edgelist(gg)
    tmp     = matrix(0,n,n)
    tmp[ed] = runif(m)

    return(tmp)
    
}

test_ll_convergance <- function( old.ll, new.ll ){

    maxdelta = abs( abs(old.ll) - abs(new.ll) )

    return(maxdelta)
    
}

test_param_convergance <- function( old.theta, theta ){

    maxdelta = 0
    
    N <- length(names(old.theta))
    for( i in 1:N ){
        if( is.matrix(old.theta[[i]]) ){            
            nr = nrow(old.theta[[i]])
            nc = ncol(old.theta[[i]])
            for( r in 1:nr ){
                for( c in 1:nc ){
                    deltap = abs( theta[[i]][r,c] - old.theta[[i]][r,c])
                    if( deltap > maxdelta ){ maxdelta = deltap }
                }                
            }
        } else {
            nr = length(old.theta[[i]])
            for( r in 1:nr ){
                deltap = abs( theta[[i]][r] - old.theta[[i]][r])
                if( deltap > maxdelta ){ maxdelta = deltap }
            }
        }
    }
    

    return(maxdelta)         

}

em_conditions <- function( max.steps, steps, ll, delta.diff, maxdelta, emTESTS ){

    store.delta.N = length(delta.diff)

    #--- test number of iterations taken
    if( steps > max.steps ){
        emTESTS[1] = 1
        cat("Stopping because steps >", max.steps, "!\n")
    }

    #--- test if our current ll is nan?
    if( is.nan(ll) ){
        emTESTS[2] = 1
        cat("Stopping because ll = nan!\n")
    }

    #--- test if old.ll < new.ll?
    store.delta = grep(TRUE, is.na(delta.diff))
    if( length(store.delta) > 0 ){
        delta.diff[min(store.delta)] = maxdelta
    } else {
        if( (sum(diff(delta.diff) > 0) == (store.delta.N-1)) ){
            emTESTS[3] = 1
            cat("Stopping because diff.ll is increasing!\n")
        } else {
            tmp                      = rep(NA,store.delta.N)
            tmp[1:(store.delta.N-1)] = delta.diff[2:store.delta.N]
            tmp[store.delta.N]       = maxdelta
            delta.diff               = tmp
        }
    }
               
    return(new("CND",
               emTESTS=emTESTS,
               delta.diff=delta.diff))
    
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



adj.log <- function( scale.l=NULL, arg.l=NULL ){

    result = 0
    
    if( !is.null(scale.l) && !is.null( arg.l ) ){

        if( (scale.l == 0) && (arg.l == 0) ){
            result = 0
        } else{
            result = scale.l*log(arg.l)
        }        
    }

    return(result)
    
}

# return the log-likelihood of the posterior, i.e. P(A,theta | data)
# See ref [2] eqn S4
posterior <- function( gg, QQ, meas, obs, theta ){

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

     for( i in 1:n ){
        for( j in 1:n ){

            if( i < j ){
                                
                #qij = QQ[i,j]
                aij = Aij[i,j]               

                ll = ll + adj.log( scale.l=aij,     arg.l=rho     ) +
                          adj.log( scale.l=(1-aij), arg.l=(1-rho) )
                
                for( k in 1:modes ){

                    nij.m = Nij[1,k]
                    eij.m = Eij[[k]][i,j]
                    
                    alpha.m = alpha[1,k] 
                    beta.m  = beta[1,k]

                   ll = ll + adj.log( scale.l=(aij*eij.m),           arg.l=alpha.m    ) +
                             adj.log( scale.l=aij*(nij.m-eij.m),     arg.l=(1-alpha.m)) +
                             adj.log( scale.l=(1-aij)*eij.m,         arg.l=beta.m     ) +
                             adj.log( scale.l=(1-aij)*(nij.m-eij.m), arg.l=(1-beta.m) )

                }                    
            }

        }
     }


    return(ll)
    
}


edge_odds_ratio <- function( gg, QQ, meas, obs, theta ){

    n   <- length(V(gg))
    Aij <- get.adjacency(gg)

    modes <- length(obs)
    
    rho   = theta[[which(names(theta)=="rho")]]
    alpha = theta[[which(names(theta)=="alpha")]]
    beta  = theta[[which(names(theta)=="beta")]]

    norm           = list()
    norm[[1]]      = matrix(0, n, n)
    names(norm)[1] = "rho"

    term1         = list()
    for( k in 1:modes ){
        term1[[k]]      = matrix(0, n, n)
        names(term1)[1] = sprintf("term1_%d",k)
    }

    term2         = list()
    for( k in 1:modes ){
        term2[[k]]      = matrix(0, n, n)
        names(term2)[1] = sprintf("term2_%d",k)
    }
    
    for( i in 1:n ){
        for( j in 1:n ){

            #if( i < j ){

                norm[[1]][i,j] = rho / (1-rho)

                for( k in 1:modes ){
                    
                    nij.m = meas[1,k]#[i,j] #for the moment 
                    eij.m = obs[[k]][i,j]

                    alpha.m = alpha[1,k] 
                    beta.m  = beta[1,k]

                    term1[[k]][i,j] = ( alpha.m/beta.m )^(eij.m)
                    term2[[k]][i,j] = ( (1-alpha.m)/(1-beta.m) )^(nij.m-eij.m)

                }
            #}                
            
        }
    }

    return(new("OR",norm=norm, term1=term1, term2=term2))
    
}

em <- function( Adj, Nij, Eij, Qij, params, modes, tol, max.steps,
               initPARAMS, fixPARAMS, constPARAMS,
               store.delta.N, e.first, conv.params ){

    
    #--- No: of nodes in graph    
    n <- length(V(Adj))
    
    #--- No: of edges in graph                  
    m <- length(E(Adj))

    modes <- modes
    
    #--- Init-step: assign random values to the model's parameters
    if( is.null(initPARAMS) ){
        Qij    <- init_qij( gg=Adj, QQ=Qij )
        params <- m.step( n=n, QQ=Qij, meas=Nij, obs=Eij, theta=params)
    } else {
        e.first = FALSE
        params <- init_params(initPARAMS, fixPARAMS, constPARAMS) 
        Qij    <- e.step( n=n, QQ=Qij, meas=Nij, obs=Eij, theta=params)
    }

    #--- print expected param values
    cat("init params: \n")
    print_params( params, modes )
    
    #--- em alg. steps    
    steps = maxdelta = 1

    #--- em test condition flags  
    emTESTS = rep(0,length=3)


    if( max.steps < store.delta.N ){ store.delta.N = max.steps }
    
    delta.diff = rep(NA,length=store.delta.N)
    
    
    #run the em algo.
    while( maxdelta > tol ) {
    
        #--- record old ll for run
        if( !conv.params ){
            old.ll     <- posterior( gg=Adj, QQ=Qij, meas=Nij, obs=Eij, theta=params )
        }
        
        if( e.first ){        
        #--- expectation step
        Qij        <- e.step( n=n, QQ=Qij, meas=Nij, obs=Eij, theta=params)
        }
        
        #--- maximazation step
        old.params <- params
        params     <- m.step( n=n, QQ=Qij, meas=Nij, obs=Eij, theta=params)

        if( !e.first ){        
        #--- expectation step
        Qij        <- e.step( n=n, QQ=Qij, meas=Nij, obs=Eij, theta=params)
        }        
                
        #--- record new ll for run
        new.ll     <- posterior( gg=Adj, QQ=Qij, meas=Nij, obs=Eij, theta=params )
        
        
        #--- test convergence of parameters        
        if( conv.params ){
            maxdelta = test_param_convergance( old.params, params )
        } else {
            #--- test convergence of ll
            maxdelta = test_ll_convergance( old.ll, new.ll )
        }

        #--- print expected param values
        cat(" EM step ", steps, " , max change = ", maxdelta, "\n")
        cat(" params: \n")
        print_params( params, modes )
        
            
        #--- EM STOPPING CONDITIONS ---# 
        steps = steps + 1
        tests = em_conditions( max.steps=max.steps, steps=steps, ll=new.ll,
                               delta.diff=delta.diff, maxdelta=maxdelta,
                               emTESTS=emTESTS )
        if( sum(tests@emTESTS) > 0 ){ break; }
        emTESTS    = tests@emTESTS
        delta.diff = tests@delta.diff                       
        #--- EM STOPPING CONDITIONS ---# 

        
    }#end em algo

    converage = ifelse( maxdelta > tol, 1, 0) 
    
    #--- return result of running em alg
    return(new("EM",
               params=params,
               Qij=Qij,
               steps=steps,
               ll=new.ll,
               delta.diff=maxdelta,
               converage=converage,
               emTESTS=emTESTS))
    
}
 
run.em <- function( Adj, obs, meas, tol=1e-5, max.steps=1e3, restarts=10,
                   initPARAMS=NULL, fixPARAMS=NULL, constPARAMS=NULL,
                   store.delta.N=2, e.first=TRUE, conv.params=TRUE ){

    #--- No: of nodes in graph    
    n <- length(V(Adj))
    
    #--- No: of edges in graph                  
    m <- length(E(Adj))

    #--- No: of 
    modes <- length(obs)

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
    #Eij <- matrix(as.vector(unlist(obs)),ncol=n,byrow=FALSE)
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

    rho   = matrix(NA, ncol=1,     nrow=1)
    alpha = matrix(NA, ncol=modes, nrow=1)
    beta  = matrix(NA, ncol=modes, nrow=1)
        
    #--- the model parameters
    params    <- list()
    
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
         if( is.null(fixPARAMS) ){
             fixPARAMS <- list()
             fixPARAMS[[1]] = matrix(0,ncol=1, nrow=1)
             names(fixPARAMS)[1] = "rho"
             fixPARAMS[[2]] = matrix(0,ncol=modes, nrow=1)
             names(fixPARAMS)[2] = "alpha"
             fixPARAMS[[3]] = matrix(0,ncol=modes, nrow=1)
             names(fixPARAMS)[3] = "beta"
         }             
    }

    if( is.null(constPARAMS) ){
        constPARAMS <- list()
        constPARAMS[[1]] = matrix(NA,ncol=1,nrow=2)
        names(constPARAMS)[1] = "rho"
        constPARAMS[[2]]  = matrix(NA,ncol=modes,nrow=2)
        names(constPARAMS)[2] = "alpha"
        constPARAMS[[3]] = matrix(NA,ncol=modes,nrow=2)
        names(constPARAMS)[3] = "beta"
    }

    
    #--- Qij
    Q = matrix(0, n, n)

    foundMin = rep(FALSE, length=3)
    f=1
    
    res <- em( Adj=Adj, Nij=Nij, Eij=Eij, Qij=Q, params=params, modes=modes,tol=tol, max.steps=max.steps,
               initPARAMS=initPARAMS, fixPARAMS=fixPARAMS, constPARAMS=constPARAMS,
               store.delta.N=store.delta.N, e.first=e.first, conv.params=conv.params )
    
    #--- the parameters to retrun
    save.steps     = res@steps;
    save.params    = res@params;
    save.Q         = res@Qij;
    save.ll        = res@ll;
    save.diff      = res@delta.diff;
    save.converage = res@converage;

    cat("init.ll ", save.ll, "\n")
    
    for( r in 1:restarts ){

        res <- em( Adj=Adj, Nij=Nij, Eij=Eij, Qij=Q, params=params, modes=modes, tol=tol, max.steps=max.steps,
                  initPARAMS=initPARAMS, fixPARAMS=fixPARAMS, constPARAMS=constPARAMS,
                  store.delta.N=store.delta.N, e.first=e.first, conv.params=conv.params )

        #--- if we didn't converage, because we ran out of steps, should we retart res?
        if( res@emTESTS[1] == 1 ){
            cat("ran out of steps, try restarting...")
            res <- em( Adj=Adj, Nij=Nij, Eij=Eij, Qij=res@Qij, params=res@params, modes=modes, tol=tol, max.steps=max.steps,
                      initPARAMS=res@params, fixPARAMS=fix_all_params(res@params), constPARAMS=constPARAMS,
                      store.delta.N=store.delta.N, e.first=e.first, conv.params=conv.params )
            cat(" done.\n")
            res@steps = max.steps + res@steps
        }
        
        cat("saved.ll ", save.ll, " new.ll ", res@ll, "...")

        if( res@emTESTS[3] == 0 ){
            if( abs(res@ll) < abs(save.ll) ){
                cat(" will save new.ll result.\n")
                #--- the parameters to retrun
                save.steps     = res@steps;
                save.params    = res@params;
                save.Q         = res@Qij;
                save.ll        = res@ll;
                save.diff      = res@delta.diff;
                save.converage = res@converage;
            } else {
                if( abs(res@ll) == abs(save.ll) ){ foundMin[f] = TRUE; f=f+1; }
                cat(" will continue.\n")
            }
        }

        if( sum(foundMin) == length(foundMin) ){
            cat("keep finding the minimum so stop!\n")
            break;
        }
        
    }

    #--- return result of running em alg
    return(new("runEM",
               Eij=Eij,
               Nij=Nij,
               params=save.params,
               Qij=save.Q,
               steps=save.steps,
               ll=save.ll,
               delta.diff=save.diff,
               converage=save.converage))
    
}

    
