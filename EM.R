# References:
# M. E. J. Newman, Network sturcture from rich but noisy data, nature physics, 2018.
# example EM algo structure: http://www.di.fc.ul.pt/~jpn/r/EM/EM.html

# Generate set of random measurements from test graphs 
cat('Generate oberserved and measured datastes...\n')
#source('Generate_Rndm_Measurements.R')
cat('done!\n')

library(methods)

#--- return type for em algo
setClass(Class="runEM",representation(
                           Eij="matrix",
                           Nij="matrix",
                           params="list",
                           Qij="matrix",
                           steps="numeric",
                           ll="numeric")
         )			   


#--- return type for em algo
setClass(Class="EM",representation(
                           params="list",
                           Qij="matrix",
                           steps="numeric",
                           ll="numeric")
         )			   

#--- return type for M.step
setClass(Class="M",representation(
                       params="list",
                       ll="numeric")
         )			   


fix_params <- function( theta, pnames=c("alpha","rho"), fixedValue ){

    for( i in 1:length(pnames)){
        indx <- which(names(theta)==pnames[i])
        if( length(indx) != 0 ){
            theta[indx] = fixedValue
        }
        
    }

    return(theta)
    
}

#--- init model parameters to random values
init_params <- function( theta ){

    N <- length(names(theta))
    for( i in 1:N ){
        theta[[i]] <- runif(1)
    }
           
    return( theta )
}

init_qij <- function( QQ ){

    NROW = nrow(QQ)
    NCOL = ncol(QQ)
    
    return(matrix(runif(NROW*NCOL),NROW,NCOL))
    
}

e.step <- function( GG, QQ, meas, obs, theta){

    NN <- length(V(GG))

    #Adj <- get.adjacency(GG)
    
    h_rho   = theta[[which(names(theta)=="rho")]]
    h_alpha = theta[[which(names(theta)=="alpha")]]
    h_beta  = theta[[which(names(theta)=="beta")]]

    for( i in 1:NN ){
        for( j in 1:NN ){

            if( i < j ){

                eij = nij = num = dem = val = 0
                
                eij     = obs[i,j]
                nij     = meas[i,j]
                num     = h_rho * (h_alpha^eij) *(1 - h_alpha)^(nij-eij)
                dem     = num + (1 - h_rho)*(h_beta^eij) * (1 - h_beta)^(nij-eij)

                if( dem != 0 && !is.na(dem) ){
                    val     = num / dem
                    if( !is.na(val) ){
                        QQ[i,j] = val
                        QQ[j,i] = val
                    }
                }

                
            }
        }
    }

    return(QQ)
    
}

m.step <- function( GG, QQ, meas, obs, theta){

    NN  <- length(V(GG))
    Adj <- get.adjacency(GG)

    
    h_alpha   = 0
    h_beta    = 0
    h_rho     = 0
    
    num_alpha = 0
    num_beta  = 0
    num_rho   = 0
    
    dem_alpha = 0
    dem_beta  = 0
    dem_rho   = 0

    new_ll    = 0
    
    for( i in 1:NN ){
        for( j in 1:NN ){

            if( i < j ){
                                
                eij = nij = qij = aij = term1 = term2 = 0;
                
                eij = obs[i,j]
                nij = meas[i,j]
                qij = QQ[i,j]
                aij = Adj[i,j]
                
                num_alpha = num_alpha + (eij * qij)
                dem_alpha = dem_alpha + (nij * qij)

                num_beta  = num_beta  + (eij * (1-qij))                
                dem_beta  = dem_beta  + (nij * (1-qij))
                
                num_rho   = num_rho   + qij

                term1  = qij^aij
                term2  = (1-qij)^(1-aij)
                new_ll = new_ll + log( term1*term2 )
                
            }

            
        }
    }

    if( dem_alpha != 0 ){
        h_alpha = num_alpha / dem_alpha
    }

    if( dem_beta != 0 ){
        h_beta  = num_beta / dem_beta 
    }        
    
    dem_rho = choose( NN, 2 ) 
    #OR, i.e. just the number of elements in the upper triangle,
    #    excluding the diagonal 
    #dem_rho = (NN*(NN-1))/2
    
    if( dem_rho != 0 ){
        h_rho   = num_rho / dem_rho
    }
        
    theta[[which(names(theta)=="rho")]]   = h_rho
    theta[[which(names(theta)=="alpha")]] = h_alpha
    theta[[which(names(theta)=="beta")]]  = h_beta
   
    #--- return result of running em alg
    return(new("M",
               params=theta,
               ll=(new_ll)))#NOTE: -new_ll, inverts Qij, i.e. 1-Qij??
    
    
}


prosterior <- function( GG, QQ ){

    NN  <- length(V(GG))
    Adj <- get.adjacency(GG)
    
    ll = 0;
    
    for( i in 1:NN ){
        for( j in 1:NN ){

            if( i > j ){
                term1 = QQ[i,j]^Adj[i,j]
                term2 = (1-QQ[i,j])^(1-Adj[i,j])
                ll    = ll + log( term1*term2 )
                }
                    
        }
    }

    return(ll)
    
}

Mbind <- function(nij, eij, qij, adj){

    indx <- upper.tri(nij)
    oo <- cbind(as.numeric(nij[indx]),as.numeric(eij[indx]),as.numeric(qij[indx]),
                as.numeric(adj[indx]))
    return(oo)
    
}

em <- function( Adj, Nij, Eij, Qij, params, max.ll, fix.params, fix.params.val, tol, max.steps, g.steps ){

    #--- Init-step: assign random values to the model's parameters
    params <- init_params( params )
    #params <- fix_params( params, pnames=c("beta"),fixedValue=(1-params[[which(names(params)=="alpha")]] ))
  #params <- fix_params( params, pnames=c("rho"),fixedValue=0.5 )
    
  Qij    <- init_qij( Qij )
    
  cat("init params ", unlist(params),"\n")

  gradient <- rep(0,length=g.steps)  
  g = 1
    
  #--- em alg. steps    
  steps = maxdelta = 1

  #--- em test flags  
  emTEST1 = 0
  emTEST2 = 0
  emTEST3 = 0
    
    #run the em algo.
    while( maxdelta > tol ) {

        #--- expectation step
        Qij        <- e.step( Adj, Qij, Nij, Eij, params)
        old.params <- params

        #--- maximazation step
        m.res      <- m.step( Adj, Qij, Nij, Eij, params)
        params     <- m.res@params
        new.ll     <- m.res@ll

        #--- recorded the highest ll for run
        #if( steps == 1 ){
        #    saved.params = params
        #    saved.Qij    = Qij
        #    saved.ll     = new.ll
        #    saved.steps  = steps
        #} else {
        #    if( new.ll > saved.ll ){
        #        saved.params = params
        #        saved.Qij    = Qij
        #        saved.ll     = new.ll
        #        saved.steps  = steps
        #    }
            
        #}

        
        #--- test convergence of parameters        
        maxdelta = max(abs(unlist(old.params)-unlist(params)),na.rm=T)
                 
        cat("maxdelta ", maxdelta, " LL ", new.ll ," params ", unlist(params),"\n")


        #--- EM STOPPING CONDITIONS ---# 
        #--- record gradient of ll's for the last g.steps
        gradient[g] = as.numeric(new.ll)

        ##--- test gradient of the ll's for the last g.steps
        if( g %% g.steps == 0 ){
            if( (gradient[g] - gradient[1]) < 0 ){
                cat("Stopping because gradient is still negative!\n")
                emTEST1 = 1
                break
            } else {
                g = 1
                gradient    <- rep(0,length=g.steps)
                gradient[g] <- new.ll 
            }
        }

        g = g + 1
        #---
        
        
        #--- test number of iterations taken
        steps = steps + 1
        if( steps > max.steps ){
            emTEST2 = 1
            cat("Stopping because steps >", max.steps, "!\n")
            break
        }

        #--- test if our current ll < the max.ll
        if( new.ll < max.ll ){
            emTEST3 = 1
            cat("Stopping because ll <", max.ll, "!\n")
            break
        }
        #--- EM STOPPING CONDITIONS ---# 

        
    }#end em algo

    #--- if we stopped naturally, i.e. converage to highest ll, use this run.
    #if( (maxdelta < tol) && new.ll > saved.ll ){
    #    saved.params = params
    #    saved.Qij    = Qij
    #    saved.ll     = new.ll
    #    saved.steps  = steps
    #} 
    
    #--- return result of running em alg
    return(new("EM",
               params=params,
               Qij=Qij,
               steps=steps,
               ll=new.ll))
    
}
    
run.em <- function( Adj, meas, obs, Nmeas=c(4,8,16), study=1, errRATE, tol=1e-3, max.steps=1e3, restarts=10, g.steps=20 ){

#--- No: of nodes in graph    
NN <- length(V(Adj))
#--- No: of edges in graph                  
MM <- length(E(Adj))

nmeas <- Nmeas
nmeas <- nmeas[study]    
    
#--- format edge measurements     
#Nij <- list2matrix( meas[study], NN )
#Nij = ifelse(Nij > 0, nmeas,0) ##test
Nij <- matrix(nmeas,NN,NN)
    
#--- format edge counts
Eij <- list2matrix( obs[study], NN )

#Eij <- get.adjacency(Adj)
#Eij = as.matrix(Eij)
#Eij[Eij == 1] = nmeas
    
#Eij <- matrix(0,NN,NN)
#for( m in 1:MM ){
#    indx <- as.vector(unlist(ends(Adj,m)))
#    Eij[indx[1],indx[2]] = obs[(m+1),(study+1)]
#    Eij[indx[2],indx[1]] = obs[(m+1),(study+1)]
# }


#--- prior probablity, the existance of an edge in any position
rho = 0

#--- true-positive rate, prob observing an edge where one truely exists 
alpha = 0

#--- false-positive rate, prob observing an edge where none exists
beta = 0

#--- log-likeihood
ll   = -1e9
    
#--- the model parameters
params <- list()

params[[1]]      <- rho
names(params)[1] <- "rho"
params[[2]]      <- alpha
names(params)[2] <- "alpha"
params[[3]]      <- beta
names(params)[3] <- "beta"

#--- Qij
Q = matrix(0, NN, NN)
    
    res <- em( Adj=Adj, Nij=Nij, Eij=Eij, Qij=Q, params=params, max.ll=ll,
              fix.params=c("rho","alpha"), fix.params.val=(1-errRATE),
              tol=tol, max.steps=max.steps, g.steps=g.steps )
    
#--- the parameters to retrun
save.steps = res@steps;
save.params= res@params;
save.Q     = res@Qij;
save.ll    = res@ll;    

cat("init.ll ", save.ll, "\n")
    
    for( r in 1:restarts ){

        res <- em( Adj=Adj, Nij=Nij, Eij=Eij, Qij=Q, params=params, max.ll=save.ll,
                  fix.params=c("rho","alpha"), fix.params.val=(1-errRATE),
                  tol=tol, max.steps=max.steps, g.steps=g.steps )

        cat("old.ll ", save.ll, " new.ll ", res@ll, "\n")
        
        if( res@ll > save.ll && is.finite(res@ll) ){
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

    
#
#--- run the em alg.
run=1

if( run ){

oo     <- matrix(0,ncol=2,nrow=length(errRates))
oo[,1] <- errRates 

s= 1

saved <- list()
    
for( p in 1:length(errRates) ){

    
cat("\n")
cat("errRate = ", errRates[p],"\n")
    
res <- run.em( Adj=Aij[[p]], meas=Ecount[[p]], obs=Ecount[[p]], Nmeas=c(8),study=s, errRATE=errRates[p], tol=1e-4, max.steps=100,g.steps=10 )

    saved[[p]] <- res
    names(saved)[p] <- sprintf("errRate_%#.2f",errRates[p])
    
    adj <- get.adjacency(Aij[[p]])
    
    xx     <- Mbind(res@Nij,res@Eij,res@Qij,adj)
    RECALL <- table(xx[xx[,4]==1,3] > 0.5)/M

    recall <- NA
    indx <- which(names(RECALL) == "TRUE")
    if( length(indx) != 0 ){
        recall <- as.vector(RECALL[indx])
    }
    
    cat("errRate = ", errRates[p], ", params = [" , unlist(res@params), "], steps = ", res@steps , " LL = ", res@ll , ", recall = ", recall,"\n")

    oo[p,2] <- recall
    
}

        
}
