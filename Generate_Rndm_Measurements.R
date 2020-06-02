# References:
# M. E. J. Newman, Network sturcture from rich but noisy data, nature physics, 2018.
# example EM algo structure: http://www.di.fc.ul.pt/~jpn/r/EM/EM.html


#--- reset
#rm(list=ls())

#--- imports/headers
#require(igraph)
  
#--- set random number seed
setSeed <- function(SEED){
    if( is.null(SEED) ){
        SEED = as.integer(Sys.time())
        set.seed(SEED)
    } else {
        set.seed(SEED)
    }
}


list2matrix <- function( LL, NN, byROW=FALSE ){

    MM <- matrix(as.vector(unlist(LL)),ncol=N,byrow=byROW)
    return(MM)

}


realEdge <- function( GG, ED ){

    indx <- get.edge.ids(GG,ED)
    return(indx)
    
}

generateMeasurements <- function( GG, GGname, Nrand, errorRate, Nmeas=c(4,8,16), SEED=NULL  ){

    setSeed(SEED)
    
    Ntot  <- length(V(GG))
    Mtot  <- length(E(GG))

    #alpha <- 1 - errorRate
    
    Mrand <- floor(Mtot * errorRate)
    
    Lmeas <- length(Nmeas)    

    Gen   <- 1000
    
    SETmeas <- list()    
    
    for( s in 1:Lmeas ){

        Ntt <- matrix(0,Ntot,Ntot)

        for( g in 1:Gen ){

        if( Mrand == 0 ){
            rd <- GG
        } else if( Mrand == Mtot ){
            rd <- sample_gnm(n=Ntot,m=Mtot,directed=F,loops=F)
        } else {

            #--- introduces false negative errors, i.e. 
            #    randomly remove edges from graph
            rd <- delete_edges(GG,sample(E(GG),Mrand))
    
            
            #--- introduces false postive errors, i.e. 
            #    spurious edges (i.e. where none existed before) into the graph
            rd <- rd + igraph::edges(sample(V(rd),2*Mrand,replace=T))
                
        }
            
        RAD <- as.matrix(get.adjacency(rd))        

        count = Nmeas[s]

            while( count > 0 ){

                si <- sample(1:Ntot,1,replace=T)
                sj <- sample(1:Ntot,1,replace=T)

                if( si != sj ){
                    aij  = RAD[si,sj]
                    prob = ifelse(aij>0,1,0)
                    if( (prob == 1) && (Ntt[si,sj] < Nmeas[s]) ){
                        Ntt[si,sj] = Ntt[si,sj] + prob
                        Ntt[sj,si] = Ntt[si,sj]
                        count = count - 1
                    }
                }

            }
            
            rm(rd)
            
        }

        #for( i in 1:Ntot ){
        #    for( j in 1:Ntot ){

        #        if( i < j ){

        #            if( realEdge( GG, c(i,j)) != 0 ){
        #                prob = 1
        #                Ntt[i,j] = Nmeas[s]*Gen
        #                Ntt[j,i] = Ntt[i,j]
        #            }
        #        }
                
        #    }
        #}
        
        #SETmeas[[s]] <- round(Ntt/Gen)
        SETmeas[[s]] <- Ntt
        names(SETmeas)[s] <- sprintf("%s_Rnd_%d_Nmeas_%s",GGname, Nrand, Nmeas[s])

    }
       
   
    return(SETmeas)
    
}

Test <- function( GG, MM ){

   NN   <- length(V(GG))
   adj  <- as.matrix(get.adjacency(GG))
   indx <- upper.tri(adj)
   mm1  <- list2matrix(MM,NN)

    xx   <- cbind(adj[indx],mm1[indx])

    return(xx)
    
}


  
