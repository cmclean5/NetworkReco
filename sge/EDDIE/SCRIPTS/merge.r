#--- Merge result files ---#

rm(list=ls())

#library(igraph)

#load corresponding graph
#---Directories
DIR    <- vector(length=3)
DIR[1] <- "../../datasets/"
DIR[2] <- "RESULTS"
DIR[3] <- "OUT"

files    <- vector(length=1)
files[1] <- "recall.csv"  

subdirs  = list.files(path=sprintf("../%s/",DIR[2]));
nstudies = length(subdirs);


 for( f in 1:length(files) ){

    st1 = sprintf("../%s/%s/%s",DIR[2],subdirs[1],files[f]);	
    if( file.exists(st1) && file.info(st1)$size!=0 ){
     tb = read.table(st1,header=T,sep="\t");
    }

     for( s in 2:nstudies ){

      st1 = sprintf("../%s/%s/%s",DIR[2],subdirs[s],files[f]);	

      if( file.exists(st1) && file.info(st1)$size!=0 ){
      	  temp <- read.table(st1,header=T,sep="\t");
          tb <- rbind(tb,temp)
      }
     }

     outfile <- file(sprintf("%s/%s",DIR[3],files[f]),"w");
     write.table(tb, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
     close(outfile);

     rm(tb)
}



