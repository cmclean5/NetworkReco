require(ggplot2)

reshape.df <- function( df ){

    N    <- length(df[,1])    
    tmp  <- cbind( df[,1], df[,2], rep("N4",length=N))
    tmp2 <- cbind( df[,1], df[,3], rep("N8",length=N))
    tmp3 <- cbind( df[,1], df[,4], rep("N16",length=N))
    res  <- rbind(tmp,tmp2,tmp3)
    
    colnames(res) <- c("errorRate","Recall","Nmeas")

    return(res)
    
}

file <- "recall.csv"

df <- read.delim(sprintf("%s",file),sep="\t",header=T)
df <- reshape.df( df )
DF <- as.data.frame(df)

colours <- c('lawngreen','firebrick2',"royalblue2")

gplot <- ggplot(DF,aes(x=as.numeric(as.vector(DF[,1])),y=as.numeric(as.vector(DF[,2])),color=as.vector(DF[,3]),group=as.vector(DF[,3])))+
    geom_point()+
    geom_line()+
    labs(x="error Rate",y="Recall",title="")+
    theme(            
        axis.title.x=element_text(face="bold",size=rel(2)),
        axis.text.x =element_text(face="bold",size=rel(2)), 
        axis.title.y=element_text(face="bold",size=rel(2)),
        axis.text.y =element_text(face="bold",size=rel(2)), 
        legend.title=element_text(face="bold",size=rel(1.5)),
        legend.text=element_text(face="bold",size=rel(1.5)),
        legend.key=element_blank())+
    theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
          panel.grid.minor = element_line(colour="grey40",size=0.1),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype="solid",fill=NA))+
    scale_color_manual("",breaks=levels(factor(DF[,3])),values=c(colours))+
    scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2))+
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2))

    ggsave("simulatedData.pdf", width = 20, height = 20, units = "cm")




