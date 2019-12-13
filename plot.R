require(ggplot2)

df <- read.delim("datapoints_TestEM.csv",sep="\t",header=T)

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

    ggsave("TestEM.pdf", width = 20, height = 20, units = "cm")




