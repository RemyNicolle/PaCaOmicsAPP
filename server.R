library(shiny)
library(beeswarm)


#source("getExpGroupsFunction.R")
load("shinyData.RData")

momcolors= c( basal="#ff7f00",classical="#377eb8",Outlier="#636363")
tumostromcol=c(tumor="#89000B",stroma="#148D4F")


#H2M<<-cit.load("/datacit/00_DATABANKS/HumanMouseHomologs/H2M.RData")

#"Tumor","Stroma","TumorVsStroma"
# Mrnadiff,Mrna,Hrnadiff,Hrna,fullclassif,conorm,humanmousefac,StromaTumDiff

#selectionmethod=list(totalReads=sum,SD=sd)

L=list(
Tumor=list(diff=Hrnadiff,exp=Hrna,cl=fullclassif,col=momcolors,genecol="GeneName"),
Stroma=list(diff=Mrnadiff,exp=Mrna,cl=fullclassif,col=momcolors,genecol="GeneName"),
TumorVsStroma=list(diff=StromaTumDiff,exp=conorm,cl=humanmousefac,col=tumostromcol,genecol="HumanGene"))


genebyclass=function(gene,dataname){
    DATA=L[[dataname]]
    
    
    
    if(gene %in% DATA$diff[,DATA$genecol]){
        selectgenesym=gene
        selectgeneEnsID=rownames(DATA$diff)[which(DATA$diff[,DATA$genecol]==selectgenesym)]
    }else{
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("No gene called",gene,"found."),cex = 1.6, col = "black")
        return(NULL)

    }
    
    
    if(length(selectgeneEnsID)>1){
        selectgeneEnsID=rownames(DATA$diff)[which.min(DATA$diff[selectgeneEnsID,"logFC"])]
        selectgenesym=DATA$diff[selectgeneEnsID,DATA$genecol]
    }
  
  
  if(dataname=="TumorVsStroma"){
  ylabel=paste("co-normalised expression",sep="")
  }else{
  ylabel=paste("Normalised log-count (",selectgeneEnsID,")",sep="")
  }
  beeswarm(as.numeric(DATA$exp[selectgeneEnsID,names(DATA$cl)])~DATA$cl,
  bty="l",pch=16,
  xlab=NA,
  ylab=ylabel,
  main=selectgenesym)
  
  
  boxplot(as.numeric(DATA$exp[selectgeneEnsID,names(DATA$cl)])~DATA$cl,add=T,axes=F,outline=F,
  col=DATA$col,
  sub=paste("log fold-change:",signif(DATA$diff[selectgeneEnsID,"logFC"],3),
  "\nFDR:",signif(DATA$diff[selectgeneEnsID,"adj.P.Val"],3)))
  
  beeswarm(as.numeric(DATA$exp[selectgeneEnsID,names(DATA$cl)])~DATA$cl,
  pch=16,
  add=T,axes=F,
  xlab=NA,
  ylab=NA,
  main=NA,cex=1.5)
  
}


shinyServer(function(input, output) {
    
    
    
    output$main_plot <- renderPlot({
        genebyclass( gene=toupper(input$gene),input$data)
        
    })
})
