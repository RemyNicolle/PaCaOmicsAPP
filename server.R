library(shiny)
library(beeswarm)


#source("getExpGroupsFunction.R")
load("shinyData.RData")

momcolors= c( basal="#ff7f00",classical="#377eb8")#,Outlier="#636363")
tumostromcol=c(tumor="#148D4F",stroma="#89000B")

#methdat,methdiff

RNAseqAllDat=list(
Tumor=list(diff=Hrnadiff,exp=Hrna,cl=noutlierclassif,col=momcolors,genecol="GeneName"),
Stroma=list(diff=Mrnadiff,exp=Mrna,cl=noutlierclassif,col=momcolors,genecol="HumanHomolog_Symbol"),
TumorVsStroma=list(diff=StromaTumDiff,exp=conorm,cl=humanmousefac,col=tumostromcol,genecol="HumanGene")
)


RNAseqgeneplot=function(gene,dataname){
    DATA=RNAseqAllDat[[dataname]]
    
    
    
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
  
  if(dataname=="TumorVsStroma"){
      difflabel="tumor versus stroma"
  }else{
      difflabel="basal versus classical"
  }
  
  beeswarm(as.numeric(DATA$exp[selectgeneEnsID,names(DATA$cl)])~DATA$cl,
  bty="l",pch=16,
  xlab=NA,
  ylab=ylabel,
  main=selectgenesym)
  
  
  boxplot(as.numeric(DATA$exp[selectgeneEnsID,names(DATA$cl)])~DATA$cl,add=T,axes=F,outline=F,
  col=DATA$col,
  sub=paste(difflabel,"log fold-change:",signif(DATA$diff[selectgeneEnsID,"logFC"],3),
  "\nFDR:",signif(DATA$diff[selectgeneEnsID,"adj.P.Val"],3)))
  
  beeswarm(as.numeric(DATA$exp[selectgeneEnsID,names(DATA$cl)])~DATA$cl,
  pch=16,
  add=T,axes=F,
  xlab=NA,
  ylab=NA,
  main=NA,cex=1.5)
  
}


selectCpG=function(gene,type,group){
    
    if(gene==""| gene ==" " | is.null(gene)){
        return(NULL)
    }else if(gene %in% methdiff[,"symbol"]){
        selectgeneCPGID=rownames(methdiff)[which(methdiff[,"symbol"]==gene)]
    }else if(tolower(gene) %in% rownames(methdiff)){
        return(tolower(gene))
    }else{
        return(NULL)
        
    }
    
    
    if(type!="any"){
        
        selectgeneCPGID=selectgeneCPGID[which(methdiff[selectgeneCPGID,"Type"]==type)]
        
    }
    
    if(group!="any"){
        

        selectgeneCPGID=selectgeneCPGID[which(methdiff[selectgeneCPGID,"GROUP"]==group)]
    }
    
    if(length(selectgeneCPGID)==0){
        return(NULL)
    }else if(length(selectgeneCPGID)==1){
            return(selectgeneCPGID)
    }else{
        return(selectgeneCPGID[order(methdiff[selectgeneCPGID,"adj.P.Val"])])
    }
    
    
    
    
}

Methgeneplot=function(cpgs,gene,type,group){
    
    
    if(is.null(gene) |is.null(cpgs)){
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("No CpG for",gene,"found with island type",
        type,"and with",group,"position"),cex = 1.6, col = "black")
        return(NULL)
        
        
    }
    
        beeswarm(as.numeric(methdat[cpgs[1],names(noutlierclassif)])~noutlierclassif,
    bty="l",pch=16,
    xlab=NA,
    ylab=paste("Methylation Beta-value (",cpgs[1],")"),
    main=paste(gene,"(",type,",",group,")"))
    
    
    boxplot(as.numeric(methdat[cpgs[1],names(noutlierclassif)])~noutlierclassif,
    add=T,axes=F,outline=F,
    col=momcolors,
    sub=paste( "Basal versus classical FDR:",signif(methdiff[cpgs[1],"adj.P.Val"],3)))
    
    beeswarm(as.numeric(methdat[cpgs[1],names(noutlierclassif)])~noutlierclassif,
    pch=16,
    add=T,axes=F,
    xlab=NA,
    ylab=NA,
    main=NA,cex=1.5)
    
}




shinyServer(function(input, output) {
    
    
    
    output$rnaseq_plot <- renderPlot({
        RNAseqgeneplot( gene=toupper(input$gene),input$data)
        
    })
    
    
    selectcpg=reactive({
        selectCpG(toupper(input$methgene),input$island,input$group)

    })
    
    output$methylation_plot <- renderPlot({
        
        Methgeneplot(selectcpg()
                ,toupper(input$methgene),input$island,input$group)
        
    })
    
    output$cpglist <- renderTable({
        df=methdiff[selectcpg(),c("symbol","Type","GROUP")]
        colnames(df)=c("Gene symbol","Island status","gene position")
        df=data.frame(CpG=rownames(df),df)
        df
    })

})
